# cg_md

`cg_md` is a C++17 command-line tool that orchestrates a full Martini 3
coarse-grained molecular dynamics pipeline for studying desmin protein
aggregation: coarse-graining an all-atom PDB with `martinize2`, clustering
copies with Packmol, solvating/ionizing the system, running energy
minimization and NVT/NPT equilibration with GROMACS, and finally driving
production MD either for a fixed length, adaptively (stopping once
aggregation-relevant sampling metrics converge), or under PLUMED metadynamics
biased by a pre-trained `mlcolvar` collective variable.

`cg_md` does not link against GROMACS, PLUMED, or any MD library. It is a thin
process orchestrator: it builds command lines for `gmx`, `martinize2`,
`insane`, `packmol`, and `plumed`, and shells out to them. That also means the
core logic (string/file/topology/index parsing) can be unit tested with
nothing but a C++ compiler — see [Testing](#testing).

## Running a batch

To run a whole manifest of simulations rather than one command:

- `tools/batch_host.sh` launches the container and checks the bind mount and
  the SIMD level before committing to the run;
- `tools/preflight.sh <manifest.tsv>` validates every entry — including the
  coarse-graining and bead count — before any compute is spent;
- `tools/run_batch.sh <manifest.tsv>` sequences the entries and survives their
  failures.

Each of those scripts documents its own procedure in its header. What the
2026-08-03 batch failure was, and the two rules that came out of it, are in
[archivio/_incident-2026-08-03/README.md](archivio/_incident-2026-08-03/README.md)
and `tools/lib/logdir.sh`.

## Building

```sh
cmake -S . -B build
cmake --build build
```

This produces two executables under `build/`:

- `cg_md` — the pipeline tool itself.
- `cg_md_tests` — the unit test suite (see [Testing](#testing)).

## Running

```sh
build/cg_md --project-dir /path/to/project --n-prot 10 --protomer desmin_head --stage all
```

Run `cg_md --help` for the full option list. Every path option (input PDB,
force fields, mapping directories, output directories) is resolved relative
to `--project-dir` unless given as an absolute path. `--stage` selects which
part of the pipeline to (re-)run: `all`, `cg` (coarse-grain one protomer and
stop — see below), `prepare`, `em`, `nvt`, `npt`,
`production`, `adaptive`, `metadynamics` (alias `metad`), `fes` (re-run just
the metadynamics finalization), `center` (re-center an existing trajectory),
or `report` (re-render the offline figures). `--dry-run` prints every command
it would execute without running anything or requiring the project's input
files to already exist.

## Restart continuity (adaptive chunks and metadynamics batches)

Both long-running production modes are split into segments so a
multi-microsecond run can be spread over several scheduler jobs. Re-invoking
the stage **resumes**; it never restarts. Continuity has three independent
parts, and the pipeline is responsible for all three:

| What carries over | How | Where |
|---|---|---|
| System state — coordinates, velocities, box, thermostat/barostat integrals, RNG state | `gmx grompp -c <prev>.gro -t <prev>.cpt` | `AdaptiveSampler::run()`, `MetadynamicsRunner::prepare_batch()` |
| Simulation clock | `tinit` in the segment's own `.mdp`, advanced by one segment length | `write_md_chunk_mdp()` / `write_metad_batch_mdp()` |
| Accumulated metadynamics bias | `RESTART` in the segment's `plumed.dat`, pointing at the same `HILLS` file | `write_metad_plumed_dat()` |

`grompp -t` restores the full physical state but always restarts GROMACS's
step counter at zero, which is why the advancing `tinit` is not optional: it
is the only thing that keeps the concatenated COLVAR time axis continuous
instead of having every segment re-emit `0 .. segment_length`.

Each segment writes its own numbered `.mdp` (`system/md_chunk_007.mdp`,
`system/md_metad_002.mdp`) rather than overwriting a shared one, so the exact
input that produced a given `.tpr` stays on disk.

A segment counts as complete only if its log says `Finished mdrun` **and** it
left behind the `.gro` and `.cpt` the next segment needs; the resume scan
stops at the first incomplete segment, so a hole in the middle invalidates
everything after it rather than silently splicing two trajectory pieces that
were never physically connected.

The segmentation parameters are recorded in a `*.ledger` file on the first
invocation (`adaptive.ledger`, `metadynamics/run.ledger`) and re-checked on
every later one. Resuming with a different `--adaptive-chunk-us` /
`--metad-chunk-us` / `--md-dt-ps` / `--metad-walkers` / `--n-prot` is refused
rather than accepted: `tinit` is computed as `segment × segment_length`, so a
changed length would keep the physical state from the last checkpoint while
stamping it with a time that was never reached — inventing or discarding a
slab of trajectory with no error anywhere downstream.

If a segment dies part-way, whatever it already appended to `HILLS`,
`COLVAR` and `COLVAR_monitor` is discarded before it is re-run. Otherwise the
re-simulated stretch would deposit its hills twice, and `plumed sum_hills` has
no duplicate filter. The rewind is by **record count** — how long each file was
when the previous segment finished, recorded in `plumed_records.ledger` — and
not by timestamp, for the reason in the next paragraph. Runs older than that
ledger fall back to counting clock restarts in the file, which is exact too.

**PLUMED's `time` column is not GROMACS's.** PLUMED derives time from its own
step counter, which `mdrun` resets on every invocation and which never sees the
`.mdp`'s `tinit`, so each segment's COLVAR starts again at 0 while the matching
`.edr` carries absolute times. Measured on `5x1-108_1`: `md_chunk_001.edr` runs
100 000–200 000 ps and `COLVAR_chunk_001.dat` runs 0–100 000 ps. The physics is
continuous; only the label restarts.

Every reader therefore puts the segments back on one increasing axis — a
timestamp that jumps *backwards* by more than half a frame starts a new segment
and is shifted to continue the previous one, while a repeated boundary frame is
dropped as the duplicate it is. `ColvarTable::read_colvars()` (C++),
`report_io.Clock` (the offline report) and `cvgen_estimators.rebase_clocks()`
(CV training) implement the same rule; the number of frames dropped as
duplicates is reported as `duplicate_frames_dropped`.

```sh
# 10 us of metadynamics as 20 x 0.5 us batches; re-run the same command
# after each job and it picks up where it left off, bias included.
build/cg_md --project-dir . --stage metad --metadynamics \
    --metad-total-us 10.0 --metad-chunk-us 0.5 ...
```

## Adaptive sampling: the stop rule (`--adaptive-*`)

`--adaptive` runs production in chunks of `--adaptive-chunk-us` and, after
every chunk, evaluates the whole trajectory so far (`SamplingMonitor`,
`adaptive_sampling_metrics.jsonl`). The run stops early when every gate
below holds, and otherwise at `--adaptive-max-total-us`. Reaching the cap is
not a failure: for a construct that never assembles it is the sampling one
wants for the negative result.

None of the thresholds is part of the run identity, so a run can be resumed
with tighter or looser ones.

The gates are of two kinds, and the distinction matters more than any
individual value.

**Exploration gates** answer "has X been seen yet". Each is a cumulative
count or a min/max range - a non-decreasing function of simulated time - so
it can only ever say "not yet", never "enough". Their job is to stop a run
from being declared converged while it is still all monomers.

| Option | Measures | Default | `sims.tsv` |
|---|---|---|---|
| `--adaptive-min-total-us` | time floor | 0.5 | 0.2 |
| `--adaptive-min-lcc-unique` | distinct largest-oligomer sizes seen (3 = a trimer has formed, 4 = a tetramer) | 3 | 3 for `5x1-108`, 4 for `5x1-263` |
| `--adaptive-min-contact-patterns` | distinct pairwise contact graphs seen | 10 | 10 |
| `--adaptive-min-cn-transitions`, `--adaptive-min-bidirectional-events` | crossings between thirds of the smoothed `cn_total` range, and the smaller of the up/down counts | 10, 1 | 20, 10 |
| `--adaptive-min-cn-range`, `--adaptive-min-rg-range` | span of `cn_total` / of the COM radius of gyration | 20, 0.5 | 0, 0 (off) |

**Convergence gates** answer "has more sampling stopped changing the
answer". These can go from satisfied back to unsatisfied, and they are the
ones the literature on sampling quality actually uses.

| Option | Measures | Default | `sims.tsv` | Basis |
|---|---|---|---|---|
| `--adaptive-max-pattern-growth` | contact patterns over the whole run divided by those over its first half | 1.10 | 1.10 | the cumulative-cluster-count plateau (Daura et al. 1999; Smith, Daura & van Gunsteren 2002) |
| `--adaptive-max-pattern-jsd` | Jensen-Shannon divergence, in bits, between the pattern frequencies of the two halves | 0.05 | 0.02 | the cluster-population comparison of Grossfield & Zuckerman 2009; 0 = same distribution, 1 = disjoint. For scale: a dominant pattern moving from 50% to 66% of the frames between halves is 0.02 bits, to 75% is 0.05 |
| `--adaptive-min-effective-samples` | independent samples of `cn_total` from the block-averaging plateau, N_eff = var / SE_max^2 | 20 | 50 | Flyvbjerg & Petersen 1989 for the estimator; Grossfield & Zuckerman: fewer than ~20 independent samples is unreliable |
| `--adaptive-min-assembly-events`, `--adaptive-event-residence-ps` | growth AND shrink events of the largest oligomer, each way, counting only changes that persist for the residence time | 20, 100 ps | 30, 100 ps | Poisson counting: the relative error on a rate from N events is 1/sqrt(N) - 30 gives 18%; the residence separates a collision from an association |
| `--adaptive-min-time-over-its` | run length over the slowest implied timescale of the descriptor set, which must also sit on a plateau against the lag (`scripts/cv_readiness.py`, cg_cvgen's own scan) | 10 | 10 | Sinitskiy & Pande 2018: the slowest timescale an analysis can report is bounded by the sampling, so an ITS of the order of the run length is that bound, not a measurement. This is the gate that says the run can train a CV: cg_cvgen picks its lag from the same plateau. 0 turns it off |

Why the growth ratio is not enough on its own: Caflisch's peptide-folding
work found that the total number of clusters never plateaus - rarely visited
ones keep turning up one at a time - while the populated ones converge early.
The divergence is insensitive to that tail, and `populated_contact_patterns`
(patterns holding at least 1% of the frames) is written to the metrics file
so the two can be read together.

Why the event count gates and the `cn_total` transition count does not: the
three `cn_total` states are thirds of a range that itself grows with the run,
so a "transition" there is not a physical event. Growth and shrinkage of the
largest oligomer are, and a run that assembled once and never came apart has
zero of the latter however many of the former it shows.

Why the tetramer gate differs between constructs: `1-108` is the head domain
alone, which has no reason to tetramerise, and requiring a tetramer there
would only ever postpone the stop to the cap; `1-263` carries coil 1B, the
A11 tetramer interface, and a run of it that has not made a tetramer has not
seen the biology it is there for.

Why the timescale gate matters more here than in a study that stops at the
unbiased run: these trajectories are the training data for the DeepTICA CV
that drives the metadynamics. The first real 207 ns of `5x1-108_1` gave a
slowest implied timescale of 110-126 ns - the run length - and a CV trained
on it would have learnt the run length. After every chunk
`scripts/cv_readiness.py` writes `cv_readiness.txt` in the run directory
with the timescale, the plateau and the ratio, so the number is on record
before anything is trained on it.

The script borrows the scan from `cg_cvgen/scripts/cvgen_estimators.py`
(numpy only) and looks for it in `$CG_CVGEN_SCRIPTS`, then next to this
project (`../cg_cvgen/scripts`, the layout of the repository, the Leonardo
build and the container), then `$CGMD_SRC/cg_cvgen/scripts`. Where none of
those exists - a Docker session that bind-mounts `cg_md` alone - the gate
cannot pass, the log says so after every chunk, and the run goes to its cap;
set the variable or pass `--adaptive-min-time-over-its 0`.

### Which numbers are published and which are ours

Only one threshold above is a number taken as such from a paper. The rest
rest on a published *principle* and a value chosen here; the distinction is
stated so that nobody quotes a choice as a citation.

| threshold | what the literature says, verbatim | status of the number |
|---|---|---|
| 20 independent samples (default; 50 in `sims.tsv`) | Grossfield & Zuckerman 2009: "As a conceptual rule-of-thumb, any estimate for the average of an observable which is found to be based on fewer than ~20 statistically independent configurations (or trajectory segments) should be considered unreliable." | **published**; 50 is that floor times 2.5 |
| block averaging as the estimator | Flyvbjerg & Petersen 1989, "Error estimates on averages of correlated data" | published method; no threshold is stated there |
| run length >= 10 x slowest implied timescale | Sinitskiy & Pande 2018: "A rule of thumb claiming that the slowest implicit time scale captured by an MSM should be comparable by the order of magnitude to the aggregate duration of all MD trajectories used to build this MSM has been known in the field." | published principle (an ITS of the order of the run length is the bound); the factor 10 is one order of magnitude of margin, **our choice** |
| pattern growth ratio <= 1.10 | Smith, Daura & van Gunsteren 2002: "a cluster analysis of the simulation trajectories is identified as a very effective method for judging the convergence of the simulations"; Caflisch's peptide work: the total cluster count does not plateau within 12.6 us while the populated clusters converge within 2 us | published method; 10% is **our choice**, and the JSD gate exists because the count alone is known to flatter |
| JSD between halves <= 0.02 bits | the segment-comparison test of Grossfield & Zuckerman; JSD as the ensemble-comparison statistic of ENCORE (Tiberti et al. 2015) and PENSA (Vögele et al. 2025) | published statistic; **no paper fixes a threshold**; 0.02 bits is a dominant pattern moving from 50% to 66% between halves |
| 30 assembly events each way, 100 ps residence | Poisson counting: the relative error on a rate estimated from N events is 1/sqrt(N); Pan et al. 2019 (PNAS) and the Strodel protocol (2022) count association and dissociation events but state no minimum | the 1/sqrt(N) is arithmetic; 30 (18%) and 100 ps (~400 ps of effective time at Martini's factor of 4) are **our choices** |
| contact: `cn_i_j >= 1` with r0 = 0.7 nm | Strodel protocol: "two proteins were considered to be in contact with each other if the minimum distance with respect to any two atoms from either protein was below 0.4 nm" (all-atom); Martini studies use 0.6-1.2 nm between beads | the coarse-grained analogue of a published definition; 0.7 nm sits inside the range used in the Martini literature |
| 0.2 us floor, 10 us cap | Strodel protocol: "MSMs require converged MD data, which usually implies tens of microseconds of MD sampling" (all-atom); Martini dynamics run about four times faster than atomistic ("a standard conversion factor of 4", Martini FAQ) | consistent with the published scale (10 us of Martini is ~40 us effective); the cap itself is **a budget** |

### Where the previous values came from

The values in `sims.tsv` before this rule (`--adaptive-min-cn-range 300
--adaptive-min-contact-patterns 120 --adaptive-min-cn-transitions 200
--adaptive-min-bidirectional-events 100 --adaptive-min-rg-range 2.0`) were
set from the first real 0.1 us chunk of `5x1-108_1`, at roughly two to twenty
times what that chunk showed (`cn_range` 154, 10 transitions, 2 patterns, 2
sizes), after the original thresholds turned out to be met 13-51x over by
that same chunk. They were guards against stopping at the time floor, not
convergence criteria, and two of them do not transfer between constructs:
`cn_total` scales with the square of the beads per interface, so 300 is two
simultaneous interfaces for the 108-residue head and a single encounter for
the 526-residue dimer; and the 2.0 nm `rg` range was read off `rg_global`,
which at the time jumped by 9.9 nm whenever a chain crossed a periodic
boundary. Both are now off (0) and the tetramer gate carries their meaning.

### What is deliberately not a stop criterion

Agreement between independent replicas - the strongest test there is
(Grossfield & Zuckerman) - cannot be evaluated from inside one run. The five
seeds of each construct exist for it; `analyze_offline.py aggregation` on
each and a comparison of their size distributions is the check to make
before quoting a number.

## Reproducible replicas (`--seed`)

For statistically independent repeats of the same setup (e.g. N metadynamics
replicas, or a handful of independent adaptive runs to pool COLVAR data for
CV training), pass a different `--seed` per replica:

```sh
build/cg_md --project-dir /path/to/project --n-prot 5 --protomer desmin_1-263 --seed 1 ...
build/cg_md --project-dir /path/to/project --n-prot 5 --protomer desmin_1-263 --seed 2 ...
```

`--seed` controls every stochastic step in the pipeline in one place:
Packmol's initial spatial packing (`seed` keyword - different starting
geometries), GROMACS's initial Maxwell-Boltzmann velocity draw (`gen_seed`,
NVT only), and the V-rescale thermostat's stochastic term (`ld-seed`, every
stage with a thermostat). The default, `-1`, matches Packmol/GROMACS's own
"pick a fresh pseudo-random seed" convention - i.e. today's non-reproducible
behaviour, unchanged unless you opt in.

**How the seed is spread across segments.** GROMACS does not keep a stateful
RNG stream that a checkpoint could restore: since 5.0 its stochastic
integration uses the counter-based ThreeFry2x64 engine, where a random number
is a pure function of `(ld-seed, step index, atom index, domain)`. Because
`grompp -t` restarts the step counter at zero for every segment, two segments
sharing an `ld-seed` would see the same counters and replay *bit-identical*
thermostat noise.

`cg_md` therefore writes a distinct `ld-seed` per segment — NVT gets
`seed + 0`, NPT `seed + 1`, fixed-length production `seed + 2`, adaptive
chunk *k* `seed + 100 + k`, metadynamics batch *k* `seed + 100000 + k`, and
walker *w* of that batch `+ 10000 w` on top - every walker starts from the
same `npt.cpt`, so with one shared `ld-seed` four walkers would replay one
trajectory four times. A replica is still exactly reproducible from its
`--seed`, and consecutive segments are still statistically independent. With the default `--seed -1`
grompp draws a fresh seed for every segment anyway.

(Note also that `grompp -t` carries over coordinates, velocities and box, but
*not* Nosé-Hoover / Parrinello-Rahman coupling variables — those need
`grompp -e`. This pipeline uses V-rescale and C-rescale, which are stochastic
and memoryless, so no such state exists and the segment chaining is exact.)

For a normal-then-metadynamics workflow, the recommended pattern is: run the
adaptive/unbiased stage once (or a few times with different seeds, if its
sampling metrics are only marginally above the `--adaptive-min-*` thresholds
and you want more diverse training data), train the CVs once with `cg_cvgen`,
then run the metadynamics stage N times with N different `--seed` values,
reusing the same trained CV model/params each time.

## Phosphorylated residues (SEP/TPO/PTR) and `--phospho`

There is no published Martini 3 parameter set for phosphorylated amino acids
(as of this writing, martini3IDP/vermouth only ships parameters for
lipidation PTMs). Phosphoserine support here is a hand-derived, best-effort
parametrization — not independently validated against an atomistic reference
— built as a vermouth **modification** (not a whole new residue "Block"; see
[vermouth's tutorials 6](https://vermouth-martinize.readthedocs.io/en/stable/tutorials/6_adding_residues_links/)
and [7](https://vermouth-martinize.readthedocs.io/en/latest/tutorials/7_adding_modifications/index.html)
for the underlying mechanism). Treating SEP as "phosphorylated SER" rather
than a standalone residue means it automatically inherits martini3IDP's real,
secondary-structure-aware backbone bonded parameters, instead of requiring a
hand-written (and easy to get subtly wrong) set of backbone links.

Required files, under the directories `--martinize-ff-dir`/`--martinize-map-dir`
point at (`Config::martinize_ff_dir`/`martinize_map_dir`, defaulting to
`force_fields`/`mappings` relative to `--project-dir`):

- `force_fields/charmm/modification.ff` — input-side (all-atom) recognition:
  declares the `SER-phos` modification's extra atoms (P, O1, O2, O3) attached
  to SER's `OG`.
- `force_fields/martini3IDP/modification.ff` — output-side CG parameters: the
  `SER-PO4` modification replaces SER's neutral `TP1` side-chain bead with a
  charged `Q5n` bead (matching the bead *type* martini3IDP already uses for
  GLU's full-size anionic side chain — chosen over ASP's smaller `SQ5n`
  because the phosphate group is bulkier than a carboxylate), charge set to
  **-2** (the phosphate's dominant protonation state at physiological pH;
  Martini's `charge` column is independent of bead `type`, so this doesn't
  need a dedicated "double-charge" bead type).
- `mappings/SEP.mapping` — the all-atom-to-CG mapping between the two
  modifications above (a different file format/extension than a Block-style
  `.map` file - see the vermouth tutorial 7 link above).

**Important - `-ff-dir` layout.** vermouth always recognizes input-structure
residues against a force field literally named `charmm`, regardless of
`-ff`. `-ff-dir` must therefore be the *parent* directory holding `charmm/`
and `martini3IDP/` as sibling subdirectories - pointing it directly at
`.../martini3IDP` hides the `charmm/` sibling and SEP silently fails to be
recognized at all (see `Workflow::check_inputs()`, which checks for this
specific misconfiguration). Verify with:
`martinize2 -ff-dir <dir> -list-blocks | grep -i sep`.

**Important - input PDB residue naming.** vermouth's modification mechanism
requires the PDB residue to be named after its *unmodified parent* (`SER`,
not `SEP`) — residue recognition/repair matches by name against the base
force field first, and only then explains extra atoms via a modification; a
residue still named `SEP` is never recognized at all. `Preparer::coarse_grain()`
therefore writes a throwaway `<protomer>_phosphorepaired.pdb` copy (SEP→SER,
TPO→THR, PTR→TYR in the fixed-width resName column only) before invoking
martinize2 for `--phospho` runs — the user's original PDB is never modified.
It also passes `-bonds-from both`, since name-based bond guessing alone knows
nothing about the nonstandard phosphate atoms.

If your PDB contains phosphoserine (residue name `SEP`), phosphothreonine
(`TPO`), or phosphotyrosine (`PTR`), run with `--phospho`:

```sh
build/cg_md --project-dir . --protomer desmin_head_pS27 --phospho ...
```

`Workflow::check_inputs()` fails fast, before any external tool runs, if:

- the input PDB is missing;
- the PDB contains a `SEP`/`TPO`/`PTR` residue but `--phospho` was not given;
- `--phospho` was given but the force-field/mapping directories above cannot
  be found, or `--martinize-ff-dir` has no `charmm/` subdirectory (see above).

Only SEP has a matching `modification.ff`/`.mapping` pair. To add TPO or PTR:
follow the same pattern (a `*-phos`/`*-PO4`-style modification pair plus a
`.mapping` file), reusing THR's/TYR's own existing side-chain bead as the
unmodified anchor the same way SER-PO4 does.

## Run diagnostics (`diagnostics.jsonl`)

Every dynamics stage with a thermostat (NVT, NPT, production, each adaptive
chunk, metadynamics) appends one JSON-Lines record to
`<result_dir>/diagnostics.jsonl` (or `<result_dir>/metadynamics/diagnostics.jsonl`
for metadynamics) after it completes - see `RunDiagnostics.hpp`. This is
**purely observational**: it never influences simulation control flow or the
adaptive stop criteria (that's `SamplingMonitor`/`adaptive_sampling_metrics.jsonl`,
a deliberately separate concern), and a diagnostics failure is caught and
logged as a warning rather than aborting the run.

Two categories are reported:

- **Energy/health** (from the stage's own `.edr` via `gmx energy`): mean,
  standard deviation, and linear drift-per-ns of temperature (every stage),
  plus pressure/density (NPT, production, adaptive chunks, metadynamics -
  skipped for NVT, which has no barostat) and total energy. A large
  fluctuation in pressure/density is normal for a small system at
  equilibrium; a nonzero **drift** in any of these is the actual red flag
  (thermostat/barostat miscoupling, timestep too large, an unstable
  interaction).
- **Structural** (from the stage's own COLVAR, when one exists - production,
  adaptive chunks, and metadynamics via an auxiliary monitoring-only
  `COLVAR_monitor` PRINT that does not feed the CV/bias): per-protomer and
  global radius-of-gyration mean/stddev/drift, total contact number, the size
  of the largest oligomer and the number of oligomers per frame, a histogram
  of every cluster size observed and of the per-frame largest cluster (which
  sums to the frame count, so it is directly a distribution over aggregation
  states), and the per-pair contact occupancy - the inter-protomer contact
  map, which says *which* interfaces dominate rather than only how many
  contacts existed.

  Time axes come from PLUMED's own `time` column, so a segment starting at
  t = 100 ns reports drift against real simulated time.

All energy observables for a stage are extracted in a **single** `gmx energy`
invocation rather than one per quantity.

## Non-bonded parameters (`epsilon_r`)

Every `.mdp` — energy minimization included — is written from one
`nonbonded_block()`, so no stage can end up with different electrostatics
from another. The block sets `epsilon_r` explicitly (`--epsilon-r`, 15 by
default).

That is not a stylistic choice. Martini's coarse-grained water bead cannot
screen charges the way explicit water does, so the model compensates with a
uniform relative dielectric of 15 applied everywhere, and the whole force
field is parameterized against it. GROMACS's own default is `epsilon_r = 1`;
omitting the line makes every electrostatic interaction fifteen times too
strong, which shows up as over-stable salt bridges, excessive ion
condensation and grossly overestimated charge-driven association — most
damagingly in exactly the comparison this pipeline exists for, phosphorylated
versus unphosphorylated protomers. Use `--epsilon-r 2.5` only with
polarizable Martini water, which does its own explicit screening.

`epsilon_rf = 0` (infinite dielectric beyond the cut-off), the cut-offs, and
`coulomb-modifier` / `vdw-modifier` are likewise written out rather than left
to defaults, because Martini treats the cut-off treatment as part of the
parameterization.

## Metadynamics and the CV contract

The DeepTICA CV that METAD biases is trained (by `cg_cvgen`) on the columns of
an *unbiased* production COLVAR. The biased run therefore has to hand the
network the same descriptors, computed the same way, in the same order — any
divergence either aborts PLUMED with a size mismatch or, worse, silently
evaluates the model on inputs it has never seen.

That is enforced structurally rather than by convention:

- **One descriptor builder.** `build_descriptors()` emits the groups,
  `WHOLEMOLECULES`, per-protomer COM, pairwise COM distances, pairwise
  `COORDINATION` (with `NN`/`MM` written explicitly rather than inherited
  from a PLUMED default), `cn_total`, per-chain `rg1..rgN` and `rg_com`.
  Both `write_plumed_dat()` and `write_metad_plumed_dat()` use it.
- **One canonical order.** `d_i_j…`, `cn_i_j…`, `rg1..rgN`, `rg_com`,
  `cn_total` — exactly what the unbiased COLVAR PRINTs and exactly what
  `cg_cvgen`'s default `--feature-regex` selects.
- **The model states its own inputs.** `MetadynamicsRunner` reads
  `feature_cols`, `sigma`, `grid_min`/`grid_max` and `n_cvs` from the trained
  model's `cv_params.pkl` in one pass and cross-checks them against
  `--metad-nodes`. `write_metad_plumed_dat()` refuses to emit a file naming a
  descriptor the block does not define.

`WHOLEMOLECULES` matters more than it looks: PLUMED receives raw periodic
coordinates, so a protomer straddling a box face is reconstructed as two
half-molecules on opposite sides — its COM lands mid-box and its Rg is
roughly half a box vector. Every distance, contact and gyration built on that
(and hence the bias itself) would be meaningless.

Two further defaults exist for reasons worth stating:

- **The bias lives on a grid** (`GRID_MIN`/`GRID_MAX`/`GRID_BIN`, bounds taken
  from the trained CV's sampled range). Gridless METAD re-evaluates every
  Gaussian ever deposited at every MD step, so its cost grows linearly with
  simulated time: at `PACE=500` and `dt = 0.01 ps` that is ~2·10⁵ hills per
  microsecond, and a multi-microsecond run ends up spending most of its wall
  time inside PLUMED. `--no-metad-grid` opts out.
- **`CALC_RCT` is on**, producing `metad.rbias`. This is the only route from a
  biased trajectory to a free energy along an observable that was *not*
  biased — oligomer size, contact number, radius of gyration — and it requires
  a grid, which is why `--no-metad-grid` without `--no-metad-rct` is rejected
  rather than silently dropping the capability.

With multiple walkers, each walker writes only the hills *it* deposited to its
own `HILLS` (the sharing happens in memory over MPI), so `finalize()` passes
every walker's file to `plumed sum_hills`. Summing one of them would
reconstruct 1/N of the bias and understate every barrier.

## Permutation invariance of the CV features

Identical protomers are interchangeable, so a configuration and the same
configuration with two chains relabelled are **one physical state** — but they
produce completely different `d_i_j` vectors. A CV trained on the raw,
index-ordered vector therefore gives one state many different values, and a
metadynamics bias built on it is incoherent.

Both halves of the pipeline now remove that degeneracy the same way (the
approach of Samantray et al. / TICAgg): each permutable block is sorted.

- `cg_cvgen` sorts `d_i_j`, `cn_i_j` and `rg_i` row-wise before training, and
  records `permutation_invariant` in `cv_params.pkl`.
- `cg_md` emits matching PLUMED `SORT` actions and feeds `PYTORCH_MODEL` the
  components `sorted_d.1…`, `sorted_cn.1…`, `sorted_rg.1…`, `rg_com`,
  `cn_total`, in exactly that order.
- `MetadynamicsRunner` reads the flag from the model's own parameter file, so
  a CV trained on raw features cannot be silently driven with sorted ones.

Blocks are sorted *independently* rather than co-sorted by a shared key. That
loses the pairing between a distance and its own contact number, but it is
what PLUMED's `SORT` can reproduce at biasing time — and a feature vector
`cg_md` cannot rebuild exactly would be worse than one carrying slightly less
information. `--no-permutation-invariant` restores the old behaviour, which is
only correct for a single protomer pair.

The COLVAR itself stays **unsorted**: it is the raw record, and
`RunDiagnostics` and the report need the per-pair columns.

## Secondary structure (`--dssp`)

The default is an all-coil `-ss` string. That is right for an intrinsically
disordered fragment and wrong for anything containing a coiled coil or folded
domain, whose helices would get coil bonded parameters and be held together
only by the elastic network. `--dssp` makes martinize2 derive the assignment
from the input structure instead.

Either way, note what Martini cannot do here: the assignment is **fixed at
setup**, so secondary structure cannot change during the simulation. β-sheet
formation — the central observable of all-atom amyloid studies — is not
observable in this model. What is observable is association, its strength, and
which interfaces form.

## A relaxed starting structure (`--stage relax`)

Packing copies of the raw input means every copy starts in the same, possibly
strained conformation, and they aggregate while still relaxing — so part of
what assembles is an artefact of the starting geometry. `--stage relax`
simulates one protomer in solvent for `--relax-us`, clusters the trajectory
(`gmx cluster`, gromos), and installs the centroid of the most populated
cluster as the structure Packmol replicates — written to its own
`<protomer>_cg_relaxed.pdb`, so the raw coarse-grained `<protomer>_cg.pdb` it
came from is left untouched and a later re-run of the coarse-graining cannot
silently throw the relaxation away.

Run it between `prepare` and the packing:

```sh
build/cg_md ... --stage prepare      # coarse-grain
build/cg_md ... --stage relax        # relax + cluster + install representative
build/cg_md ... --stage all          # pack, solvate, equilibrate, produce
```

## Residue-level interface map (`--stage contact-map`)

`scripts/contact_map.py` reads the *trajectory* (mdtraj) rather than the
COLVAR, because two quantities cannot be had from the per-protomer
descriptors:

- **Inter-chain residue-residue contact probability** — the per-protomer
  occupancy says *how much* two chains touch; this says *where*. For a
  phosphorylation study that is the difference between "the variants associate
  differently" and "the variants use a different interface".
- **NPMI shape** — ratio of the smallest to the largest principal moment of
  inertia, the descriptor used in transition-network analyses: 1 for a compact
  globule, small for an elongated aggregate.

Adaptive runs now concatenate their chunk trajectories into
`md_all_center.xtc` (whole molecules, centred) at the end of the run, so this
analysis works on them exactly as on a fixed-length production run.

## What counts as a "protomer" (`--merge-chains`)

Everything in this pipeline — the per-protomer index groups, the inter-protomer
distances and contacts, `--n-prot`, `--atoms-per-prot` — is defined per **CG
molecule**. Choosing what that molecule is is a modelling decision, not a
formality.

For an intermediate-filament fragment the physical unit is the parallel
coiled-coil **dimer**, not a single chain: coil1A/coil1B exist only as a
two-chain coiled coil, so an isolated chain has no partner and its helix would
be held together purely by the artificial elastic network. Feed the dimeric PDB
and merge the chains:

```sh
build/cg_md ... --merge-chains A,B --seq-length <residues in BOTH chains>
```

`--n-prot` then counts dimers, and every descriptor becomes inter-dimer — which
is the association that forms a tetramer, the first real step of assembly.

A multi-chain PDB **without** `--merge-chains` is rejected up front: the
topology patching includes exactly one `.itp` and writes one `[ molecules ]`
entry, so unmerged chains would leave an undefined moleculetype and fail deep
inside `grompp` with a message that says nothing about the cause.

The intrinsically disordered head domain is the exception — it is not helical
and does not coiled-coil, so a single-chain construct is a legitimate model of
its own self-association propensity.

## Dimerization PMF (`--stage pmf`)

A free energy read off an N-chain box is not a binding free energy: it carries
a translational-entropy term set by the box concentration, so it is only
comparable between systems simulated at *identical* concentration, and the
statistics rest on however many association events N chains happen to produce.

`--stage pmf` runs the same well-tempered metadynamics machinery — batching,
restart, grid, `CALC_RCT`, `sum_hills`, FES convergence, crash recovery — but
biases the inter-chain centre-of-mass distance `d_1_2` of a **two-chain**
system, with a flat-bottom `UPPER_WALLS` confining the pair to a finite
sampling volume. One pair associating and dissociating repeatedly gives far
better statistics, and the result is corrected to the standard state, so it is
comparable across variants outright.

```sh
build/cg_md --project-dir . --stage pmf --n-prot 2 --protomer desmin_head \
    --pmf-wall-nm 6.0 --metad-total-us 10.0 --metad-chunk-us 0.5
```

`--pmf` implies `--metadynamics`, sets one CV, and derives sigma/grid from the
`--pmf-*` options — no trained model or `cv_params.pkl` is involved. Output
lands in `runs/<system>/pmf/`.

### The correction, and why it is not optional

`plumed sum_hills` returns `F(r) = -kT ln P(r)` where `P(r)` is the probability
*density in r*. That density already contains the `4πr²` phase-space factor, so
`F(r)` is **not** the pair potential of mean force — a pair with no interaction
at all still shows `F(r)` falling as `-2kT ln r`. Reading a binding free energy
off the depth of `F(r)` therefore mixes in the translational entropy.

The pair potential is recovered by removing the Jacobian and referencing it to
the unbound plateau:

```
w(r) = F(r) + 2kT ln r  −  ⟨F + 2kT ln r⟩_plateau
Ka   = 4π ∫₀^rb r² exp(−w(r)/kT) dr        [nm³]
ΔG°  = −kT ln(Ka / V°),   V° = 1.66054 nm³  (1 mol/L)
```

Integrating only over the bound region `r < rb` is what makes the answer
independent of the sampling volume the wall defines. The report emits `w(r)`
next to the raw `F(r)`, and ΔG° as a function of `rb` — a ΔG° that plateaus
with the cut-off is trustworthy; one that keeps sliding is not.

`analyze_run.py`'s implementation is validated against analytical square wells
and against the null case (a non-interacting pair confined to exactly the
standard-state volume must give ΔG° = 0); it reproduces both to ~0.002 kT.

### Limits worth knowing before you use it

- **The coordinate.** A centre-of-mass distance is sound for compact domains.
  For elongated (e.g. coiled-coil) constructs it conflates orientation with
  separation — side-by-side and end-to-end arrangements share the same `r` —
  and the PMF converges poorly and reads ambiguously.
- **The wall must fit the box.** Beyond half the shortest box vector the COM
  separation is aliased by the minimum-image convention and the unbound
  plateau is meaningless. `Config::validate()` warns when `--pmf-wall-nm` looks
  too large for the requested box, but it cannot know the real solvated box.
- **The force field.** Martini 3 over-stabilizes protein–protein interactions
  for flexible proteins in solution. Thomasen et al. (*Nat Commun* 15:6645,
  2024) correct it either by strengthening protein–water by ~10% or by
  weakening protein–protein to λ_PP = 0.88, and recommend the second; this
  pipeline uses the first (`--martini-lambda-pw`). Note also that martini3IDP
  already addresses IDP compaction through its bonded parameters and is
  presented as an alternative to λ rescaling rather than a base for it — see
  `tools/make_lambda_itp.py`. Absolute ΔG° will be too negative; differences
  between variants are more robust but not immune, since a phosphate changes
  the charge and the imbalance is not uniform.
- **Dimerization is not aggregation.** This ranks self-association propensity;
  it does not describe how many chains assemble or which higher-order
  architecture they form. The many-chain metadynamics stage answers that, via
  the per-pair contact map and the oligomer-size distribution.

## Offline report (`report/`)

`scripts/analyze_run.py` runs automatically at the end of production,
adaptive and metadynamics stages (`--no-report` to disable, `--stage report`
to re-render) and writes figures plus a `summary.json` into
`<run_dir>/report/`. Like `RunDiagnostics` it is purely observational: a
missing matplotlib or a parsing failure is a warning, never a run failure.

`analyze_run.py` is only the entry point; the work sits in sibling modules
it imports by name, all in `scripts/`: `plumed_clock.py` (COLVAR/HILLS/FES
files onto one clock), `report_io.py` (JSON-Lines, thinning, saving),
`report_metrics.py` and `pmf_metrics.py` (the numbers; pure numpy),
`report_health.py`, `report_structure.py`, `report_adaptive.py`,
`report_metad.py` and `report_pmf.py` (one figure family each).

Unbiased / adaptive runs get: per-stage run health (mean ± sd and drift of
temperature, pressure, density, total energy); `cn_total`, `rg_com`,
largest-oligomer size and oligomer count versus time; the oligomer size
distribution; the mean inter-protomer contact map; per-chain Rg traces and
distributions; the `-kT ln P(cn_total, rg_com)` landscape; and a
convergence panel with the block-averaging error curve (the plateau of that
curve, not the naive N^-1/2 value, is the honest error bar on a correlated
series; `summary.json` records its maximum and the block length it came
from) plus the running mean. Adaptive runs additionally get the stop-criteria values and
a pass/fail grid per chunk.

Metadynamics runs get: the well-tempered hill-height decay (which must fall
towards zero); the final FES; the FES convergence curve from
`sum_hills --stride` (`--metad-fes-stride N`); the CV and bias trajectories;
the same structural panels computed from `COLVAR_monitor`; and the
`c(t)`-reweighted free energy along `cn_total` and `rg_com`.

## Offline analysis of copied data (`analyze_offline.py`)

`analyze_run.py` above needs a run directory in the layout the pipeline
builds. `scripts/analyze_offline.py` takes file paths and globs instead, so it
works on whatever was rsynced off Leonardo, in whatever layout, on a laptop
with numpy and matplotlib and nothing else installed. It never writes into a
run, and it is not invoked by anything.

```sh
# how the chains assembled
python3 scripts/analyze_offline.py aggregation \
    --colvar 'copied/COLVAR_chunk_*.dat' --n-prot 5 --out-dir report

# whether the CV is sound, and whether it was used where it was trained
python3 scripts/analyze_offline.py cv \
    --cv-dir cvs/5x1-108 --colvar 'metad/walker*/COLVAR' \
    --monitor 'metad/walker*/COLVAR_monitor' --out-dir report

# whether the metadynamics converged and the reweighting is worth anything
python3 scripts/analyze_offline.py metad \
    --hills 'metad/walker*/HILLS' --colvar 'metad/walker*/COLVAR' \
    --monitor 'metad/walker*/COLVAR_monitor' --fes metad/fes.dat \
    --fes-convergence 'metad/fes_convergence/fes_*.dat' --out-dir report
```

The globs can be quoted (the script expands them, in natural order) or left
to the shell; the file list is the same either way.

### `aggregation` -> `aggregation_summary.json`, figures 40-42

| What | Why it is there |
|---|---|
| Size distribution, number- and mass-weighted | `<s>_n` counts oligomers, `<s>_w` counts chains. A box that is mostly free monomers with one aggregate has `<s>_n` near 1 and a large `<s>_w`; quoting only the first hides the aggregate. |
| `-kT ln (N_s / N_1)` | The same distribution on an energy axis, referenced to the free chain. Depends on the box and on N: it ranks sizes within a run, it is not a standard-state `dG`. |
| Contact lifetimes, continuous and intermittent | The continuous one ends at any break and is short by construction at a contact cut-off; the intermittent one tolerates recrossings. They answer different questions and the pair is the useful output. |
| First-passage time to each size | When a size was first reached, and `null` when it never was - which is not the same as "at the end". |
| `censoring.fraction_at_n_prot` | How often the largest oligomer was every chain in the box. The distribution cannot go past N, so a large value means the upper tail is measuring the box. |

### `cv` -> `cv_quality.json`, figures 50-52

| What | The failure it catches |
|---|---|
| Implied-timescale plateau, and whether the selected lag sat in it | The lag fell back to raw VAMP-2, which always favours the shortest one. |
| VAMP-2 train vs test, Chapman-Kolmogorov error, spectral gap | A CV that fits the training replicas and nothing else. |
| Training coverage | The biased run left the region the network was fitted on. A DeepTICA CV is a neural network: outside its training box the output is an extrapolation with no reason to be monotone or bounded, and metadynamics pushes towards exactly the rare regions the training set is thinnest in. |
| Grid coverage | Frames outside the METAD grid (PLUMED stops with "Extrapolating from the grid"), or a grid so oversized that every `sum_hills` bin is coarser than it needs to be. |
| Component independence | TICA components are uncorrelated by construction; a large off-diagonal means the exported network is not the fitted one, usually a feature vector handed to PLUMED in a different order. |
| Rank correlation against `cn_total` and `rg_com` | What the CV actually tracks. Joined on the time column, because `COLVAR` and `COLVAR_monitor` come from two PRINT actions with independent strides. |

### `compare` -> `replica_comparison.json`, figure 43

```sh
python3 scripts/analyze_offline.py compare \
    --summaries 'runs/5x1-108_*/report/aggregation_summary.json' --out-dir report
```

Pairwise Jensen-Shannon divergence between the size distributions of
independent replicas, and the pooled distribution with a between-replica
standard error per size. This is the test the stop rule cannot make and the
five seeds per construct exist for; `--max-jsd` (0.05) is the largest
pairwise divergence still reported as agreement, the same scale as the
within-run gate.

### `metad` -> `metad_quality.json`, figures 60-62

| What | The failure it catches |
|---|---|
| Hill-height decay per walker | Necessary but not sufficient: a bias can grow smoothly in a basin the walker never leaves. |
| Recrossings of the CV | The sufficient part. A free-energy difference between two regions is only measured if the trajectory went back and forth between them, and the round-trip count is the sample size for it. |
| FES RMSD to the final surface | Compared only where the surface is populated (`--fes-window`, default 30 kJ/mol) and after removing the mean offset. Empty grid otherwise moves by tens of kJ/mol between blocks for reasons that are not physics. |
| Kish effective sample size | Reweighting weights span orders of magnitude. When a handful of frames carry the weight, the reweighted histogram is those frames however many were recorded. |
| Reweighted `F(s)` with a block error band | The PLUMED block-analysis estimator over the weighted samples. Blocks are taken inside one walker at a time, so none straddles the join between two. |
| `c(t)` monotonicity | `c(t)` cannot decrease. A backward step means the rows are not in time order - a restart read without rebasing the clock, or stale rows from a batch that died. |
| Walker balance | A walker that deposited far fewer hills than the rest was stalled; one with none never ran, and the run had fewer walkers than it was charged for. |

The estimators are numpy only, no torch, no mdtraj: `aggregation_clusters.py`,
`aggregation_kinetics.py`, `aggregation_shape.py` (gathered by
`aggregation_metrics.py`), `cv_features.py`, `cv_stats.py`, `cv_quality.py`,
`metad_reweight.py`, `metad_quality.py`. The figures are in
`figures_aggregation.py`, `figures_cv.py`, `figures_metad.py`; the file
handling in `offline_io.py`; and there is one small driver per subcommand
(`offline_aggregation.py`, `offline_cv.py`, `offline_metad.py`,
`offline_compare.py`).
`tests/test_analysis.py` checks every estimator against arithmetic worked out
separately and then runs the three subcommands end to end.

## Architecture

| Component | File | Responsibility |
|---|---|---|
| `Config` | `Config.hpp`, `Config.cpp` + `ConfigOptions.cpp` (argv) + `ConfigHelp.cpp` (`--help`) + `ConfigValidate.cpp` (consistency checks) + `ConfigIdentity.cpp` (what makes a run the same run) | Parses argv, holds every tunable parameter, computes derived paths. Never touches the filesystem or spawns processes. |
| `Shell` | `Shell.hpp/.cpp` | Runs external commands (`std::system`/`popen`), honoring `--dry-run` and verbose logging. |
| `Workflow` | `Workflow.hpp/.cpp` + `WorkflowChecks.cpp` (directories, `residuetypes.dat`, `check_inputs()`) | Top-level orchestrator: directory bootstrapping, `check_inputs()` validation, `.mdp` generation, and stage dispatch to the collaborators below. |
| `GromacsDriver` | `GromacsDriver.hpp/.cpp` | Shared wrappers for the `gmx` subcommands used everywhere (`grompp`, `mdrun`, `energy`, `trjconv`). |
| `Preparer` | `Preparer.hpp`, `PreparerCoarseGrain.cpp` (martinize2) + `PreparerPackmol.cpp` + `Preparer.cpp` (solvation, ions, topology) + `PreparerChecks.cpp` (the `verify_*()` guards) | All-atom PDB → solvated, ionized, EM-ready CG system: `martinize2`, Packmol, box/solvation, ion placement, topology patching. |
| `SimulationRunner` | `SimulationRunner.hpp/.cpp` | Fixed-length legs: EM, NVT/NPT equilibration, fixed-length production, trajectory centering. |
| `AdaptiveSampler` | `AdaptiveSampler.hpp/.cpp` + `AdaptiveSamplerFinalize.cpp` (concatenate, centre, report) | Chunked production MD with early stopping once aggregation-sampling metrics converge (see `SamplingMonitor`). |
| `RelaxRunner` | `RelaxRunner.hpp`, `RelaxRunner.cpp` (the single-protomer run) + `RelaxRunnerSelect.cpp` (clustering, the representative and the distinct copies) | `--stage relax`: one protomer in solvent, clustered, the centroid installed as the starting conformation. |
| `MetadynamicsRunner` | `MetadynamicsRunner.hpp`, `MetadynamicsRunner.cpp` (the stage) + `MetadynamicsBatch.cpp` (one batch, restart) + `MetadynamicsGuards.cpp` (walker-sharing and restart guards) + `MetadynamicsCvParams.cpp` (reads `cv_params.pkl`) + `MetadynamicsFinalize.cpp` (`sum_hills`, centring, report) + `MetadynamicsInternal.hpp/.cpp` (file names and helpers the four share) | PLUMED metadynamics production driven by a trained `mlcolvar` CV, optionally multi-walker. |
| `MdpWriter` | `MdpWriter.hpp/.cpp` (the stage writers) + `MdpProduction.cpp` / `MdpSpec.hpp` (the segment block they all build on) | Renders GROMACS `.mdp` files from `Config`. |
| `PlumedWriter` | `PlumedWriter.hpp`, `PlumedDescriptors.cpp` (what is measured) + `PlumedWriter.cpp` (the unbiased `plumed.dat`) + `PlumedWriterMetad.cpp` (the metadynamics one: bias, walls, grid) | Renders `plumed.dat` (inter-protomer distances/contacts/Rg, optionally the metadynamics CV). |
| `SamplingMonitor` | `SamplingMonitor.hpp/.cpp` + `SamplingEstimators.hpp/.cpp` (the numerics on plain vectors) + `CvReadiness.hpp/.cpp` (the timescale gate's input) | Reads COLVAR files and computes the metrics `AdaptiveSampler` uses to decide when to stop. |
| `RunDiagnostics` | `RunDiagnostics.hpp/.cpp` | Purely observational health/structural metrics (see "Run diagnostics" above) - never affects control flow. |
| `ColvarTable` | `ColvarTable.hpp/.cpp` + `ColvarTruncate.cpp` (rewinding a file to the last completed batch) | Shared PLUMED COLVAR parser, used by `SamplingMonitor` and `RunDiagnostics`. |
| `ClusterAnalysis` | `ClusterAnalysis.hpp/.cpp` | Shared connected-component (oligomer) analysis over a pairwise-contact graph, used by `SamplingMonitor` and `RunDiagnostics`. |
| `Reporter` | `Reporter.hpp/.cpp` | Locates and invokes `scripts/analyze_run.py` for the offline figures/summary. Never affects control flow. |
| `TimeSeriesStats` | `TimeSeriesStats.hpp/.cpp` | Mean/stddev/linear-drift reduction of a scalar time series, and a minimal `.xvg` reader for GROMACS energy output. |
| `IndexBuilder` | `IndexBuilder.hpp/.cpp` | Builds/reads/writes GROMACS `.ndx` index files, including per-protomer groups. |
| `TopologyEditor` | `TopologyEditor.hpp/.cpp` | Patches a `martinize2`-generated `.top` file with Martini includes and the protomer count. |
| `GroUtils` | `GroUtils.hpp/.cpp` | Normalizes ion residue/atom names in `.gro` files to Martini's `NA+`/`CL-` spelling. |
| `FileUtils`, `StringUtils` | `*.hpp/.cpp` | Small, dependency-free file and string helpers used throughout. |

`Preparer`, `SimulationRunner`, `AdaptiveSampler`, and `MetadynamicsRunner` each
hold references to the single `Config`, `Shell`, and `GromacsDriver` instances
owned by `Workflow` — nothing is duplicated or copied.

## Testing

```sh
cmake --build build --target cg_md_tests
build/cg_md_tests        # or: ctest --test-dir build
```

`cg_md_tests` covers the pure-logic modules (one `tests/test_<module>.cpp`
per module: `SamplingMonitor`, `SamplingEstimators`, `CvReadiness`,
`ColvarTable`, `PlumedWriter`, `MdpWriter`, `Config` identity, and the rest)
with a small, dependency-free harness (`tests/test_framework.hpp`) and needs nothing beyond
the standard library — no GROMACS, `martinize2`, Packmol, or PLUMED
installation required. It does not exercise `Preparer`/`SimulationRunner`/
`AdaptiveSampler`/`MetadynamicsRunner`/`GromacsDriver`, since those only
assemble and shell out to external commands; use `cg_md --dry-run` to inspect
the exact command lines a given configuration would produce.

The offline analysis modules have a Python suite, which needs numpy and uses
matplotlib if it is there:

```sh
python3 tests/test_analysis.py
```

`test_analysis.py` runs the `tests/test_analysis_*.py` files in turn and adds
up their tallies; each of those runs on its own too. They cover the
aggregation, CV and metadynamics estimators case by case, then drive the
`analyze_offline.py` subcommands over COLVAR, HILLS and `sum_hills` files
written for the purpose - including one whose clock restarts half way
through, the shape a batched run really has - and `cv_readiness.py` on a
series with a known slowest timescale.

The batch tooling in `tools/` has its own suite:

```sh
./tests/test_tools.sh    # inside the container: bash 4+ needed
```

`test_tools.sh` sources the parts in `tests/tools/` in order, in one shell,
so they share the stub `cg_md`, the temporary directory and the tally. They
cover `tools/lib/{manifest,rundir,logdir,resources,entry}.sh` and drive
`tools/run_batch.sh` end to end against a stub `cg_md`, including the failure
that lost the 2026-08-03 batch: the log directory being deleted mid-run. See
`tools/lib/logdir.sh` for that rule and
[archivio/_incident-2026-08-03/README.md](archivio/_incident-2026-08-03/README.md)
for the preserved evidence.

Both suites, plus the build, run with one command from the Mac:

```sh
./tools/build.sh
```

### Validating a manifest before running it

`--stage cg` coarse-grains a single protomer and stops. That is enough to check
the input PDB exists, that its residue count matches `--seq-length`, that its
chains match `--merge-chains`, and that the bead count martinize2 produces
matches `--atoms-per-prot` — all in seconds, instead of after the relaxation,
Packmol, solvation and EM. `tools/preflight.sh` runs it for every distinct
protomer in a manifest, and `./tools/batch_host.sh --check` runs that.
