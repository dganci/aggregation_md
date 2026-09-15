# Visualizing trajectories in VMD

Martini 3, rhombic-dodecahedron box, 5 free protomers. A trajectory is only
usable for visualization if it satisfies **both** properties at once:

* **whole** — no molecule split across a cell face;
* **continuous** — no atom teleports between consecutive frames.

Checking only one of them is how you end up with a file that looks fixed and
isn't. Both checks are in section 5.

## The two rules

**1. `-pbc whole` must come first.** `-pbc nojump` moves each atom independently
to the image nearest its own previous position — it never looks at bonds. An
atom sitting on a cell face can be assigned to a different image than the rest
of its molecule, and once frame 0 is broken it stays broken for the whole run.
Only `-pbc whole` uses the bonded topology to reconstruct molecules, which is
why it needs a `.tpr`.

**2. `-pbc nojump` must come last.** `-center`, `-pbc mol` and `-ur compact`
wrap coordinates back into the cell, re-creating exactly the jumps `nojump`
removed. Running `nojump` then `center` does not combine the two effects — the
second call cancels the first.

Measured on `runs/5x1-108_1_ca3b5833/` (5 chains, 1000 frame transitions;
"broken" = consecutive backbone beads further apart than 1 nm):

| file | pbc order | max COM jump | jumps > 6 nm | max BB-BB | broken frames |
|---|---|---|---|---|---|
| `..._center.xtc`      | `mol,center` last | 26.99 nm | 351 | 0.43 nm | 0 |
| `..._center_mol.xtc`  | `mol` last        | 19.63 nm | 143 | 0.43 nm | 0 |
| `..._nojump_mol.xtc`  | `nojump` last, no `whole` |  0.61 nm | **0** | **19.61 nm** | **51/51** |
| `vis_prot_nojump.xtc` | `whole` then `nojump` | **0.61 nm** | **0** | **0.45 nm** | **0** |

The `_center` files are whole but jump. `_nojump_mol` is continuous but every
frame has a chain split in half — the failure mode rule 1 describes. Only the
last row is both.

## 0) Container

```bash
docker run --rm -it --platform linux/amd64 \
  --mount type=bind,source="<path/to/aggregation_md/cg_md>",target=/data \
  dganci/aggregation-md:alma9-metatomic bash

cd /data/runs/<RUN>
G=/opt/gromacs-plumed/bin/gmx_mpi
```

## 1) Concatenate the chunks FIRST

`nojump` re-references itself at the start of every input file, so running it per
chunk leaves a discontinuity at each chunk boundary.

```bash
$G trjcat -f md_chunk_000.xtc md_chunk_001.xtc -o all.xtc -cat
```

Skip only if you are visualizing a single chunk.

## Branch A — continuous trajectory (what you normally want)

Verified end to end on chunk 000. Three calls, in this order.

```bash
R=/data/runs/<RUN>

# A1  whole molecules (needs the .tpr for bonds) + strip solvent
echo 0 | $G trjconv -s md_chunk_000.tpr -f all.xtc -n index.ndx \
    -o whole_prot.xtc -pbc whole
    # 0 = AllProteins  ->  1180 beads, ~5 MB instead of ~280 MB

# A2  reference structure: frame 0 of the same pipeline
echo 0 | $G trjconv -s md_chunk_000.tpr -f all.xtc -n index.ndx \
    -dump 0 -pbc whole -o $R/vis_prot_ref.gro

# A3  remove the jumps — nothing after this
echo 0 | $G trjconv -s $R/vis_prot_ref.gro -f whole_prot.xtc \
    -o $R/vis_prot_nojump.xtc -pbc nojump
```

**A3 must take `-s vis_prot_ref.gro`, never the `.tpr`** — and this is a
correctness requirement, not just bookkeeping. The obvious reason is that after
A1 the trajectory holds 1180 atoms so the reference has to match; but the rule
holds even when the counts *do* match, i.e. when you keep the full system and
the `.tpr` would be accepted. Measured on chunk 000, largest gap between
consecutive backbone beads:

| reference | frame spacing | frame 0 | later frames |
|---|---|---|---|
| `.tpr`     | 100 ps  | 4.21 A | **190-194 A** |
| `.tpr`     | 1000 ps | 4.21 A | **190-196 A** |
| `ref.gro`  | 100 ps  | 4.21 A | 4.02-4.19 A |
| `ref.gro`  | 1000 ps | 4.21 A | 4.09-4.20 A |

Same input trajectory, same stride, only `-s` differs: with the `.tpr` the
chains come apart a few tens of frames in and stay broken. The `.tpr` holds the
grompp input configuration, not frame 0 of this trajectory, and feeding nojump a
reference that is a different configuration is what breaks it. Always build the
reference with `-dump 0 -pbc whole` from the very trajectory you are about to
unwrap.

A3 also drops `-n`; the single group offered is called `System` there — it is
the 1180 protein beads.

**Stop at A3.** No `-center`, no `-pbc mol`, no `-ur compact`.

### A4) Optional smoothing for a movie

```bash
# -fit asks for TWO groups: the least-squares fit group, then the output group.
# A single "echo 0" leaves the second prompt with an empty stdin and trjconv
# dies with "Fatal error: Cannot read from input".
printf '0\n0\n' | $G trjconv -s $R/vis_prot_ref.gro -f $R/vis_prot_nojump.xtc \
    -o fit.xtc -fit progressive
```

`-fit progressive` aligns each frame to the previous one — continuous by
construction, the smoothest result. **Never use it for analysis**, it destroys
the real motion.

`-fit rot+trans` over all 5 chains is meaningless while they are dispersed: it
fits a rigid body to a cloud that has no rigid-body motion. If you want it, fit
on one chain (`Protein1`) only.

## Branch B — keep the chains visually together

Only if you specifically want a view where the chains do not drift apart.
Centre on one chain and wrap the rest to their minimum image:

```bash
printf '11\n0\n' | $G trjconv -s md_chunk_000.tpr -f all.xtc -n index.ndx \
    -pbc mol -ur compact -center -o centered.xtc
    # 11 = Protein1 (centring), 0 = AllProteins (output)
```

`Protein1` is fixed and continuous; the other four **will** jump when they cross
a cell boundary. That is unavoidable — with freely diffusing molecules no single
unwrapping keeps every chain simultaneously continuous in time *and* in the same
periodic image as its neighbours. Branch A gives continuity, branch B gives
proximity; you cannot have both.

`-pbc cluster` is not an option unless the protomers have actually aggregated —
the GROMACS docs state it "will only give meaningful results if you in fact have
a cluster". At the end of chunk 000 the five COMs needed four mutually
inconsistent image shifts (no cluster), and `cluster` re-decides them every
frame, which is where the flicker comes from.

## 2) Index groups — two traps

Always pass `-n index.ndx`. Without it trjconv falls back to the tpr's default
groups, which are not the ones below.

```
 0  AllProteins    1180     9  Protein        1180     16  SideChain    1180  (!)
 1  Backbone          0    10  Protein-H      1180     18  Solvent...  42169
 2  C-alpha           0    11  Protein1        236     19  System      43349
 4  MainChain         0    12  Protein2        236     20  W           41261
```

* `Backbone`, `C-alpha`, `MainChain*` are **empty** — `make_ndx` builds them from
  atomistic atom names, which do not exist in Martini. Selecting one silently
  produces an empty trajectory.
* `SideChain` and `Protein-H` hold all 1180 beads, not side chains. Equally
  meaningless. Use `AllProteins` or `Protein1`..`Protein5`.

`AllProteins` happens to be atoms 1..1180, contiguous, so a stripped trajectory
still lines up with the full-system index. Do not rely on that in other systems.

## 3) VMD

```bash
vmd vis_prot_ref.gro vis_prot_nojump.xtc
```

**Use the `.gro` from step A2, never `runs/<RUN>/md_chunk_000.gro`.** Two
reasons, either one fatal:

* *Atom count.* The pipeline strips solvent, so the `.xtc` holds 1180 beads
  while the run's `.gro` holds 43349. VMD refuses the trajectory outright —
  `ERROR) Incorrect number of atoms (1180) in coordinate file` — and leaves you
  with a single static frame.
* *Wrong frame.* `md_chunk_000.gro` is what mdrun writes at the **end** of the
  run; its title carries no `t=`, while `vis_prot_ref.gro` says
  `t= 0.00000 step= 0`. Even at matching atom counts it would place the final
  configuration at position 0 of the animation, in a different periodic image
  from the rest of the trajectory.

The rule: the structure file must come from the same trjconv pipeline as the
trajectory — same index group, same pbc treatment, `-dump 0`.

To keep the solvent, re-run A1-A3 selecting `System` instead of `AllProteins`
everywhere; the matching reference is then a full-system `-dump 0 -pbc whole`
`.gro`. Note that over 100 ns the water diffuses far outside the cell in a
nojump trajectory — for looking at solvent, branch B is usually the better view.

The `.gro` loads as frame 0 and the `.xtc` supplies its own frame 0, so the
animation has one duplicated frame at the start (1002 for 1001);
`animate delete beg 0 end 0 top` removes it.

**Martini bonds.** VMD guesses bonds by distance, which is wrong for CG beads and
makes molecules look exploded or flickering — easy to mistake for a PBC problem.
Load the real topology:

```tcl
source cg_bonds-v5.tcl
cg_bonds -tpr /data/runs/<RUN>/md_chunk_000.tpr
```

`cg_bonds-v5.tcl` is on the Martini site; `MartiniGlass` (JCIM 2025) is the newer
Python generator for the same purpose.

**Do not wrap.** With a branch-A trajectory never call `pbc wrap` — it undoes A3.
Only draw the cell:

```tcl
pbc box -center bb -style dashed
```

Starting from a raw `.xtc` instead, PBCTools does the same job VMD-side:
`pbc join fragment -all` then `pbc unwrap -sel "name BB SC1"`.

**Residual jerkiness is sampling, not PBC.** `nstxout-compressed = 10000` with
`dt = 0.01` is one frame per 100 ps, and Martini dynamics run ~4x fast. Either
set Graphical Representations -> Trajectory -> **Smoothing Window Size** to 3-5
(a display filter, it does not touch the data), or re-run a short 5-10 ns chunk
with `nstxout-compressed = 500` for a genuinely fluid movie.

## 4) Caveat for analysis, not visualization

`gmx trjconv -pbc nojump` implements the heuristic "HLAT" unwrapping scheme,
which Bullerjahn et al. (JCTC 2023, doi:10.1021/acs.jctc.3c00308) show
"occasionally unwraps particles into the wrong box, which results in an
artificial speed up of the particles" under a barostat. These runs use
`pcoupl = C-rescale`, so the trajectories above are fine to look at but must
**not** be used to compute MSD or diffusion coefficients. For that, unwrap the
*wrapped* trajectory with their displacement-based TOR scheme instead.

## 5) Checking a trajectory

Two independent checks — a file can pass one and fail the other.

**Continuity** (per-frame COM displacement per chain; should be well under 1 nm):

```bash
python - <<'PY' > /tmp/prot.ndx
for c in range(5):
    print(f"[ P{c+1} ]")
    ids=[str(i) for i in range(c*236+1,(c+1)*236+1)]
    for k in range(0,len(ids),15): print(" ".join(ids[k:k+15]))
PY
echo "0 1 2 3 4" | $G traj -f vis_prot_nojump.xtc -s vis_prot_ref.gro \
    -n /tmp/prot.ndx -com -ng 5 -ox com.xvg
```

**Wholeness** (largest gap between consecutive backbone beads; should be ~0.45 nm,
never metres):

```bash
echo 0 | $G trjconv -s vis_prot_ref.gro -f vis_prot_nojump.xtc -dt 2000 -o probe.gro
python - <<'PY'
import math
L=open('probe.gro').read().splitlines(); i=0; worst=0
while i < len(L)-1:
    n=int(L[i+1]); body=L[i+2:i+2+n]
    for c in range(5):
        bb=[(float(a[20:28]),float(a[28:36]),float(a[36:44]))
            for a in body[c*236:(c+1)*236] if a[10:15].strip()=='BB']
        worst=max(worst, max(math.dist(bb[k],bb[k+1]) for k in range(len(bb)-1)))
    i+=2+n+1
print(f"max consecutive BB-BB distance: {worst:.2f} nm")
PY
```
