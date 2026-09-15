# cg_cvgen

Standalone C++ front-end for generating DeepTICA collective variables from an unbiased PLUMED COLVAR file.

The package is intentionally separate from the MD runner. It consumes the unbiased trajectory descriptors and produces the files needed by the biased MetaD/OPES workflow.

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Runtime Python dependencies

```bash
python3 -m pip install -r requirements.txt
```

The C++ executable is dependency-light; the DeepTICA training is executed by the versioned backend in `scripts/cvgen_backend.py`.

## Minimal run

```bash
./build/cg_cvgen \
  --input-colvar ../runs/10xdesmin_head_pS31/COLVAR \
  --output-dir ../runs/10xdesmin_head_pS31/CVs \
  --lags 5,7,10,12 \
  --n-cvs 3 \
  --equilibration-time-ps 50 \
  --max-epochs 200
```

### Training on an adaptive run

An adaptive `cg_md` run does not produce one COLVAR but a numbered series that
together form a single continuous trajectory. `--input-colvar` accepts a glob;
the backend concatenates the files in simulated-time order and drops the frame
PLUMED reprints at each chunk boundary (a zero-length time gap would otherwise
corrupt the time-lagged pairs DeepTICA is built on):

```bash
./build/cg_cvgen \
  --input-colvar '../runs/10xdesmin_head_pS31/COLVAR_chunk_*.dat' \
  --output-dir ../runs/10xdesmin_head_pS31/CVs
```

Quote the pattern so the shell does not expand it.

### How the lag is chosen

`--lags` is a fallback. By default the backend first runs a linear-TICA
scan over a geometric ladder of lags (5 to 8000 frames, capped at a tenth of
the shortest replica) and reads the slowest implied timescale off each. The
lag has to sit where that timescale stops depending on it, and that region
moves with the system, so no fixed list is right everywhere.

The plateau is the longest run of consecutive lags whose timescales agree to
within 3% of the run's own median; between runs of equal length the one at
shorter lags wins, since it has more pairs behind it. The scan is written to
`its_scan.csv` beside the model, with the plateau and the selected lags
marked, and `analyze_offline.py cv` in `cg_md/scripts` draws it.

Two things to read off that file before trusting the model. If no plateau was
found, the training fell back to `--lags` and said so in the log: the run is
too short for its own slowest mode. And if the slowest timescale is of the
same order as the run that produced it, it is the bound the data can support
rather than a measurement (Sinitskiy & Pande 2018); `cg_md`'s adaptive
production checks exactly this after every chunk (`cv_readiness.txt`) and
does not stop until the run is `--adaptive-min-time-over-its` times longer
than that timescale.

## Where the Python lives

`scripts/cvgen_backend.py` is the entry point the binary runs; the work is one
module per stage, each importable on its own: `cvgen_data.py` (reading the
COLVARs onto one clock), `cvgen_features.py` (the permutation-invariant
feature vector and the check that it is the whole one), `cvgen_plateau.py`
(the implied-timescale scan), `cvgen_train.py` (one DeepTICA model per lag),
`cvgen_score.py` (ranking the lags), `cvgen_export.py` (TorchScript model,
`cv_params`, manifest). torch, lightning, mlcolvar and pandas are imported on
demand by `cvgen_deps.py`, so `cvgen_clock.py`, `cvgen_metrics.py` and
`cvgen_plateau.py` - and the tests - need only numpy.

## Outputs

```text
CVs_torchscript.pt     TorchScript model for PLUMED PYTORCH_MODEL
cv_params.pkl          MetaD-compatible CV parameters
cv_params.json         Human-readable parameters
cv_manifest.json       Stable handoff contract for the MD package
feature_schema.json    Ordered feature list used during training
its_scan.csv           Implied-timescale scan, plateau and selected lags
lag_scores.csv         Lag ranking and diagnostics
lag_scores.json        Full diagnostics
CVs_model_lagX.pt      PyTorch checkpoint
model_lagX.pkl         Pickled DeepTICA model
```

## Handoff to the biased MD package

Pass the generated directory or explicit files:

```bash
./build/cg_md \
  --stage metadynamics \
  --metadynamics \
  --metad-model ../runs/10xdesmin_head_pS31/CVs/CVs_torchscript.pt \
  --metad-cv-params ../runs/10xdesmin_head_pS31/CVs/cv_params.pkl
```

There is no `--cv-dir`: `cg_md` takes the two files individually, or falls back
to a fixed location. With `--metad-model` / `--metad-cv-params` omitted it looks
for them under the run's own metadynamics directory:

```
runs/<system>_<runid>/metadynamics/CVs/<protomer>/CVs_torchscript.pt
runs/<system>_<runid>/metadynamics/CVs/<protomer>/cv_params.pkl
```

(`Config::metadModelPath()` / `metadCvParamsPath()`; for `--stage pmf` the
`metadynamics` component becomes `pmf`, though a PMF needs no model at all.)
Since `<runid>` is a hash of the run's parameters, that path is only convenient
when the CVs were trained for that exact run — otherwise pass the two files
explicitly, which is the usual case.

The system never invents neural CVs: if neither the options nor the fallback
path yield the files, the stage stops rather than proceeding unbiased.

`cv_params.pkl` is the contract: `cg_md` reads `feature_cols`, `sigma`,
`grid_min`/`grid_max` and `n_cvs` from it in one pass, uses `feature_cols`
verbatim as the `PYTORCH_MODEL` argument list, and refuses to start if the
model names a descriptor the system's `plumed.dat` does not define or if
`n_cvs` disagrees with `--metad-nodes`. Consequently the CV must be trained on
a COLVAR from the *same* `--n-prot` / descriptor set that will be biased.
