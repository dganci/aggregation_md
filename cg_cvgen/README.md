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

## Outputs

```text
CVs_torchscript.pt     TorchScript model for PLUMED PYTORCH_MODEL
cv_params.pkl          MetaD-compatible CV parameters
cv_params.json         Human-readable parameters
cv_manifest.json       Stable handoff contract for the MD package
feature_schema.json    Ordered feature list used during training
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

If the MD package supports a conventional CV directory, the equivalent is:

```bash
./build/cg_md --stage metadynamics --metadynamics --cv-dir ../runs/10xdesmin_head_pS31/CVs
```

The system never invents neural CVs. If explicit files are omitted, the MD package must find `CVs_torchscript.pt` and `cv_params.pkl` in the configured/default CV directory.
