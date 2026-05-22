# aggregation_md

A C++ toolkit for protein aggregation molecular dynamics simulations, supporting deep learning-based slow collective variable generation and metadynamics-enhanced sampling.

Installation:
```bash
git clone https://github.com/dganci/aggregation_md
cd aggregation_md/cg_md/build
cmake ..
cmake --build . -j
```

Examples of usage:
```bash
i) Non-adaptive version
./cg_md --project-dir .. --n-prot 10 --protomer desmin_head_pS31 --seq-length 108 --atoms-per-prot 237 --elastic-units 90:102 --packmol-cluster-radius-A 145 --packmol-tolerance-A 4.0 --packmol-radius-A 5.0 --box-margin-nm 1.8 --solvation-mode insane --salt-M 0.15 --nvt-dt-ps 0.01 --nvt-nsteps 100000 --npt-dt-ps 0.01 --npt-nsteps 1000000 --md-dt-ps 0.01 --md-total-us 0.0001 --plumed-stride 200 --thermostat legacy

ii) Adaptive version
./cg_md --project-dir .. --n-prot 10 --protomer desmin_head_pS31 --phospho --seq-length 108 --atoms-per-prot 237 --elastic-units 90:102 --packmol-cluster-radius-A 145 --packmol-tolerance-A 4.0 --packmol-radius-A 5.0 --box-margin-nm 1.8 --solvation-mode insane --salt-M 0.15 --nvt-dt-ps 0.01 --nvt-nsteps 100000 --npt-dt-ps 0.01 --npt-nsteps 1000000 --md-dt-ps 0.01 --plumed-stride 200 --thermostat legacy --adaptive --adaptive-chunk-us 0.5 --adaptive-max-total-us 10.0 --adaptive-min-total-us 2.0 --adaptive-min-cn-range 10 --adaptive-min-rg-range 0.3 --adaptive-min-contact-patterns 5 --adaptive-min-lcc-unique 3
```
