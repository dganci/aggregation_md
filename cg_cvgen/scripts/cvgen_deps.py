"""The training stack, imported on demand."""
import os
import sys


class _Deps:
    pd = torch = nn = mlcolvar = L = EarlyStopping = Callback = create_timelagged_dataset = None


deps = _Deps()


def load():
    _mlc = os.environ.get("CG_CVGEN_MLCOLVAR_PATH")
    import pandas
    import torch
    import mlcolvar
    import lightning.pytorch
    from lightning.pytorch.callbacks import EarlyStopping, Callback
    from mlcolvar.utils.timelagged import create_timelagged_dataset
    deps.pd, deps.torch, deps.nn, deps.mlcolvar = pandas, torch, torch.nn, mlcolvar
    deps.L, deps.EarlyStopping, deps.Callback = lightning.pytorch, EarlyStopping, Callback
    deps.create_timelagged_dataset = create_timelagged_dataset

    if _mlc:
        if _mlc not in sys.path:
            sys.path.insert(0, _mlc)
        from mlcolvar_copy.utils.timelagged import create_timelagged_dataset as patched
        deps.create_timelagged_dataset = patched
        print(f"time-lagged pairs from the patched mlcolvar in {_mlc} (walker-aware)")
