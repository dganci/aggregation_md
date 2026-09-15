# Upstream mlcolvar generates this file at build time; this vendored copy was
# taken from the source tree without it, which made the package unimportable
# ("No module named 'mlcolvar_copy._version'") and left the local patch below
# unreachable.
#
# What the patch adds: `walker` support in utils/timelagged.py, i.e. the ability
# to say which frames belong to which trajectory so that no time-lagged pair
# spans two independent replicas. The mlcolvar released on PyPI (1.3.1, the one
# installed in the image) does NOT accept that argument.
#
# Not used by default. cg_cvgen imports the installed mlcolvar unless
# CG_CVGEN_MLCOLVAR_PATH points here - see cvgen_backend.py. Switching the
# library that does the training is a decision to take deliberately, not a
# side effect.
__version__ = "1.3.1+local.walker"
