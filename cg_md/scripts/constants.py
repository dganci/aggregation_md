"""Two things every metric module needs and none should own."""

import numpy as np

KB_KJ_PER_MOL_K = 0.008314462618

trapezoid = getattr(np, "trapezoid", None) or np.trapz
