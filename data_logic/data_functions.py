# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: February, 2016.

# LIBRARIES
import numpy as np
import scipy.stats.stats as sss


# ==============
# DATA FUNCTIONS
# ==============

# FUNCTION: data_switcher
def data_switcher(a_switch_log, a_switch_zscore, rep):
    if a_switch_log:
        rep = log_function(rep)
    if a_switch_zscore:
        rep = zscore_function(rep)
    return rep


# FUNCTION: log
def log_function(rep):
    rep = map(lambda x: np.array([np.log(item) for item in x]), rep)
    return list(rep)


# FUNCTION: zcore
def zscore_function(rep):
    rep = map(lambda x: np.transpose(np.array([sss.zscore(item) for item in np.transpose(x)])), rep)
    return list(rep)
