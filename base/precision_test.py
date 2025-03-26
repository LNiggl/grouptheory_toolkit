import numpy as np
from base import groups as g
from base import objects as o
from base import representations as r
from base import testing as t

nums = [0,1/np.sqrt(2),1/np.sqrt(4),1/np.sqrt(6),1/np.sqrt(8),1/np.sqrt(12),1/np.sqrt(16),1/np.sqrt(24),1/np.sqrt(32),1/np.sqrt(48),1/np.sqrt(96)]
t.test_largest_deviations_from_set("D:/Master/Masterarbeit/results/twopi/data/twopipqr_vecdata.npy",nums)