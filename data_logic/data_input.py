# LIBRARIES
import numpy as np

"""
This set of replicate data is produced by GenData.m (John&Norris)
with p = 0.9 (correlation) and myseed = 0.

BE1, BE2, BE3 and BE4 are highly correlated.
BE5, BE6, BE7 and BE8 are highly correlated.
BE9, BE10, BE11 and BE12 are highly correlated.

Nine biological entities and five times
"""

a_name = "Stat Simulated Set B (3 Reps)"
a_protein_names = ["BE1", "BE2", "BE3", "BE4", "BE5", "BE6", "BE7", "BE8", "BE9", "BE10", "BE11", "BE12"]

# A1, A2 and A3 are 10X12 matrices, each row represents Time And each column represents a gene.
a1 = np.matrix([[-0.2661, -0.2626, -0.6433, -0.3175, 0.3228, 0.2400, 0.2066, -0.0202, -0.6299, -0.7205, -0.2354, 0.1648],
               [0.6691, 0.7789, 0.8018, 0.6711, 1.9862, 1.1671, 1.3754, 1.4281, 0.6615, 0.6932, 1.2144, 1.4087],
               [-0.6429, -0.3077, -0.2047, -0.2941, -0.1656, -0.2159, -0.2763, -0.8702, -0.6379, -0.8224, -1.0898, -0.6721],
               [0.1612, 0.2670, 0.0785, 0.2360, 1.1760, 1.2615, 1.6320, 1.1136, -0.4862, -0.2013, -0.7602, -1.0125],
               [0.1346, 0.0021, 0.0201, 0.1648, 0.6069, 0.9519, 1.1732, 1.3939, -0.7848, -0.2662, -0.6466, -1.2426],
               [0.8943, 1.0763, 1.0509, 1.5519, -0.1021, -0.1147, -0.2551, -0.3097, -0.6358, -0.7667, -0.7451, -0.5511],
               [-0.0514, -0.1906, 0.3919, 0.5753, -0.0334, -0.0854, -0.2285, -0.8379, 0.2725, 0.0424, 0.3760, -0.0489],
               [0.3345, -0.1412, 0.2092, 0.2328, 2.0073, 1.3516, 1.7328, 1.7239, -0.2620, -0.0963, -0.0756, -0.1581],
               [-0.7090, -0.5052, -0.0511, -0.1281, -0.5949, -0.3089, -0.0123, 0.1103, 0.8824, 0.6439, 0.5111, 0.4841],
               [-0.4162, -0.7906, -1.3446, -1.4964, -2.4507, -1.8441, -1.6639, -1.4802, 0.7134, 0.5876, 0.5868, 0.6281]])

a2 = np.matrix([[-0.7612, -0.7883, -0.6785, -0.9226, -0.6354, -0.0823, -0.2506, -0.1936, 1.3030, 0.8455, 1.2651, 0.7278],
               [-0.1785, -0.0862, 0.3871, 0.1094, 0.7278, 0.6072, 0.0491, 0.4727, 0.3305, -0.1081, -0.2017, 0.0821],
               [0.5407, -0.1403, 0.4301, 0.6738, 0.5464, 0.2520, 0.2069, 0.4790, 0.2844, -0.5119, -0.4460, -0.6602],
               [1.5499, 2.1117, 2.4988, 1.8139, -1.2409, -1.3755, -0.8864, -1.0047, 3.4543, 3.4543, 2.9854, 3.0441],
               [-0.4955, -0.1515, -0.3941, -0.0431, -0.9351, -0.6828, -0.8531, -0.1337, -0.3966, -0.0363, -0.1961, 0.1383],
               [-0.9278, -1.3548, -1.5932, -1.4410, -0.4490, -0.8777, -0.4198, -0.7993, 0.9905, 1.0479, 0.1260, 0.1513],
               [-0.0929, -0.1074, -0.1561, 0.6938, 0.3870, -0.4728, -0.4245, -0.6551, 0.1956, -0.1754, -0.1391, -0.1474],
               [-0.2849, -0.1246, -0.4313, 0.3400, 0.5782, -0.0697, -0.2645, -1.0503, 1.0305, 0.4124, -0.1001, -0.1666],
               [-0.6317, -0.1261, -0.7340, -0.2645, 1.0818, 0.9110, 1.0916, 1.0355, -0.2227, -0.5560, -0.6784, -0.4526],
               [0.2362, 0.2499, 1.1525, 1.1781, -1.3257, -1.4328, -0.8431, -0.8595, 0.7892, 0.9708, 1.3949, 1.7278]])

a3 = np.matrix([[-0.3037, -0.1319, -0.3548, -0.8087, -0.4033, -0.5092, -0.9379, -0.9475, 1.4610, 1.3728, 1.3414, 1.0966],
                [0.1352, 0.3214, -0.0928, 0.6711, 1.0562, 1.3608, 1.7692, 1.4852, -0.3569, -0.5944, -0.4883, -0.7396],
                [0.9714, 1.5624, 1.7980, 1.3441, 0.6450, 1.0425, 1.2587, 1.4933, -1.8839, -1.0112, -1.1005, -0.7493],
                [-0.2349, 0.0986, 0.4457, 0.2786, 0.3080, -0.0251, 0.1724, 0.0127, 0.8407, 0.8135, 1.0810, 0.9515],
                [-1.9821, -2.2717, -2.0839, -2.0991, -0.3238, -0.5609, 0.3326, -0.0750, -0.2390, -0.1517, -0.3918, -0.4087],
                [-0.4556, -0.5551, -0.8019, -0.3289, 1.4006, 0.5875, 0.9213, 0.3090, 0.5715, 0.9193, 0.6079, 0.7242],
                [0.3840, -0.0040, 0.2416, 0.1417, -0.7404, -0.7390, -0.6517, -1.1019, 0.5963, 0.7686, 0.0148, -0.0481],
                [1.8277, 1.6887, 1.2728, 1.2458, 0.3785, 0.1710, 0.0240, 0.0360, -0.9723, -0.6729, -0.4827, -1.2915],
                [-0.3325, 0.5144, 0.4440, 0.1577, 0.2143, 0.1138, 0.5748, 0.0893, -0.1129, 0.1674, 0.5880, 0.7355],
                [2.3348, 2.1127, 1.8327, 1.3711, 0.1261, -0.1367, -0.1846, -0.0825, -0.3075, 0.3674, 0.1265, -0.0358]])


# The Stat Simulated Set A times are artificial, indicated here by 1, 2, 3, 4 and 5:
a_times = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Set switches for internally selected transforms
# Do NOT apply log transform
# Do apply the Zscore transform
a_switch_log = False
a_switch_zscore = True
