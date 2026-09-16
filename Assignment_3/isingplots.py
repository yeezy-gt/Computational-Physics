import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_3/ising_T2_N40_init_random.dat', delimiter=" ",header=None)

time, m, e = df

plt.plot(m)