import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# df = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_2/exprandom.dat', header=None)
# binwidth = 0.1
# hist, bin_edges = np.histogram(df,bins=int((df.max()-df.min())/binwidth))
# bin_center = (bin_edges[:-1] + bin_edges[1:])/2
# plt.scatter(bin_center, hist)
# plt.xlabel('random number')
# plt.ylabel("frequency")
# plt.show()

df = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_2/gaus.dat', header=None)
binwidth = 0.1
hist, bin_edges = np.histogram(df,bins=int((df.max()-df.min())/binwidth))
bin_center = (bin_edges[:-1] + bin_edges[1:])/2
plt.scatter(bin_center, hist)
plt.xlabel('random number')
plt.ylabel("frequency")
plt.show()


