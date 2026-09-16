import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_2/random.dat', header=None)

# #Question 2a
# hist, bin_edges = np.histogram(df,bins=1000, density = True)
# bin_center = (bin_edges[:-1] + bin_edges[1:])/2
# plt.scatter(bin_center, hist)
# plt.ylim(0,1.5)
# plt.xlabel("random number")
# plt.ylabel("frequency")
# plt.title("Probability Distribution(bins = 1000)")
# plt.show()

#Question 2b
# plt.scatter(df[:-1],df[1:])
# plt.xlabel("$x_{i}$")
# plt.ylabel("$x_{i+1}$")
# plt.title("Scatter plot to check correlation")
# plt.show()


# df = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_2/corr.dat', header=None)
# plt.scatter(range(1, 101),df)
# plt.xlabel('k')
# plt.ylabel("correlation function between $x_{i}$ and $x_{i+k}$")
# plt.ylim(-0.5,0.5)
# plt.show()
