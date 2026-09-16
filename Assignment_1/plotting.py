import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import norm
import statistics
from scipy.optimize import curve_fit 
import math


# #creating plots with bin width 0.5, 1, 2
# binwidth = 5
# df = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_1/sum_10k.dat', header=None)
# print(df)
# hist, bin_edges = np.histogram(df,bins=int((df.max()-df.min())/binwidth),density = True)
# bin_center = (bin_edges[:-1] + bin_edges[1:])/2


# mean = statistics.mean(df[0])
# sd = statistics.stdev(df[0])
# x_axis = np.arange(mean + 5* sd, mean - 5* sd, 1)

# mu, sigma = norm.fit(df)

# x = np.arange(4880, 5120, 0.5)
# y = norm.pdf(x, mu, sigma)
# plt.plot(x,y, "-")
# print(x)
# print(y)

# plt.scatter(bin_center, hist)
# plt.xlabel("Sum")
# plt.ylabel("Frequency")
# plt.title("Histogram of sum of 10000 random numbers with bin size {}".format(binwidth))
# plt.show()
# print(mean)
# print(sd)


#creating plots with bin width 0.5, 1, 2
binwidth = 10
df = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_1/dist_sum_random_walk_1m1m.dat', header=None)
print(df)
hist, bin_edges = np.histogram(df,bins=int((df.max()-df.min())/binwidth),density = True)
bin_center = (bin_edges[:-1] + bin_edges[1:])/2


mean = statistics.mean(df[0])
sd = statistics.stdev(df[0])
x_axis = np.arange(mean + 5* sd, mean - 5* sd, 1)

mu, sigma = norm.fit(df)

x = np.arange(-200, 200, 0.5)
y = norm.pdf(x, mu, sigma)
plt.plot(x,y, "-")


plt.scatter(bin_center, hist)
plt.xlabel("Sum")
plt.ylabel("Frequency")
plt.title("Histogram of sum of 100000 random walks 100000 steps".format(binwidth))
plt.show()
print(mean)
print(sd)



# #creating normalised plots
# df2 = pd.read
# _csv('./norm_sum_10k.dat', header=None)
# plt.hist(df2, bins = 100)
# plt.title("Normalised distribution")
# plt.xlabel("Normalised sum")
# plt.ylabel("Frequncy")
# plt.show()

# df3 = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_1/dist_sum_random_walk.dat', header=None)
# plt.hist(df3, bins=int((df3.max()-df3.min())/10), color = "white", ec="red", density= True)
# plt.xlabel("Walk")
# plt.ylabel("Frequency")
# plt.title("10000 Random Walks Histogram with bin size 10")
# plt.show()

# df4 = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_1/dist_sum_random_walk_1m.dat', header=None)
# plt.hist(df4, bins=int((df4.max()-df4.min())/10), color = "white", ec="red", density= True)
# plt.xlabel("Walk")
# plt.ylabel("Frequency")
# plt.title("100000 Random Walks Histogram with 10000 steps bin size 10")
# plt.show()

# df5 = pd.read_csv('/Users/yash/Documents/Computational Physics/Assignment_1/dist_sum_random_walk_1m1m.dat', header=None)
# plt.hist(df5, bins=int((df5.max()-df5.min())/10), color = "white", ec="red", density= True)
# plt.xlabel("Walk")
# plt.ylabel("Frequency")
# plt.title("100000 Random Walks Histogram with 100000 steps bin size 10")
# plt.show()




############################


# hist, bin_edges = np.histogram(sum, bins=400)
# bin_center = (bin_edges[:-1] + bin_edges[1:])/2
# plt.scatter(bin_center, hist)

# plt.title("Sum of 10000 random numbers between 0 and 1")
# plt.xlabel("Sum")
# plt.ylabel("Frequency”)







