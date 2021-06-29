import numpy as np

table = np.genfromtxt('kymotracks.txt', skip_header=1, delimiter=';')

position = table[:, 1]
time = table[:, 2]
