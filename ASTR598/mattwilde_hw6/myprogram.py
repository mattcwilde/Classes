from __future__ import print_function
from functools import partial

import numpy as np
import mergesort
import timeit

N = [100,10000,1000000,100000000]
# N = [10,20,30,40]
time = []

for n in N:
	arr = np.random.randint(low=0, high=n, size=n)
	ave_t = np.mean(timeit.Timer(partial(mergesort.mergesort,arr)).repeat(repeat=3, number=1))
	print(n, ave_t, 'seconds')
	time.append(ave_t)



import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(8,8))

# const = 2.5e-6
theory = N*np.log10(N)

ax.loglog(N,time,lw=1,label=r'Computed Time')
ax.loglog(N,theory,lw=1,label=r'NlogN')

ax.set_xlabel('Number of Ints')
ax.set_ylabel('Time [seconds]')
ax.set_title('Mergesort Scaling')
ax.legend(loc='upper left')
fig.savefig('mergesort_time.png')