import numpy as np
import matplotlib.pyplot as plt

# Per-cell reaction rate expression
# r_X = k1 * (k2^n + [X]^n) / [X]^n
# r_X / k1 = (1 + [X]^n / k2^n) / ([X]^n / k2^n)

k1 = 1 #5
k2 = 1 #50
n = 1 #3

X = np.linspace(0.1, 10, 100) ** (1/n) * k2
rate = k1 * (k2 ** n + X ** n) / X ** n

plt.figure(constrained_layout=True)
plt.plot((X / k2) ** n, rate / k1, '-k', lw=2)

plt.axhline(1, ls='--', color='r', lw=2)
plt.hlines(2, xmin=-0.5, xmax=1, ls='--', color='b', lw=2)
plt.vlines(1, ymin=-0.5, ymax=2, ls='--', color='b', lw=2)

plt.xlim(left=-0.5)
plt.ylim(bottom=-0.5)
plt.xlabel(r'$[X]^n / k_2^n$', fontsize=16)
plt.ylabel(r'$r_X / k_1$', fontsize=16)
plt.tick_params(labelsize=14)

plt.annotate(r'$r_X = k_1 \cdot \frac{k_2^n + [X]^n}{[X]^n}$', xy=(0.55, 0.85), xycoords='axes fraction',
             fontsize=20)
plt.annotate('1', xy=(0.125, -0.062), color='b', xycoords='axes fraction', fontsize=14)
plt.annotate('1', xy=(-0.04, 0.105), color='r', xycoords='axes fraction', fontsize=14)

plt.show()
