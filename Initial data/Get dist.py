import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Load data
sample = np.loadtxt('Real Duration.txt')


dist = stats.gamma
bounds = [
    (0.1, 10),    # Shape (α): Avoid 0, constrain near 1 for Poisson-like behavior
    (-1, 1),       # Location (loc): Small finite range (often 0 for durations)
    (0.1, 10)      # Scale (β): Strictly positive
]
fitRes = stats.fit(dist, sample, bounds)
alpha, loc, beta = fitRes.params
print(f"Fitted Gamma: α={alpha:.3f}, loc={loc:.3f}, β={beta:.3f}")


ks_statistic, p_value = stats.kstest(
    sample, 'gamma', args=(alpha, loc, beta)
)
print(f"KS Statistic: {ks_statistic:.3f}")
print(f"P-value: {p_value:.3f}")


plt.figure(figsize=(10, 6))
ecdf = stats.ecdf(sample)
ecdf.cdf.plot(plt.gca(), label='Empirical CDF')

x = np.linspace(min(sample), max(sample), 1000)
gamma_cdf = stats.gamma.cdf(x, alpha, loc, beta)
plt.plot(x, gamma_cdf, label=f'Gamma CDF (α={alpha:.3f}, β={beta:.3f})', color='red')

plt.xlabel('Duration')
plt.ylabel('Cumulative Probability')
plt.title('ECDF vs. Gamma CDF Fit')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()