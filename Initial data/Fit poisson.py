import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Load data
sample = np.loadtxt('Real Duration.txt')

# Fit exponential distribution: p(x) = λ * exp(-λx)
# In scipy: expon has parameters (loc, scale), where scale = 1/λ
loc, scale = stats.expon.fit(sample)  # fix loc=0 to match theory
lambda_est = 1 / scale

print(f"Fitted Exponential: λ = {lambda_est:.3f}")

# KS test
ks_statistic, p_value = stats.kstest(sample, 'expon', args=(loc, scale))
print(f"KS Statistic: {ks_statistic:.3f}")
print(f"P-value: {p_value:.3f}")

# Plot ECDF vs Exponential CDF
plt.figure(figsize=(10, 6))

# Empirical CDF
sample_sorted = np.sort(sample)
ecdf = np.arange(1, len(sample) + 1) / len(sample)
plt.plot(sample_sorted, ecdf, label='Empirical CDF')

# Theoretical Exponential CDF
x = np.linspace(min(sample), max(sample), 1000)
exp_cdf = stats.expon.cdf(x, loc=loc, scale=scale)
plt.plot(x, exp_cdf, label=f'Exponential CDF (λ={lambda_est:.3f})', color='red')

plt.xlabel('Duration')
plt.ylabel('Cumulative Probability')
plt.title('ECDF vs. Exponential CDF Fit')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()
