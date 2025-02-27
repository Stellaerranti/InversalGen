import numpy as np
import matplotlib.pyplot as plt

def simulate_geomagnetic_reversals(rate_per_myr, time_span_myr, min_gap_years=500):
    reversal_times = []
    current_time = 0
    min_gap_myr = min_gap_years / 1e6  # Convert years to million years
    
    while current_time < time_span_myr:
        wait_time = np.random.exponential(1 / rate_per_myr)  # Poisson process: Exponential waiting times
        if wait_time < min_gap_myr:
            continue  # Ensure minimum gap is maintained
        current_time += wait_time
        if current_time < time_span_myr:
            reversal_times.append(current_time)
    
    return np.array(reversal_times)

def plot_magnetochron(reversal_times, time_span_myr):

    fig, ax = plt.subplots(figsize=(10, 2))
    
    times = np.concatenate(([0], reversal_times, [time_span_myr]))
    colors = ["black" if i % 2 == 0 else "white" for i in range(len(times) - 1)]
    
    for i in range(len(times) - 1):
        ax.fill_betweenx([0, 1], times[i], times[i+1], color=colors[i])
    
    ax.set_xlim(0, time_span_myr)
    ax.set_yticks([])
    ax.set_xlabel("Time (Million Years)")
    ax.set_title("Simulated Geomagnetic Reversals")
    plt.show()

# User-defined parameters
time_span_myr = 1  # Total time in million years
mean_reversal_rate = 5  # Reversals per million years
min_gap_years = 20000  # Minimum gap between reversals in years

# Run simulation
reversal_times = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, min_gap_years)

# Plot results
plot_magnetochron(reversal_times, time_span_myr)
