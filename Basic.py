import numpy as np
import matplotlib.pyplot as plt

def simulate_geomagnetic_reversals(rate, time_span):

    reversal_times = []
    current_time = 0
    
    while current_time < time_span:
        wait_time = np.random.exponential(1 / rate)  # Poisson process: Exponential waiting times
        current_time += wait_time
        if current_time < time_span:
            reversal_times.append(current_time)
    
    return np.array(reversal_times)

def plot_magnetochron(reversal_times, time_span):
    fig, ax = plt.subplots(figsize=(10, 2))
    
    # Initial state (assume normal polarity at the start)
    times = np.concatenate(([0], reversal_times, [time_span]))
    colors = ["black" if i % 2 == 0 else "white" for i in range(len(times) - 1)]
    
    for i in range(len(times) - 1):
        ax.fill_betweenx([0, 1], times[i], times[i+1], color=colors[i])
    
    ax.set_xlim(0, time_span)
    ax.set_yticks([])
    ax.set_xlabel("Time (Million Years)")
    ax.set_title("Simulated Geomagnetic Reversals (Magnetochron Format)")
    plt.show()

# Parameters
mean_reversal_rate = 2  # Average reversals per million years (example value)
time_span = 10  # Total time span in million years

# Run simulation
reversal_times = simulate_geomagnetic_reversals(mean_reversal_rate, time_span)

# Plot results
plot_magnetochron(reversal_times, time_span)
