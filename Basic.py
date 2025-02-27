import numpy as np
import matplotlib.pyplot as plt
import random

plt.close("all")

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

def plot_interstate_magnetochron(reversal_times, inter_states, time_span_myr):
    fig, ax = plt.subplots(figsize=(10, 2))
    
    times = np.concatenate(([0], reversal_times, [time_span_myr]))
    colors = ["black" if i % 2 == 0 else "white" for i in range(len(times) - 1)]
    
    for i in range(len(times) - 1):
        ax.fill_betweenx([0, 1], times[i], times[i+1], color=colors[i])
    
    # Add interstate regions as grey dashed areas
    for start, end in inter_states:
        ax.fill_betweenx([0, 1], start, end, color='grey', alpha=0.5)
    
    ax.set_xlim(0, time_span_myr)
    ax.set_yticks([])
    ax.set_xlabel("Time (Million Years)")
    ax.set_title("Geomagnetic Reversals with Interstate Periods")
    plt.show()

def plot_interstate_magnetochron_with_gaps(reversal_times, inter_states, gaps_pos, fixed_gap_length, time_span_myr, variable_lengths=None):
    fig, ax = plt.subplots(figsize=(10, 2))
    
    times = np.concatenate(([0], reversal_times, [time_span_myr]))
    colors = ["black" if i % 2 == 0 else "white" for i in range(len(times) - 1)]
    
    for i in range(len(times) - 1):
        ax.fill_betweenx([0, 1], times[i], times[i+1], color=colors[i])
    
    # Add interstate regions as grey dashed areas
    for start, end in inter_states:
        ax.fill_betweenx([0, 1], start, end, color='grey', alpha=0.5)
    
    # Add red gap regions with fixed or variable lengths
    for i, center in enumerate(gaps_pos):
        length = fixed_gap_length if variable_lengths is None else variable_lengths[i]
        start = center - length / 2
        end = center + length / 2
        ax.fill_betweenx([0, 1], start, end, color='red', alpha=0.5)
        
    if (variable_lengths is not None):
        title_name = "Geomagnetic Reversals with Interstate Periods and variable Gaps"
    else:
        title_name = "Geomagnetic Reversals with Interstate Periods and fixed Gaps"
    
    ax.set_xlim(0, time_span_myr)
    ax.set_yticks([])
    ax.set_xlabel("Time (Million Years)")
    ax.set_title(title_name)
    plt.show()
    
def plot_all_in_one(reversal_times, inter_states, gaps_pos, fixed_gap_length, time_span_myr, variable_lengths):
    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    titles = ["Simulated Geomagnetic Reversals", 
              "Geomagnetic Reversals with Interstate Periods", 
              "Geomagnetic Reversals with Interstate Periods and Fixed Gaps", 
              "Geomagnetic Reversals with Interstate Periods and Variable Gaps"]
    
    for i, ax in enumerate(axes):
        times = np.concatenate(([0], reversal_times, [time_span_myr]))
        colors = ["black" if i % 2 == 0 else "white" for i in range(len(times) - 1)]
        
        for j in range(len(times) - 1):
            ax.fill_betweenx([0, 1], times[j], times[j+1], color=colors[j])
        
        # Add interstate regions as grey dashed areas in all but first plot
        if i > 0:
            for start, end in inter_states:
                ax.fill_betweenx([0, 1], start, end, color='grey', alpha=0.5)
        
        # Add red gap regions in third and fourth plots
        if i == 2:  # Fixed gap length
            for center in gaps_pos:
                start = center - fixed_gap_length / 2
                end = center + fixed_gap_length / 2
                ax.fill_betweenx([0, 1], start, end, color='red', alpha=1)
        elif i == 3:  # Variable gap lengths
            for j, center in enumerate(gaps_pos):
                start = center - variable_lengths[j] / 2
                end = center + variable_lengths[j] / 2
                ax.fill_betweenx([0, 1], start, end, color='red', alpha=1)
        
        ax.set_xlim(0, time_span_myr)
        ax.set_yticks([])
        ax.set_title(titles[i])
    
    axes[-1].set_xlabel("Time (Million Years)")
    plt.tight_layout()
    plt.show()

# User-defined parameters
time_span_myr = 5  # Total time in million years
mean_reversal_rate = 4  # Reversals per million years
min_gap_years = 10000  # Minimum gap between reversals in years
changing_state_time = 50000  # Time in years the field is in an intermediate state
fixed_gap_length = 10000 # Fixed gap length in years

mean_gap_length = 20000
mean_gap_length_std = 10000

min_gaps = 10
max_gaps = 15

gaps_num = random.randint(min_gaps, max_gaps)
gap_lengths = np.random.normal(mean_gap_length, mean_gap_length_std, gaps_num) / 1e6  # Convert to million years
gaps_pos = np.random.uniform(0,time_span_myr,gaps_num)

changing_state_time /= 1e6  # Convert to million years
fixed_gap_length /= 1e6
mean_gap_length /= 1e6
mean_gap_length_std /= 1e6

# Run simulation
reversal_times = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, min_gap_years)

# Optimize interstates calculation using vectorized operations
inter_states = np.column_stack((reversal_times - changing_state_time / 2, 
                                reversal_times + changing_state_time / 2))

# Plot results
#plot_magnetochron(reversal_times, time_span_myr)
#plot_interstate_magnetochron(reversal_times, inter_states, time_span_myr)

#plot_interstate_magnetochron_with_gaps(reversal_times, inter_states, gaps_pos, fixed_gap_length, time_span_myr)
#plot_interstate_magnetochron_with_gaps(reversal_times, inter_states, gaps_pos, fixed_gap_length, time_span_myr, variable_lengths=gap_lengths)
plot_all_in_one(reversal_times, inter_states, gaps_pos, fixed_gap_length, time_span_myr, gap_lengths)