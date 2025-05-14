import numpy as np
import matplotlib.pyplot as plt

# Load data from file
data = np.loadtxt('length_200000.0.txt')  # Replace with your filename

# Set fixed bin length (modify this as needed)
bin_length = 100000  # Example: 250-year bins

# Calculate min and max of data
data_min = np.min(data)
data_max = np.max(data)

# Generate bins with fixed width
bins = np.arange(data_min, data_max + bin_length, bin_length)

# Create histogram with density=True to get percentages
counts, bin_edges, patches = plt.hist(data, bins=bins, density=True)

# Convert density to percentage (multiply by bin width)
percentages = counts * bin_length * 100  # bin_length = fixed bin width

# Clear the initial histogram
plt.clf()

# Create bar plot with percentages and outlines
plt.bar(bin_edges[:-1], 
        percentages, 
        width=np.diff(bin_edges), 
        align='edge',
        edgecolor='black',  # Add black outlines
        linewidth=1)        # Outline thickness

xticks = np.arange(0, np.max(bins) + 200000, 200000)  # Every 200 kyr
plt.xticks(xticks, [f"{int(x/1000)}" for x in xticks])

# Customize the plot
plt.xlabel('Length (kyr)')
plt.ylabel('Percentage (%)')
#plt.title(f'Mean hiatus length 200,000 yr (Bins: {bin_length} yr)')
plt.title(f'Mean hiatus length 200 000 y')
plt.ylim(0, 100)
#plt.grid(axis='y', linestyle='--', alpha=0.5)  # Optional grid

# Show the plot
plt.show()