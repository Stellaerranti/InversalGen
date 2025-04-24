import numpy as np
from tqdm import tqdm
from scipy import stats

def simulate_geomagnetic_reversals(lambda_value, time_span_myr, reversal_number=21, min_gap_years=30000):
    """
    Simulate geomagnetic reversals using an Exponential distribution for inter-arrival times.

    Parameters:
        lambda_value (float): Rate parameter λ of the Exponential distribution.
        time_span_myr (float): Total time span in million years.
        reversal_number (int): Target number of reversals.
        min_gap_years (float): Minimum time between reversals (in years).

    Returns:
        reversal_times (np.ndarray): Times of reversals (Myr).
        magnetozones (np.ndarray): Time intervals of stable polarity.
        change_zones (np.ndarray): Time intervals of transitional states.
    """
    reversal_times = []
    current_time = 0
    min_gap_myr = min_gap_years / 1e6  # Convert years to Myr

    while len(reversal_times) < reversal_number:
        # Sample inter-arrival time from Exponential distribution
        wait_time = stats.expon.rvs(scale=1 / lambda_value)

        if wait_time < min_gap_myr:
            continue

        current_time += wait_time
        reversal_times.append(current_time)

    reversal_times = np.array(reversal_times)
    reversal_times = (reversal_times - reversal_times.min()) / (reversal_times.max() - reversal_times.min()) * time_span_myr

    for i in range(1, len(reversal_times)):
        if reversal_times[i] - reversal_times[i - 1] < min_gap_myr:
            reversal_times[i] = reversal_times[i - 1] + min_gap_myr
        if reversal_times[i] > time_span_myr:
            reversal_times[i] = time_span_myr
            reversal_times = reversal_times[:i + 1]
            break

    magnetozones = np.column_stack((reversal_times[:-1], reversal_times[1:]))
    change_zones = np.column_stack((
        reversal_times[1:-1] - changing_state_time / 2,
        reversal_times[1:-1] + changing_state_time / 2
    ))

    return reversal_times, magnetozones, change_zones

def simulate_diastem(time_span_myr, min_gap_length, max_gap_length, gap_percent):
    diastems = []
    total_gap = 0
    target_gap = time_span_myr * gap_percent
    
    while total_gap < target_gap:
        length = np.random.uniform(min_gap_length, max_gap_length)
        start = np.random.uniform(0, time_span_myr - length)
        end = start + length
        
        # Merge overlapping intervals 
        new_diastem = (start, end)
        updated_diastems = []
        added = False
        for s, e in diastems:
            if e < start or s > end:  # No overlap
                updated_diastems.append((s, e))
            else:  # Overlapping case
                new_start = min(s, start)
                new_end = max(e, end)
                new_length = new_end - new_start
                total_gap += new_length - (e - s)  # Adjust total gap length
                new_diastem = (new_start, new_end)
                added = True
        
        updated_diastems.append(new_diastem)
        diastems = sorted(updated_diastems, key=lambda x: x[0])
        
        if not added:
            total_gap += length
    
    return np.array(diastems)

def simulate_diastem_poisson(time_span_myr, average_diastem_length, gap_percent):
    """
    Simulate diastems (hiatuses) using a Poisson process for both their occurrence and duration.
    Parameters:
        time_span_myr (float): Total modeled time span (Myr).
        average_diastem_length (float): Mean diastem duration T2 (Myr).
        gap_percent (float): Total % of time to be covered by diastems.
    Returns:
        np.ndarray: Array of (start, end) tuples for each diastem.
    """
    diastems = []
    total_gap = 0.0
    target_gap = time_span_myr * gap_percent

    current_time = 0.0

    while total_gap < target_gap and current_time < time_span_myr:
        # Wait until the next diastem begins (interarrival time)
        interarrival = np.random.exponential(scale=average_diastem_length)
        current_time += interarrival

        if current_time >= time_span_myr:
            break

        # Diastem duration
        duration = np.random.exponential(scale=average_diastem_length)

        end_time = current_time + duration
        if end_time > time_span_myr:
            end_time = time_span_myr
            duration = end_time - current_time

        diastems.append((current_time, end_time))
        total_gap += duration

        # Move time to end of this diastem (not strictly necessary, but avoids overlap)
        current_time = end_time

    return np.array(diastems)

def get_lost_percent(magnetozones, diastems, change_zones):
    
    '''
    total_magnetozone_length = np.sum(magnetozones[:, 1] - magnetozones[:, 0])
    lost_length = 0
    '''
   
    total_magnetozones = len(magnetozones)
    
    fully_lost_magnetozones = sum(
        any(d_start <= m_start and d_end >= m_end for d_start, d_end in diastems)for m_start, m_end in magnetozones)
    
    fully_lost_magnetozones_percent = fully_lost_magnetozones / total_magnetozones if total_magnetozones > 0 else 0
    
    '''
    
    for m_start, m_end in magnetozones:
        for d_start, d_end in diastems:
            if d_end < m_start:
                continue
            if d_start > m_end:
                break
            lost_length += min(m_end, d_end) - max(m_start, d_start)
            
    lost_magnetozones = lost_length / total_magnetozone_length if total_magnetozone_length > 0 else 0
    '''
    
    '''
    total_change_zones_length = np.sum(change_zones[:, 1] - change_zones[:, 0])
    lost_change_zones_length = 0
    
    for c_start, c_end in change_zones:
        for d_start, d_end in diastems:
            if d_end < c_start:
                continue
            if d_start > c_end:
                break
            lost_change_zones_length += min(c_end, d_end) - max(c_start, d_start)
            
    lost_change_zones = lost_change_zones_length / total_change_zones_length if total_change_zones_length > 0 else 0
    
    '''
    
    lost_change_zones_mean_length = []
    
    for c_start, c_end in change_zones:
        zone_lost_fractions = []
    
        for d_start, d_end in diastems:
            if d_end < c_start:
                continue
            if d_start > c_end:
                break
    
            zone_lost_fractions.append((min(c_end, d_end) - max(c_start, d_start)) / (c_end - c_start))
    
        # If the zone was untouched, append 0
        if not zone_lost_fractions:
            lost_change_zones_mean_length.append(0)
        else:
            lost_change_zones_mean_length.append(sum(zone_lost_fractions))
            
              
    lost_change_zones_mean_length = np.asarray(lost_change_zones_mean_length)
    
    if lost_change_zones_mean_length.size > 0:
        lost_change_zones = lost_change_zones_mean_length.mean()
    else:
        lost_change_zones = 0
        
    total_change_zones = len(change_zones)
    fully_lost_change_zones = sum(
        any(d_start <= c_start and d_end >= c_end for d_start, d_end in diastems)
        for c_start, c_end in change_zones
    )
    
    fully_lost_change_zones_percent = fully_lost_change_zones / total_change_zones if total_change_zones > 0 else 0
    
    return fully_lost_magnetozones_percent,lost_change_zones,fully_lost_change_zones_percent
            
def save_to_file(filename, data, header=None):
    with open(filename, 'w') as f:
        if header:
            f.write(header)
        if isinstance(data, (list, np.ndarray)):
            for item in data:
                if isinstance(item, (list, tuple, np.ndarray)):
                    f.write(" ".join(map(str, item)) + "\n")
                else:
                    f.write(str(item) + "\n")
        else:
            f.write(str(data) + "\n")
            
def iter(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,min_gap_length,max_gap_length,gap_percent,reversal_number,iterations_number):
    
    lost_magnetozones_list = []
    lost_change_zones_list = []
    fully_lost_change_zones__list = []

    
    for i in tqdm(range(iterations_number), desc=f"Running Simulation for {gap_percent*100}%"):

        reversal_times, magnetozones, change_zones = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, reversal_number, min_gap_years)
        diastems = simulate_diastem(time_span_myr, min_gap_length, max_gap_length, gap_percent)

        lost_magnetozones,lost_change_zones,fully_lost_change_zones_percent = get_lost_percent(magnetozones, diastems, change_zones)
        
        lost_magnetozones_list.append(lost_magnetozones)
        lost_change_zones_list.append(lost_change_zones)
        fully_lost_change_zones__list.append(fully_lost_change_zones_percent)
        
    lost_magnetozones_list = np.asarray(lost_magnetozones_list)
    lost_change_zones_list = np.asarray(lost_change_zones_list)
    fully_lost_change_zones__list = np.asarray(fully_lost_change_zones__list)

    lost_magnetozones_list.sort()
    lost_change_zones_list.sort()
    fully_lost_change_zones__list.sort()

    magnetozones_lower_bound, magnetozones_upper_bound = np.percentile(lost_magnetozones_list, [2.5, 97.5])
    change_zones_lower_bound, change_zones_upper_bound = np.percentile(lost_change_zones_list, [2.5, 97.5])
    fully_lost_change_zones_lower_bound, fully_lost_change_zones_upper_bound = np.percentile(fully_lost_change_zones__list, [2.5, 97.5])
    
    
    summary_data = [
        "Lost Magnetozones thickness:",
        f"95% Confidence Interval: {magnetozones_lower_bound* 100:.4f}, {magnetozones_upper_bound* 100:.4f}",
        f"Mean percent of lost thickness: {lost_magnetozones_list.mean() * 100:.4f}",
        "Lost change zones thickness:",
        f"95% Confidence Interval: {change_zones_lower_bound* 100:.4f}, {change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones thickness: {lost_change_zones_list.mean() * 100:.4f}",
        "Fully lost change zones:",
        f"95% Confidence Interval: {fully_lost_change_zones_lower_bound* 100:.4f}, {fully_lost_change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones: {fully_lost_change_zones__list.mean() * 100:.4f}"
    ]

    header_data = [
        f"Time span: {time_span_myr} MYr",
        f"Number of reversals: {reversal_number - 2}",
        f"Duration of the transition state: {changing_state_time*1e6} yr",
        f"Diastem coverage: {gap_percent*100}%",
        f"Diastem length interval {min_gap_length*1e6} yr, {max_gap_length*1e6} yr",    
        f"Iterations: {iterations_number}"
    ]

    save_to_file(f"Time_span_{time_span_myr}myr gap_percent_{gap_percent} reversals_{reversal_number - 2} iterations_{iterations_number} max gap_{max_gap_length}.txt",
                 summary_data, header="\n".join(header_data) + "\n\n")
    
def iterPoisson(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number):
    lost_magnetozones_list = []
    lost_change_zones_list = []
    fully_lost_change_zones__list = []

    
    for i in tqdm(range(iterations_number), desc=f"Running Simulation for {gap_percent*100}%"):

        reversal_times, magnetozones, change_zones = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, reversal_number, min_gap_years)
        diastems = simulate_diastem_poisson(time_span_myr, average_diastem_length, gap_percent)

        lost_magnetozones,lost_change_zones,fully_lost_change_zones_percent = get_lost_percent(magnetozones, diastems, change_zones)
        
        lost_magnetozones_list.append(lost_magnetozones)
        lost_change_zones_list.append(lost_change_zones)
        fully_lost_change_zones__list.append(fully_lost_change_zones_percent)
        
    lost_magnetozones_list = np.asarray(lost_magnetozones_list)
    lost_change_zones_list = np.asarray(lost_change_zones_list)
    fully_lost_change_zones__list = np.asarray(fully_lost_change_zones__list)

    lost_magnetozones_list.sort()
    lost_change_zones_list.sort()
    fully_lost_change_zones__list.sort()

    magnetozones_lower_bound, magnetozones_upper_bound = np.percentile(lost_magnetozones_list, [2.5, 97.5])
    change_zones_lower_bound, change_zones_upper_bound = np.percentile(lost_change_zones_list, [2.5, 97.5])
    fully_lost_change_zones_lower_bound, fully_lost_change_zones_upper_bound = np.percentile(fully_lost_change_zones__list, [2.5, 97.5])
    
    
    summary_data = [
        "Lost Magnetozones thickness:",
        f"95% Confidence Interval: {magnetozones_lower_bound* 100:.4f}, {magnetozones_upper_bound* 100:.4f}",
        f"Mean percent of lost thickness: {lost_magnetozones_list.mean() * 100:.4f}",
        "Lost change zones thickness:",
        f"95% Confidence Interval: {change_zones_lower_bound* 100:.4f}, {change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones thickness: {lost_change_zones_list.mean() * 100:.4f}",
        "Fully lost change zones:",
        f"95% Confidence Interval: {fully_lost_change_zones_lower_bound* 100:.4f}, {fully_lost_change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones: {fully_lost_change_zones__list.mean() * 100:.4f}"
    ]

    header_data = [
        f"Time span: {time_span_myr} MYr",
        f"Number of reversals: {reversal_number - 2}",
        f"Duration of the transition state: {changing_state_time*1e6} yr",
        f"Diastem coverage: {gap_percent*100}%",
        f"Average diastem length, {average_diastem_length*1e6} yr",    
        f"Iterations: {iterations_number}"
    ]

    save_to_file(f"Poisson_process Time_span_{time_span_myr}myr gap_percent_{gap_percent} reversals_{reversal_number - 2} iterations_{iterations_number} max gap_{max_gap_length}.txt",
                 summary_data, header="\n".join(header_data) + "\n\n")

time_span_myr = 5  # Total time in million years

mean_reversal_rate = 0.45

min_gap_years = 30000  # Minimum gap between reversals in years
changing_state_time = 10000  # Time in years the field is in an intermediate state

min_gap_length = 0
max_gap_length = 1000

average_diastem_length=1

#gap_percent = 20

reversal_number = 22


min_gap_myr = min_gap_years / 1e6  # Convert years to million years
changing_state_time = changing_state_time/ 1e6  # Convert to million years
min_gap_length = min_gap_length/1e6
max_gap_length = max_gap_length/1e6
#gap_percent = gap_percent/100

iterations_number = 1000

gap_percent_list = [10,20,30,40,50,60,70,80,90]
#gap_percent_list = [90]
'''
reversal_times, magnetozones, change_zones = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, reversal_number, min_gap_years)
'''

for gap_percent in gap_percent_list:
    gap_percent = gap_percent/100
    iter(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,min_gap_length,max_gap_length,gap_percent,reversal_number,iterations_number)
    #iterPoisson(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number)
    #iterPoisson(time_span_myr,alpha, beta, loc,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number)

