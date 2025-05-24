import numpy as np
from tqdm import tqdm

def simulate_geomagnetic_reversals(lambda_val, span, n=21, min_gap_yr=3e4):
    min_gap = min_gap_yr / 1e6
    needed  = n
    while True:
        waits = np.random.exponential(scale=1 / lambda_val, size=needed*2)
        waits = waits[waits >= min_gap]
        if waits.size >= n:
            waits = waits[:n]
            break
        needed *= 2
    rev = np.cumsum(waits)
    # scale to [0, span]
    rev = (rev - rev[0]) / (rev[-1] - rev[0]) * span
    # one pass to enforce min_gap after scaling
    diff = np.maximum(min_gap, np.diff(rev, prepend=0))
    rev  = np.cumsum(diff)
    rev  = np.clip(rev, None, span)
    # zones (disjoint) -------------------------------------------------
    half = changing_state_time / 2.0

    start = rev[:-1] + half           # r_i  + ½Δ  (Δ = change-zone width)
    end   = rev[1:]  - half           # r_i+1 − ½Δ
    
    # keep the outermost edges intact
    start[0] = rev[0]
    end[-1]  = rev[-1]
    
    mz = np.column_stack((start, end))
    #mz   = np.column_stack((np.r_[rev[0], rev[:-1] + half],
     #                       np.r_[rev[:-1] - half, rev[-1]]))
    cz   = np.column_stack((rev[1:-1] - half, rev[1:-1] + half))
    return rev, mz, cz


def _merge(intervals: np.ndarray) -> np.ndarray:
    """Union of [start, end) intervals.  `intervals` must be (!,2)."""
    if len(intervals) == 0:
        return intervals
    idx = np.argsort(intervals[:, 0])
    s   = intervals[idx, 0]
    e   = intervals[idx, 1]
    e   = np.maximum.accumulate(e)
    keep = np.empty_like(s, dtype=bool)
    keep[0] = True
    keep[1:] = s[1:] > e[:-1]
    out = np.column_stack((s[keep], e[keep]))
    return out

def simulate_diastem(time_span_myr: float,
                          min_len: float,
                          max_len: float,
                          gap_percent: float,
                          batch_size: int = 1024,
                          rng: np.random.Generator | None = None
                         ) -> np.ndarray:
    """
    Fast & unbiased generator of disjoint diastems whose *total length*
    >= gap_percent * time_span_myr and whose *individual* lengths are
    i.i.d. U(min_len, max_len).

    Statistics match the original slow merge-loop; speed is 10-20×.
    """
    if rng is None:
        rng = np.random.default_rng()

    target = time_span_myr * gap_percent
    covered = 0.0

    diastems = np.empty((0, 2), float)
    gaps     = np.array([[0.0, time_span_myr]], float)       # start with one free patch

    # --- main loop ---------------------------------------------------
    while covered < target and gaps.size:
        # 1 ▸ choose patches for the batch
        sizes   = gaps[:, 1] - gaps[:, 0]
        probs   = sizes / sizes.sum()
        pick_ix = rng.choice(len(gaps), size=batch_size, p=probs)

        # 2 ▸ draw lengths and keep only those that fit
        rooms   = sizes[pick_ix]
        lengths = rng.uniform(min_len, max_len, size=batch_size)
        ok      = lengths <= rooms
        if not ok.any():
            continue                       # rare when only tiny gaps remain

        pick_ix = pick_ix[ok]
        rooms   = rooms[ok]
        lengths = lengths[ok]

        # 3 ▸ draw starts for the accepted intervals
        g0  = gaps[pick_ix, 0]
        starts = rng.uniform(g0, g0 + rooms - lengths)
        ends   = starts + lengths

        new_intervals = np.column_stack((starts, ends))
        diastems = _merge(np.vstack((diastems, new_intervals)))
        covered  = (diastems[:, 1] - diastems[:, 0]).sum()

        # 4 ▸ rebuild the list of free gaps for the next batch
        if covered < target:
            # complement: [0, span] \ big_union
            lefts  = np.r_[0.0, diastems[:, 1]]
            rights = np.r_[diastems[:, 0], time_span_myr]
            mask   = rights > lefts
            gaps   = np.column_stack((lefts[mask], rights[mask]))

    # Optional last tweak: if covered overshot the target, trim the last interval
    if covered > target:
        overshoot = covered - target
        diastems[-1, 1] -= overshoot

    return diastems

'''
def simulate_diastem(time_span_myr, min_gap_length, max_gap_length, gap_percent):
    diastems = []
    target_gap = time_span_myr * gap_percent

    while True:
        # 1. draw a new interval
        length = np.random.uniform(min_gap_length, max_gap_length)
        start  = np.random.uniform(0, time_span_myr - length)
        end    = start + length

        # 2. insert and merge
        diastems.append((start, end))
        diastems.sort()
        merged = []
        for s, e in diastems:
            if merged and s <= merged[-1][1]:          # overlap/adjacent
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        diastems = merged

        # 3. recompute true coverage
        total_gap = sum(e - s for s, e in diastems)
        if total_gap >= target_gap or total_gap >= time_span_myr:
            break

    return np.asarray(diastems, dtype=float)
'''

'''
#Old/ Don't delete'
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
'''


'''
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
'''

def simulate_diastem_poisson(time_span_myr: float,
                             average_diastem_length: float,
                             gap_percent: float,
                             rng: np.random.Generator | None = None,
                             oversample: float = 2.0
                            ) -> np.ndarray:
    """
    Vectorised version of the Poisson diastem simulator.

    Parameters
    ----------
    time_span_myr : float
        Total modelled span (Myr).
    average_diastem_length : float
        Mean of BOTH the inter-diastem waiting time (T1)
        and the diastem duration (T2) in the classic two-stage
        Poisson model (Myr).
    gap_percent : float
        Fraction of the time span that must be covered by diastems
        (0 ≤ gap_percent ≤ 1).
    rng : np.random.Generator, optional
        Pass your own RNG for reproducibility.
    oversample : float, optional
        Safety factor for the batch size estimate. 2.0 is plenty
        unless you are asking for gap_percent → 1.0.

    Returns
    -------
    diastems : np.ndarray, shape (n, 2)
        Sorted, non-overlapping intervals [[start₀, end₀], …].
    """
    if rng is None:
        rng = np.random.default_rng()

    # --- how many intervals do we *expect* to need? ------------------
    target_gap   = time_span_myr * gap_percent
    mean_gap_len = average_diastem_length          # E[T2]
    n_expect     = int(np.ceil(target_gap / mean_gap_len * oversample)) + 8

    while True:                                    # usually runs ONCE
        # draw inter-arrivals and durations in one shot
        inter = rng.exponential(scale=average_diastem_length, size=n_expect)
        dur   = rng.exponential(scale=average_diastem_length, size=n_expect)

        # Poisson process ⇒ cumulative sum gives the start times
        start = np.cumsum(inter)
        keep  = start < time_span_myr              # ignore ones past the end
        if not keep.any():                         # very narrow span
            # fall back to one interval that fills the remainder
            return np.array([[0.0, min(target_gap, time_span_myr)]])
        start = start[keep]
        dur   = dur[keep]

        end   = np.minimum(start + dur, time_span_myr)
        dur   = end - start                       # clip duration

        # cumulative covered time
        covered = np.cumsum(dur)
        hit_idx = np.searchsorted(covered, target_gap, side='right')

        if hit_idx < covered.size:
            # trim to the first interval that *reaches* the target
            start = start[:hit_idx + 1]
            end   = end[:hit_idx + 1]

            # shrink the very last interval so total == target_gap exactly
            overshoot = covered[hit_idx] - target_gap
            if overshoot > 0.0:
                end[-1] -= overshoot              # cannot make it negative

            return np.column_stack((start, end))

        # Not enough coverage ⇒ double the batch and loop once more
        n_expect *= 2


def get_lost_percent(magnetozones, diastems, change_zones,min_remaining_myr):
    
    def _percent_too_short(zones, threshold):
            if threshold is None or len(zones) == 0:
                return 0.0
    
            too_short = 0
            for zs, ze in zones:
                remaining = ze - zs
                for ds, de in diastems:
                    if de <= zs or ds >= ze:          # no overlap
                        continue
                    remaining -= min(ze, de) - max(zs, ds)
                    if remaining <= 0:
                        break
                if remaining < threshold:
                    too_short += 1
            return too_short / len(zones)
        
    
    

    
    lost_total_magnetozone_length = []
    
    for m_start, m_end in magnetozones:
        zone_lost_fractions = []
        
        for d_start, d_end in diastems:
            if d_end < m_start:
                continue
            if d_start > m_end:
                break
            
            zone_lost_fractions.append((min(m_end, d_end) - max(m_start, d_start))/(m_end - m_start))
            
        if not zone_lost_fractions:
            lost_total_magnetozone_length.append(0)
        else:
            lost_total_magnetozone_length.append(sum(zone_lost_fractions))
    
    lost_total_magnetozone_mean_length = np.asarray(lost_total_magnetozone_length)
    
    if lost_total_magnetozone_mean_length.size > 0:
        lost_magnetozone_length = lost_total_magnetozone_mean_length.mean()
    else:
        lost_magnetozone_length = 0
    
    
    total_magnetozones = len(magnetozones)
    
    fully_lost_magnetozones = sum(
        any(d_start <= m_start and d_end >= m_end for d_start, d_end in diastems)for m_start, m_end in magnetozones)
    
    fully_lost_magnetozones_percent = fully_lost_magnetozones / total_magnetozones if total_magnetozones > 0 else 0
    
    
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

    too_short_magnetozones_percent = _percent_too_short(
        magnetozones, min_remaining_myr
    )
    
    too_short_change_zones_percent = _percent_too_short(
        change_zones, min_remaining_myr
    )

    return fully_lost_magnetozones_percent,lost_change_zones,fully_lost_change_zones_percent, too_short_magnetozones_percent, too_short_change_zones_percent, lost_magnetozone_length


'''
def get_lost_percent(magnetozones, diastems, change_zones, min_remaining_myr):
    def _percent_too_short(zones, threshold):
        if threshold is None or len(zones) == 0:
            return 0.0
        
        too_short = 0
        zones_arr = np.asarray(zones)
        diastems_arr = np.asarray(diastems)
        
        for (zs, ze) in zones:
            remaining = ze - zs
            # Vectorized overlap calculation
            overlaps = np.minimum(ze, diastems_arr[:,1]) - np.maximum(zs, diastems_arr[:,0])
            overlaps = overlaps[overlaps > 0]  # Only positive overlaps
            remaining -= np.sum(overlaps)
            if remaining < threshold:
                too_short += 1
        return too_short / len(zones)

    # Convert to numpy arrays for vectorization
    mag_arr = np.asarray(magnetozones)
    chg_arr = np.asarray(change_zones)
    dias_arr = np.asarray(diastems)
    
    # Magnetozone calculations
    mag_lengths = mag_arr[:,1] - mag_arr[:,0]
    mag_overlaps = np.maximum(0, 
        np.minimum(mag_arr[:,1,None], dias_arr[None,:,1]) - 
        np.maximum(mag_arr[:,0,None], dias_arr[None,:,0])
    )
    mag_lost_frac = np.sum(mag_overlaps, axis=1) / mag_lengths
    lost_magnetozone_length = np.mean(mag_lost_frac) if len(mag_lost_frac) > 0 else 0
    
    # Fully lost magnetozones
    fully_lost = np.any(
        (dias_arr[:,0,None] <= mag_arr[None,:,0]) & 
        (dias_arr[:,1,None] >= mag_arr[None,:,1]),
        axis=0
    )
    fully_lost_magnetozones_percent = np.mean(fully_lost) if len(mag_arr) > 0 else 0

    # Change zone calculations (same vectorized approach)
    chg_lengths = chg_arr[:,1] - chg_arr[:,0]
    chg_overlaps = np.maximum(0,
        np.minimum(chg_arr[:,1,None], dias_arr[None,:,1]) - 
        np.maximum(chg_arr[:,0,None], dias_arr[None,:,0])
    )
    chg_lost_frac = np.sum(chg_overlaps, axis=1) / chg_lengths
    lost_change_zones = np.mean(chg_lost_frac) if len(chg_lost_frac) > 0 else 0
    
    # Fully lost change zones
    fully_lost_chg = np.any(
        (dias_arr[:,0,None] <= chg_arr[None,:,0]) & 
        (dias_arr[:,1,None] >= chg_arr[None,:,1]),
        axis=0
    )
    fully_lost_change_zones_percent = np.mean(fully_lost_chg) if len(chg_arr) > 0 else 0

    # Too short percentages
    too_short_magnetozones_percent = _percent_too_short(magnetozones, min_remaining_myr)
    too_short_change_zones_percent = _percent_too_short(change_zones, min_remaining_myr)

    return (
        fully_lost_magnetozones_percent,
        lost_change_zones,
        fully_lost_change_zones_percent, 
        too_short_magnetozones_percent, 
        too_short_change_zones_percent, 
        lost_magnetozone_length
    )
 

'''


def _vectorised_overlap(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Return, for every interval in *a*, the total overlap length with *all*
    intervals in *b*.  Overlaps are summed even when *b* intervals overlap
    each other, reproducing the behaviour of the original loops.
    """
    if a.size == 0 or b.size == 0:
        return np.zeros(len(a))

    a_start = a[:, 0][:, None]
    a_end   = a[:, 1][:, None]
    b_start = b[:, 0][None, :]
    b_end   = b[:, 1][None, :]

    overlaps = np.clip(
        np.minimum(a_end, b_end) - np.maximum(a_start, b_start),
        0.0, None
    )
    return overlaps.sum(axis=1)


def _contained(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Boolean mask: True where the *a* interval is fully inside any *b*."""
    if a.size == 0 or b.size == 0:
        return np.zeros(len(a), dtype=bool)

    a_start = a[:, 0][:, None]
    a_end   = a[:, 1][:, None]
    b_start = b[:, 0][None, :]
    b_end   = b[:, 1][None, :]

    return ((b_start <= a_start) & (b_end >= a_end)).any(axis=1)


def _too_short(overlap: np.ndarray,
               lengths: np.ndarray,
               threshold: float | None) -> float:
    """Fraction of zones whose remaining length falls below *threshold*."""
    if threshold is None or lengths.size == 0:
        return 0.0
    remaining = lengths - overlap
    return (remaining < threshold).mean()

'''
# ──────────────────────────────────────────────────────────────────────────────
# Public vectorised replacement
# ──────────────────────────────────────────────────────────────────────────────
def get_lost_percent(
    magnetozones: list[tuple[float, float]],
    diastems:      list[tuple[float, float]],
    change_zones:  list[tuple[float, float]],
    min_remaining_myr: float | None,
):
    """
    Vectorised replacement for the original get_lost_percent.

    Returns:
        fully_lost_magnetozones_percent,
        lost_change_zones,
        fully_lost_change_zones_percent,
        too_short_magnetozones_percent,
        too_short_change_zones_percent,
        lost_magnetozone_length
    """
    # --- convert inputs to (N, 2) float64 arrays --------------------------------
    mz = np.asarray(magnetozones, dtype=np.float64)
    dz = np.asarray(diastems,      dtype=np.float64)
    cz = np.asarray(change_zones,  dtype=np.float64)

    # --- magnetozones -----------------------------------------------------------
    mz_lengths  = mz[:, 1] - mz[:, 0] if mz.size else np.array([])
    mz_overlap  = _vectorised_overlap(mz, dz)
    lost_mz_frac = (
        mz_overlap / mz_lengths if mz_lengths.size else np.array([])
    )
    lost_magnetozone_length = (
        float(lost_mz_frac.mean()) if lost_mz_frac.size else 0.0
    )
    fully_lost_mz_percent = (
        float(_contained(mz, dz).mean()) if mz.size else 0.0
    )
    too_short_mz_percent = _too_short(
        mz_overlap, mz_lengths, min_remaining_myr
    )

    # --- change-zones -----------------------------------------------------------
    cz_lengths  = cz[:, 1] - cz[:, 0] if cz.size else np.array([])
    cz_overlap  = _vectorised_overlap(cz, dz)
    lost_cz_frac = (
        cz_overlap / cz_lengths if cz_lengths.size else np.array([])
    )
    lost_change_zones = (
        float(lost_cz_frac.mean()) if lost_cz_frac.size else 0.0
    )
    fully_lost_cz_percent = (
        float(_contained(cz, dz).mean()) if cz.size else 0.0
    )
    too_short_cz_percent = _too_short(
        cz_overlap, cz_lengths, min_remaining_myr
    )

    # --- output -----------------------------------------------------------------
    return (
        fully_lost_mz_percent,
        lost_change_zones,
        fully_lost_cz_percent,
        too_short_mz_percent,
        too_short_cz_percent,
        lost_magnetozone_length,
    )

'''
         
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
            
def iter(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,min_gap_length,max_gap_length,gap_percent,reversal_number,iterations_number,min_remaining_myr):
    
    lost_magnetozones_list = []
    lost_magnetozones_length_list = []
    lost_change_zones_list = []
    fully_lost_change_zones__list = []
    
    too_short_magnetozones_percent_list = []
    too_short_change_zones_percent_list = []

    
    for i in tqdm(range(iterations_number), desc=f"Running Simulation for {gap_percent*100}%"):

        reversal_times, magnetozones, change_zones = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, reversal_number, min_gap_years)
        diastems = simulate_diastem(time_span_myr, min_gap_length, max_gap_length, gap_percent)

        lost_magnetozones,lost_change_zones,fully_lost_change_zones_percent, too_short_magnetozones_percent, too_short_change_zones_percent, lost_magnetozones_length = get_lost_percent(magnetozones, diastems, change_zones, min_remaining_myr)
        
        lost_magnetozones_list.append(lost_magnetozones)
        lost_change_zones_list.append(lost_change_zones)
        lost_magnetozones_length_list.append(lost_magnetozones_length)
        fully_lost_change_zones__list.append(fully_lost_change_zones_percent)
        too_short_magnetozones_percent_list.append(too_short_magnetozones_percent)
        too_short_change_zones_percent_list.append(too_short_change_zones_percent)
        
    lost_magnetozones_list = np.asarray(lost_magnetozones_list)
    lost_change_zones_list = np.asarray(lost_change_zones_list)
    fully_lost_change_zones__list = np.asarray(fully_lost_change_zones__list)
    too_short_change_zones_percent_list = np.asarray(too_short_change_zones_percent_list)
    too_short_magnetozones_percent_list = np.asarray(too_short_magnetozones_percent_list)
    lost_magnetozones_length_list = np.asarray(lost_magnetozones_length_list)

    lost_magnetozones_list.sort()
    lost_change_zones_list.sort()
    fully_lost_change_zones__list.sort()
    too_short_change_zones_percent_list.sort()
    too_short_magnetozones_percent_list.sort()
    lost_magnetozones_length_list.sort()

    magnetozones_lower_bound, magnetozones_upper_bound = np.percentile(lost_magnetozones_list, [2.5, 97.5])
    change_zones_lower_bound, change_zones_upper_bound = np.percentile(lost_change_zones_list, [2.5, 97.5])
    fully_lost_change_zones_lower_bound, fully_lost_change_zones_upper_bound = np.percentile(fully_lost_change_zones__list, [2.5, 97.5])
    too_short_change_zones_lower_bound, too_short_change_zones_upper_bound = np.percentile(too_short_change_zones_percent_list, [2.5, 97.5])
    too_short_magnetozones_lower_bound, too_short_magnetozones_upper_bound = np.percentile(too_short_magnetozones_percent_list, [2.5, 97.5])
    lost_magnetozones_length_lower_bound, lost_magnetozones_length_upper_bound = np.percentile(lost_magnetozones_length_list, [2.5, 97.5])
    
    summary_data = [
        "Lost Magnetozones:",
        f"95% Confidence Interval: {magnetozones_lower_bound* 100:.4f}, {magnetozones_upper_bound* 100:.4f}",
        f"Mean percent of lost Magnetozones: {lost_magnetozones_list.mean() * 100:.4f}",
        "Practically lost Magnetozones:",
        f"95% Confidence Interval: {too_short_magnetozones_lower_bound* 100:.4f}, {too_short_magnetozones_upper_bound* 100:.4f}",
        f"Mean percent of Practically lost Magnetozones: {too_short_magnetozones_percent_list.mean() * 100:.4f}",
        "Lost Magnetozones thickness:",
        f"95% Confidence Interval: {lost_magnetozones_length_lower_bound* 100:.4f}, {lost_magnetozones_length_upper_bound* 100:.4f}",
        f"Mean percent of lost Magnetozones thickness: {lost_magnetozones_length_list.mean() * 100:.4f}",
        "Lost change zones thickness:",
        f"95% Confidence Interval: {change_zones_lower_bound* 100:.4f}, {change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones thickness: {lost_change_zones_list.mean() * 100:.4f}",
        "Fully lost change zones:",
        f"95% Confidence Interval: {fully_lost_change_zones_lower_bound* 100:.4f}, {fully_lost_change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones: {fully_lost_change_zones__list.mean() * 100:.4f}",
        "Practically lost change zones:",
        f"95% Confidence Interval: {too_short_change_zones_lower_bound* 100:.4f}, {too_short_change_zones_upper_bound* 100:.4f}",
        f"Mean percent of Practically lost zones: {too_short_change_zones_percent_list.mean() * 100:.4f}"
    ]

    header_data = [
        f"Time span: {time_span_myr} MYr",
        f"Number of reversals: {reversal_number - 2}",
        f"Duration of the transition state: {changing_state_time*1e6} yr",
        f"Diastem coverage: {gap_percent*100}%",
        f"Diastem length interval {min_gap_length*1e6} yr, {max_gap_length*1e6} yr",    
        f"Iterations: {iterations_number}",
        f"Min remainig year: {min_remaining_myr*1e6} yr"
    ]

    save_to_file(f"Time_span_{time_span_myr}myr gap_percent_{gap_percent} reversals_{reversal_number - 2} iterations_{iterations_number} max gap_{max_gap_length} min rem yr {min_remaining_myr*1e6}.txt",
                 summary_data, header="\n".join(header_data) + "\n\n")
    
def iterPoisson(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number,min_remaining_myr):
    
    lost_magnetozones_list = []
    lost_magnetozones_length_list = []
    lost_change_zones_list = []
    fully_lost_change_zones__list = []
    
    too_short_magnetozones_percent_list = []
    too_short_change_zones_percent_list = []

    
    for i in tqdm(range(iterations_number), desc=f"Running Simulation for {gap_percent*100}%"):

        reversal_times, magnetozones, change_zones = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, reversal_number, min_gap_years)
        diastems = simulate_diastem_poisson(time_span_myr, average_diastem_length, gap_percent)

        lost_magnetozones,lost_change_zones,fully_lost_change_zones_percent, too_short_magnetozones_percent, too_short_change_zones_percent,lost_magnetozones_length = get_lost_percent(magnetozones, diastems, change_zones, min_remaining_myr)
        
        lost_magnetozones_list.append(lost_magnetozones)
        lost_change_zones_list.append(lost_change_zones)
        lost_magnetozones_length_list.append(lost_magnetozones_length)
        fully_lost_change_zones__list.append(fully_lost_change_zones_percent)
        too_short_magnetozones_percent_list.append(too_short_magnetozones_percent)
        too_short_change_zones_percent_list.append(too_short_change_zones_percent)
        
    lost_magnetozones_list = np.asarray(lost_magnetozones_list)
    lost_change_zones_list = np.asarray(lost_change_zones_list)
    fully_lost_change_zones__list = np.asarray(fully_lost_change_zones__list)
    too_short_change_zones_percent_list = np.asarray(too_short_change_zones_percent_list)
    too_short_magnetozones_percent_list = np.asarray(too_short_magnetozones_percent_list)
    lost_magnetozones_length_list = np.asarray(lost_magnetozones_length_list)

    lost_magnetozones_list.sort()
    lost_change_zones_list.sort()
    fully_lost_change_zones__list.sort()
    too_short_change_zones_percent_list.sort()
    too_short_magnetozones_percent_list.sort()
    lost_magnetozones_length_list.sort()

    magnetozones_lower_bound, magnetozones_upper_bound = np.percentile(lost_magnetozones_list, [2.5, 97.5])
    change_zones_lower_bound, change_zones_upper_bound = np.percentile(lost_change_zones_list, [2.5, 97.5])
    fully_lost_change_zones_lower_bound, fully_lost_change_zones_upper_bound = np.percentile(fully_lost_change_zones__list, [2.5, 97.5])
    too_short_change_zones_lower_bound, too_short_change_zones_upper_bound = np.percentile(too_short_change_zones_percent_list, [2.5, 97.5])
    too_short_magnetozones_lower_bound, too_short_magnetozones_upper_bound = np.percentile(too_short_magnetozones_percent_list, [2.5, 97.5])
    lost_magnetozones_length_lower_bound, lost_magnetozones_length_upper_bound = np.percentile(lost_magnetozones_length_list, [2.5, 97.5])
    
    summary_data = [
        "Lost Magnetozones:",
        f"95% Confidence Interval: {magnetozones_lower_bound* 100:.4f}, {magnetozones_upper_bound* 100:.4f}",
        f"Mean percent of lost Magnetozones: {lost_magnetozones_list.mean() * 100:.4f}",
        "Practically lost Magnetozones:",
        f"95% Confidence Interval: {too_short_magnetozones_lower_bound* 100:.4f}, {too_short_magnetozones_upper_bound* 100:.4f}",
        f"Mean percent of Practically lost Magnetozones: {too_short_magnetozones_percent_list.mean() * 100:.4f}",
        "Lost Magnetozones thickness:",
        f"95% Confidence Interval: {lost_magnetozones_length_lower_bound* 100:.4f}, {lost_magnetozones_length_upper_bound* 100:.4f}",
        f"Mean percent of lost Magnetozones thickness: {lost_magnetozones_length_list.mean() * 100:.4f}",
        "Lost change zones thickness:",
        f"95% Confidence Interval: {change_zones_lower_bound* 100:.4f}, {change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones thickness: {lost_change_zones_list.mean() * 100:.4f}",
        "Fully lost change zones:",
        f"95% Confidence Interval: {fully_lost_change_zones_lower_bound* 100:.4f}, {fully_lost_change_zones_upper_bound* 100:.4f}",
        f"Mean percent of lost zones: {fully_lost_change_zones__list.mean() * 100:.4f}",
        "Practically lost change zones:",
        f"95% Confidence Interval: {too_short_change_zones_lower_bound* 100:.4f}, {too_short_change_zones_upper_bound* 100:.4f}",
        f"Mean percent of Practically lost zones: {too_short_change_zones_percent_list.mean() * 100:.4f}"
    ]

    header_data = [
        f"Time span: {time_span_myr} MYr",
        f"Number of reversals: {reversal_number - 2}",
        f"Duration of the transition state: {changing_state_time*1e6} yr",
        f"Diastem coverage: {gap_percent*100}%",
        f"Average diastem length  {average_diastem_length*1e6} yr",    
        f"Iterations: {iterations_number}",
        f"Min remainig year: {min_remaining_myr*1e6} yr"
    ]

    save_to_file(f"Poisson_Time_span_{time_span_myr}myr gap_percent_{gap_percent} reversals_{reversal_number - 2} iterations_{iterations_number} av gap_{average_diastem_length} min rem yr {min_remaining_myr*1e6}.txt",
                 summary_data, header="\n".join(header_data) + "\n\n")

time_span_myr = 5  # Total time in million years

mean_reversal_rate = 0.45

min_gap_years = 30000  # Minimum gap between reversals in years
changing_state_time = 10000  # Time in years the field is in an intermediate state

min_gap_length = 0
max_gap_length = 1000

average_diastem_length=0.0005

#gap_percent = 20

reversal_number = 102

min_remaining_myr = 100 #in yr


min_gap_myr = min_gap_years / 1e6  # Convert years to million years
changing_state_time = changing_state_time/ 1e6  # Convert to million years
min_gap_length = min_gap_length/1e6
max_gap_length = max_gap_length/1e6
#min_remaining_myr = min_remaining_myr/1e6
#gap_percent = gap_percent/100

iterations_number = 1000

gap_percent_list = [10,20,30,40,50,60,70,80,90]

#got to start

#gap_percent_list = [90]

#gap_percent = 90



for min_remaining_myr in [100,20]:
    min_remaining_myr = min_remaining_myr/1e6
    for reversal_number in [102,22]:
        for gap_percent in gap_percent_list:    
            gap_percent = gap_percent/100            
            iter(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,min_gap_length,max_gap_length,gap_percent,reversal_number,iterations_number,min_remaining_myr)
            #iterPoisson(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number,min_remaining_myr)
            #iterPoisson(time_span_myr,alpha, beta, loc,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number)



'''
reversal_times, magnetozones, change_zones = simulate_geomagnetic_reversals(mean_reversal_rate, time_span_myr, reversal_number, min_gap_years)
'''       
'''                                                                                                                                                                                                                 

for min_remaining_myr in [20]:
    min_remaining_myr = min_remaining_myr/1e6
    for gap_percent in gap_percent_list:    
        gap_percent = gap_percent/100            
        #iter(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,min_gap_length,max_gap_length,gap_percent,reversal_number,iterations_number,min_remaining_myr)
        iterPoisson(time_span_myr,mean_reversal_rate,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number,min_remaining_myr)
        #iterPoisson(time_span_myr,alpha, beta, loc,min_gap_years,changing_state_time,average_diastem_length,gap_percent,reversal_number,iterations_number)
  

    
'''