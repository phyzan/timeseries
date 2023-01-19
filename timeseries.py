def Psearch(events, p, pdot, bins):
    '''
    The next lines determines whether the search is about period, pdot, or full 2D search, depending on the format of the arguments 'p' and 'pdot' given
    '''
    p_search, pdot_search, full_search = False, False, False
    if type(pdot) is not list:
        p_search, pdot = True, [pdot]
    elif type(p) is not list:
        pdot_search, p = True, [p]
    else:
        full_search = True
    N_p, N_dp, t0 = len(p), len(pdot), events[0]
    '''If period search, chi_2 is a N_p by 1 matrix. If pdot search, its a 1 by N_dp matrix. If full search its N_p by N_dp'''
    chi_2 = [[0 for _ in range(N_dp)] for _ in range(N_p)] 
    '''
    The algorithm below searches for the (p, pdot) pair that corresponds to the highest chi_squared, or equivalently,
    that when folded, the pulse intensity is the gratest possible
    '''
    I, J, chi2_max = 0, 0, 0
    for i in range(N_p):
        for j in range(N_dp):
            phase = []
            for t in events:
                '''the period change due to its derivative for its time of the events'''
                p_ = p[i] + pdot[j]*(t - t0)/2
                '''translate time to phase'''
                phase.append(((t - t0) % p_)/p_)
            _, counts, _ = hist(phase, bins)
            '''equivalently linear fit f=const to the data (phase,counts)'''
            c = sum(counts)/len(counts)
            '''calculate chi_squared for this specific (p, pdot)'''
            chi_2_ij = sum([(i - c)**2/c for i in counts])
            '''append it to the chi_squared 2d array'''
            chi_2[i][j] = chi_2_ij
            if chi_2_ij > chi2_max:
                I, J, chi2_max = i, j, chi_2_ij
    '''
    I and J are such that chi_2[I][J] = max(chi_2), which also (obviously) means that best period = p[I], best pdot = pdot[J]
    '''
    if p_search:
        return p[I], [chi_2[i][0] for i in range(N_p)], I
    elif pdot_search:
        return pdot[J], chi_2[0], J
    elif full_search:
        return p[I], pdot[J], chi_2, I, J

def efold(events, bins, period, pdot = 0, times = 1):
    '''
    period = period of the lightcurve at the beginning of the events
    pdot = period derivative (constant) throughout all the events of the lightcurve
    Given these, convert time of the events to phase
    '''
    phase = []
    t0 = events[0]
    for t in events:
        p_ = period + pdot*(t - t0)/2
        phase.append(((t - t0) % p_)/p_)
    '''
    split the events in bins, and calculate the counts in each bin
    '''
    phase, counts, binsize = hist(phase, bins)

    '''
    Shift the initial phase, so that min(counts) is at phase = 0
    '''
    min_phase = phase[counts.index(min(counts))] # these 2 lines are to "move" the pulse (change the initial phase) 
    phase = [(1 + i - min_phase) % 1 if i != min_phase else 0 for i in phase] #so that min(counts) is at phase = 0, nothing special
    '''
    Calculate uncertainties: x_err = +-binsize/2, y_error = sqrt(counts) in each bin
    '''
    x_err, y_err = [binsize/2 for _ in phase], [i**0.5 for i in counts]
    '''
    Normalize counts and y_err in each bin, then sort all lists
    '''
    N = sum(counts)/len(phase)
    norm_intensity, y_err = [i/N for i in counts], [i/N for i in y_err]
    phase, norm_intensity, y_err = sort(phase, norm_intensity, y_err)
    phase_repeated = []
    '''
    Expand the pulse beyond phase = 1, if times > 1
    '''
    for i in range(times):
        phase_repeated += [j + i*1 for j in phase]
    norm_intensity *= times
    x_err *= times
    y_err *= times
    return phase_repeated, norm_intensity, x_err, y_err

def hist(events, bins):
    '''
    Given a series of events (e.g if photons hit the detector at t=1, t=1.4, t=3, t=3.1, then events = [1, 1.4, 3, 3.1])
    returns a histogram (bin_list, counts) where counts is a list of total events in each bin. Bins are split so that they represent equal timelength
    Also works when events are given in phase (e.g events = [0.1, 0.14, 0.32, 0.47, 0.78, 0.91, 0.98]), the procedure is the same.
    '''
    counts = [0 for _ in range(bins)]
    min_, max_ = min(events), max(events)
    binsize = (max_ - min_)/bins
    bin_list = [min_ + i*binsize for i in range(bins)]
    for i in events:
        if i != max_:
            index = int((i - min_)/(max_ - min_) * bins)
        else:
            index = -1
        counts[index] += 1
    return bin_list, counts, binsize

def array(start, end, n=None, h=None):
    '''creates a list of elements with start and end, depending on the given step or given number of elements'''
    if n is not None:
        h = (end - start)/(n-1)
    elif h is not None:
        n = int((end-start)/h)+1
        h = (end - start)/(n-1)
    return [start + i*h for i in range(n)]

def fwhm(x, y):
    '''
    Finds full width at half maximum
    '''

    '''(x, y) must approach gaussian for this function to make sense'''
    assert len(x) == len(y)
    half_max = max(y)/2
    I = y.index(max(y))
    i = I
    while y[i] > half_max:
        i -= 1
    x_left = x[i+1]
    i = I
    while y[i] > half_max:
        i += 1
    x_right = x[i-1]
    return (x_right - x_left)/2

def split_events(events, num_of_parts):
    ''' splits time events in parts of equal events (but not necessarily of equal time length) '''
    n = len(events)
    part = []
    for i in range(num_of_parts):
        part.append(events[int(i*n/num_of_parts):int((i+1)*n/num_of_parts)])
    return part

def sort(*args):
    ''' Sort a series of lists, sorting the first and rearranging the other elements according to the first "pilot" list'''
    ''' e.g: ([3, 1, 5, 4], [9, 0, 1, 2]) returns ([1, 3, 4, 5], [0, 9, 2, 1]'''
    assert len({len(i) for i in args}) == 1
    A = [[j for j in i] for i in args]
    pil = A[0]
    n = len(pil)
    N = len(A)
    if n == 1:
        return A
    elif n == 2:
        i, j = pil.index(min(pil)), pil.index(max(pil))
        return [[l[i], l[j]] for l in A]
    else:
        l1 = int(n/2)
        ar1 = sort(*[Ai[:l1] for Ai in A])
        ar2 = sort(*[Ai[l1:] for Ai in A])
        n1, n2 = len(ar1[0]), len(ar2[0])
        final, i, j = [[] for _ in range(N)], 0, 0
        while i < n1:
            if j < n2:
                while ar1[0][i] > ar2[0][j]:
                    for f in range(N):
                        final[f].append(ar2[f][j])
                    j += 1
                    if j == n2:   break
            for k in range(N):
                final[k].append(ar1[k][i])
            i += 1
        for f in range(N):
            final[f] += ar2[f][j:]
        return final

'''
HOW TO USE

Psearch:
    Will search for best period and/or best period derivative. The period search, for given P_dot, searches for the best period that corresponds to the begenning of the events.
    Analytically:
        1. Period search: Psearch(events, periods to search, pdot, number of bins for counts-per-phase histogram)
        2. Pdot search: Psearch(events, period at the start of the events, list of pdot's to search, number of bins...)
        3. Full search (2D search): Psearch(events, period list, pdot list, bins...)
'''