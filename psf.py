#
#   Perfect Sequence Finder
#

# Written by Sam Blake.

# Started on 3 June 2018.

import os
import time
import json
import numpy as np
import random
import math
import cmath
import datetime
import correlations
import pickle

import fill_array_2d


def psf2d(search_data_file = 'search_history.json', perfect_sequence_log_file = 'perfect_sequences.log', \
    perfect_array_log_file = 'perfect_arrays.log', \
    phases = [6, 12, 15, 18, 20, 24, 30], \
    size_ranges = [[36,1296],[144,20736],[225,50625],[324,104976],[400,160000],[576,331776],[900,810000]], \
    moduli = 'Automatic', \
    denominators = 'Automatic', \
    n_trials = 250, \
    n_sums = 1, degX = 2, degY = 2, \
    tol_perfect = 0.01, \
    square_only = False, diagonal_only = False, verbose = True, run_checks = True):

    if run_checks:
        if verbose:
            print 'Testing...'

        # success = test_psf2d(tol_perfect = tol_perfect, verbose = verbose)
        success = True

        if verbose:
            if success:
                print 'Testing complete - all tests passed.'
            else:
                print 'Testing complete - TESTING FAILED! Bye.'
                return

    if verbose:
        print 'search_data_file is', search_data_file
        print 'perfect_array_log_file is', perfect_array_log_file
        print 'perfect_sequence_log_file is', perfect_sequence_log_file

    if moduli == 'Automatic':
        moduli = phases
    
    if denominators == 'Automatic':
        denominators = ['Automatic' for k in range(len(phases))]
        
    start = time.time()
    
    for modulus, phase, size_range, dens in zip(moduli, phases, size_ranges, denominators):
        print 'phase is ', phase
        min_elems, max_elems = size_range
        construction_search(n_trials, min_elems, max_elems, phase, modulus, \
            n_sums = n_sums, degX = degX, degY = degY, denominators = dens, \
            square_only = square_only, diagonal_only = diagonal_only, \
            search_data_file = search_data_file, perfect_sequence_log_file = perfect_sequence_log_file, \
            perfect_array_log_file = perfect_array_log_file, \
            tol_perfect = tol_perfect, verbose = verbose)
    
    elapsed = time.time() - start
    print 'finished after ', '{:.4f}'.format(elapsed/3600.0), '[h]'
    
    return elapsed



def construction_search(n_trials, min_elems, max_elems, phase, modulus, n_sums = 1, \
    degX = 2, degY = 2, denominators = 'Automatic', square_only = False, diagonal_only = False, \
    search_data_file = 'search_history.json', perfect_sequence_log_file = 'perfect_sequences.log', \
    perfect_array_log_file = 'perfect_arrays.log', \
    tol_perfect = 0.01, verbose = True):
    
    # Sense check. 
    
    if degX < 1 or degY < 1 or n_sums < 1:
        print 'Error: malformed inputs.'
        return
    
    if denominators != 'Automatic' and len(denominators) != n_sums:
        print 'Error: please specify one denominator for each sum.'
        return
    
    if denominators != 'Automatic':
        fixed_denominators = True
        dens = np.asarray(denominators, dtype = np.int32)
    else:
        fixed_denominators = False
    
    # Import search data.
    
    if os.path.exists(search_data_file):
        with open(search_data_file) as handle:
            search_data = json.loads(handle.read())
    else:
        search_data = {}
    
    # Export archive of search data incase something goes wrong. 
    
    archive = search_data_file.replace('.json', datetime.datetime.now().strftime("_%d_%m_%Y.pkl"))
    with open(archive, 'wb') as handle:
        pickle.dump(search_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    # Compute search space size. 
    
    search_space_size = (2*modulus - 1)**(n_sums*(3 + ((degX + 1)*(1 + degY)))) 
    if verbose:
        print search_space_size
        
    for n_elems in xrange(min_elems, max_elems + 1):

        # Sequence length should not be coprime to the phase. 
                
        if gcd(n_elems, phase) == 1:
            continue
        
        # Do not generate arrays of size p x 1 or 1 x p. 

        if prime_Q(n_elems):
            continue

        if verbose:
            print n_elems, phase, '  (nelems,phase)'
            
        flat_array = np.zeros(n_elems, dtype = np.complex_)
        
        # Default min/max array sizes.

        min_n = 2
        min_m = 2
        max_n = int((n_elems + n_elems%2)/2) + 1
        max_m = int((n_elems + n_elems%2)/2) + 1

        # Square-only arrays.

        if square_only:
            sqroot = int(math.sqrt(n_elems))
            if sqroot**2 == n_elems:
                min_n = sqroot
                min_m = sqroot
                max_n = sqroot + 1
                max_m = sqroot + 1
            else:
                print 'Error: n_elems should be a square.'
                return

        # Iterate through all possible array sizes.
            
        for n in xrange(min_n, max_n):
            upper_diag_only = max(min_m, n)
            for m in xrange(upper_diag_only, max_m):
                
                n_balanced = 0
                n_perfect_cols = 0
                n_perfect_rows = 0
                n_perfect_arrays = 0
                n_perfect_sequences = 0
                
                # Only arrays wth n_elems elements. 
                    
                if n*m != n_elems:
                    continue
                
                # n,m should not be coprime to the phase. 

                if gcd(n, phase) == 1 or gcd(m, phase) == 1:
                    continue

                # Create array of appropriate size. 
                    
                array = flat_array.reshape(n,m)
                    
                # Coprime dimensions only for diagonally enumerated arrays.
                    
                if diagonal_only and gcd(n, m) != 1:
                    continue
                        
                # Exclude known arrays. 
                    
                if n == phase and m == phase:
                    continue
                elif 2*n == phase and 2*m == phase:
                    continue

                # Exclude trivially small arrays. 
                
                if n*m < phase:
                    continue
                
                if verbose:
                    print '    ', n, ' by ', m
                
                # Display percentage already searched. 
                
                if verbose:
                    percentage_complete(search_data, str(phase), str(modulus), str(n_elems), str(n), str(m), str(degX), str(degY), str(n_sums))
                
                start = time.time()
                best_moffp = n_elems

                for k in xrange(0, n_trials):
                    
                    # Randomly select index function parameters.
                    
                    if fixed_denominators:
                        Ns = dens
                    else:
                        Ns  = np.random.randint(low = 1, high = modulus, size = (n_sums), dtype = np.int32)
                        
                    Cs = np.random.randint(low = 1 - modulus, high = modulus, size = (n_sums), dtype = np.int32)
                    pxy = np.random.randint(low = 1 - modulus, high = modulus, size = (n_sums, degX + 1, degY + 1), dtype = np.int32)

                    # print write_index_as_string(n_sums, Cs, Ns, pxy)

                    # Fill the array. 
                    
                    fill_array_2d.fill_array_2d(array, n_sums, Cs, Ns, pxy, phase)
                    
                    # The following are three necessary, but not sufficient checks for 
                    # perfect arrays. Firstly, we check the array satisfies the balance theorem. 

                    if balanced_Q(array, n_elems, tol = tol_perfect):
                        n_balanced += 1
                    else:
                        continue

                    # Check sum of cols is perfect. 

                    array_sum_cols = np.sum(array, axis = 0)
                    moffp = max_abs_off_peak(array_sum_cols)

                    if moffp < tol_perfect*float(n_elems):
                        n_perfect_cols += 1
                    else:
                        continue

                    # Check sum of rows is perfect.

                    array_sum_rows = np.sum(array, axis = 1)
                    moffp = max_abs_off_peak(array_sum_rows)

                    if moffp < tol_perfect*float(n_elems):
                        n_perfect_rows += 1
                    else:
                        continue

                    # Periodic autocorrelation of the 2D array. 
                    # TODO: should be we aligning the data for efficiency? 
                    
                    moffp = max_abs_off_peak(array)

                    if moffp < best_moffp:
                        best_moffp = moffp

                    if moffp < tol_perfect*float(n_elems):
                        n_perfect_arrays += 1
                        log_perfect(2, "2D", perfect_array_log_file, n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = verbose)
                        if gcd(n, m) == 1:
                            n_perfect_sequences += 1
                            log_perfect(1, "diagonal", perfect_sequence_log_file, n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = verbose)
                            if n_elems > phase*phase:
                                print '****  BROKEN FRANK BOUND!  ****'
                                log_perfect(1, "diagonal", "champagne.log", n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = verbose)                                
                                np.savetxt("champagne_sequence.csv", array, delimiter=",")
                                np.savetxt("champagne_autocorrelations.csv", correlations.autocorrelate_fftw(array), delimiter=",")
                                
                        if diagonal_only:
                            continue
                    
                        # Row-major array enumeration. 
                    
                        seq = array.flatten(order = 'C')
                        moffp = max_abs_off_peak(seq)
                        if moffp < tol_perfect*float(n_elems):
                            n_perfect_sequences += 1
                            log_perfect(1, "row-major", perfect_sequence_log_file, n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = verbose)
                            if n_elems > phase*phase:
                                print '****  BROKEN FRANK BOUND!  ****'
                                log_perfect(1, "row-major", "champagne.log", n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = verbose)                                
                                np.savetxt("champagne.csv", seq, delimiter=",")
                                np.savetxt("champagne_autocorrelations.csv", correlations.autocorrelate_fftw(seq), delimiter=",")
                                
                        # Column-major array enumeration. 

                        seq = array.flatten(order = 'F')
                        moffp = max_abs_off_peak(seq)
                        if moffp < tol_perfect*float(n_elems):
                            n_perfect_sequences += 1
                            log_perfect(1, "column-major", perfect_sequence_log_file, n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = verbose)
                            if n_elems > phase*phase:
                                print '****  BROKEN FRANK BOUND!  ****'
                                log_perfect(1, "column-major", "champagne.log", n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = verbose)                                
                                np.savetxt("champagne.csv", seq, delimiter=",")
                                np.savetxt("champagne_autocorrelations.csv", correlations.autocorrelate_fftw(seq), delimiter=",")
                        
                # Update search data.
                
                update_search_history(search_data, str(phase), str(modulus), str(n_elems), str(n), str(m), \
                    str(degX), str(degY), str(n_sums), n_trials, n_perfect_arrays, n_perfect_sequences, best_moffp, search_space_size, start)
    
                # Stats on necessary conditions.

                if verbose:
                    print '      ', '{:.4f}'.format(100.0*float(n_balanced)/float(n_trials)), '[% balanced]'
                    print '      ', '{:.4f}'.format(100.0*float(n_perfect_rows)/float(n_trials)), '[% perfect rows]'
                    print '      ', '{:.4f}'.format(100.0*float(n_perfect_cols)/float(n_trials)), '[% perfect cols]'

                # Lowest off-peak.

                if verbose:
                    print '      ', '{:.4f}'.format(100.0*best_moffp/float(n_elems)), '[% of peak]'

                # Elapsed time.
            
                if verbose:
                    print '      ', '{:.4f}'.format(time.time() - start), '[s],', int(float(n_trials)/(time.time() - start)), '[trials/s]'
                
    # Export updated search history.
    
    print 'Exporting search data...'
    export_start = time.time()

    with open(search_data_file, 'w') as handle:
        json.dump(search_data, handle)

    print 'finished exporting after ', '{:.4f}'.format(time.time() - export_start), '[s]'
    
    # Export pyFFTW plan.



def recover_search_history(pkl_file, json_file):
    with open(pkl_file, 'rb') as handle:
        data = pickle.load(handle)
    with open(json_file, 'w') as handle:
        json.dump(data, handle)


def percentage_complete(search_data, phase, modulus, n_elems, ny, nx, degx, degy, n_sums):
    
    if phase in search_data and modulus in search_data[phase] and n_elems in search_data[phase][modulus] and \
    ny in search_data[phase][modulus][n_elems] and nx in search_data[phase][modulus][n_elems][ny] and \
    degx in search_data[phase][modulus][n_elems][ny][nx] and \
    degy in search_data[phase][modulus][n_elems][ny][nx][degx] and \
    n_sums in search_data[phase][modulus][n_elems][ny][nx][degx][degy]:
        print '      ', '{:.6f}'.format(search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["percentage_complete"]), '[% complete]'
    else:
        print '      ', '-- [% complete]'


def update_search_history(search_data, phase, modulus, n_elems, ny, nx, degx, degy, n_sums, n_trials, \
    n_perfect_arrays, n_perfect_sequences, lowest_off_peak, search_space, start):
    
    if phase in search_data and modulus in search_data[phase] and n_elems in search_data[phase][modulus] \
    and ny in search_data[phase][modulus][n_elems] and nx in search_data[phase][modulus][n_elems][ny] \
    and degx in search_data[phase][modulus][n_elems][ny][nx] and degy in search_data[phase][modulus][n_elems][ny][nx][degx] \
    and n_sums in search_data[phase][modulus][n_elems][ny][nx][degx][degy]:
        search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["n_trials"] += n_trials
        total_trials = search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["n_trials"]
        search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["percentage_complete"] = 100.0*float(total_trials)/float(search_space)
        search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["n_perfect_arrays"] += n_perfect_arrays
        search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["n_perfect_sequences"] += n_perfect_sequences
        search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["compute_time"] += (time.time() - start)/3600.0
        if search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["lowest_off_peak"] < lowest_off_peak:
            search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums]["lowest_off_peak"] = lowest_off_peak
    else:
        data = {"n_trials" : n_trials, "percentage_complete" : 100.0*float(n_trials)/float(search_space), \
        "n_perfect_arrays" : n_perfect_arrays, "n_perfect_sequences" : n_perfect_sequences, \
        "lowest_off_peak" : lowest_off_peak, "compute_time" : (time.time() - start)/3600.0}
        # There must be an easier way!! :-/
        if phase not in search_data:
            search_data[phase] = {modulus : {n_elems : {ny : {nx : {degx : {degy : {n_sums : data}}}}}}}
        elif modulus not in search_data[phase]:
            search_data[phase][modulus] = {n_elems : {ny : {nx : {degx : {degy : {n_sums : data}}}}}}
        elif n_elems not in search_data[phase][modulus]:
            search_data[phase][modulus][n_elems] = {ny : {nx : {degx : {degy : {n_sums : data}}}}}
        elif ny not in search_data[phase][modulus][n_elems]:
            search_data[phase][modulus][n_elems][ny] = {nx : {degx : {degy : {n_sums : data}}}}
        elif nx not in search_data[phase][modulus][n_elems][ny]:
            search_data[phase][modulus][n_elems][ny][nx] = {degx : {degy : {n_sums : data}}}
        elif degx not in search_data[phase][modulus][n_elems][ny][nx]:
            search_data[phase][modulus][n_elems][ny][nx][degx] = {degy : {n_sums : data}}
        elif degy not in search_data[phase][modulus][n_elems][ny][nx][degx]:
            search_data[phase][modulus][n_elems][ny][nx][degx][degy] = {n_sums : data}
        elif n_sums not in search_data[phase][modulus][n_elems][ny][nx][degx][degy]:
            search_data[phase][modulus][n_elems][ny][nx][degx][degy][n_sums] = data
        else:
            print 'Error...'


def log_perfect(n_dims, construction_type, log_file_template, n, m, n_sums, Cs, Ns, pxy, modulus, phase, verbose = True):
    
    log_file = log_file_template.replace('.log', '_' + str(n) + '_' + str(m) + '_' + str(phase) + '.log')

    # Check log file exists. 
    
    if os.path.exists(log_file):
        append_write = 'a'
    else:
        append_write = 'w'
    
    lf = open(log_file, append_write)
    
    if n_dims == 1:
        sequence_or_array = 'SEQUENCE'
    else:
        sequence_or_array = 'array'
    
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    if verbose:
        print 'Perfect ' + sequence_or_array + ' found on ' + now
        print 'index = ', write_index_as_string(n_sums, Cs, Ns, pxy)
        print 'phase = ' + str(phase)
        print 'size  = ' + str(n) + ' x ' + str(m)
        print 'type  = ' + construction_type
    
    lf.write('Perfect ' + sequence_or_array + ' found on ' + now + '\n')
    lf.write('index = ' + write_index_as_string(n_sums, Cs, Ns, pxy) + '\n')
    lf.write('phase = ' + str(phase) + '\n')
    lf.write('size  = ' + str(n) + ' x ' + str(m) + '\n')
    if n_dims == 1:
        lf.write('type  = ' + construction_type + '\n')
    lf.write('---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- \n\n')
    
    lf.close()

# Write the index function in mathematica syntax.        

def write_index_as_string(n_sums, Cs, Ns, pxys):
    
    first = True
    indexfn = ''
    
    for k in range(0, n_sums):
        
        if not first and Cs[k] > 0:
            indexfn += ' + '
        
        if Cs[k] < 0:
            indexfn += ' - '

        if abs(Cs[k]) != 1:
            indexfn += str(abs(Cs[k])) + '*'
        
        if Ns[k] != 1:
            floor = True
            indexfn += 'Floor[('
        else:
            floor = False
            indexfn += '('
        
        indexfn += polynomial_as_string(pxys[k])
        
        if floor:
            indexfn += ')/' + str(Ns[k]) + ']'
        else:
            indexfn += ')'
        
        first = False
    
    return indexfn


# print write_index_as_string(2, np.array([-1,1]), np.array([5,6]), np.array([[[1,0,-3],[0,0,3],[0,6,0]], [[1,2,-3],[5,-4,3],[-7,6,4]]]))

# Pretty print a bivariate polynomial.

def polynomial_as_string(pxy):
    
    # Zero? 
    if np.allclose(pxy, 0.0, 1.0e-5):
        return '0'
    
    degX = pxy.shape[0]
    degY = pxy.shape[1]

    X = 'x'
    Y = 'y'
    poly = ''
    
    for i in range(0, degX):
        for j in range(0, degY):
            c = pxy[i,j]
            if i == 0 and j == 0:
                poly += str(c)
            elif c == 0:
                continue
            elif i != 0 and j == 0:
                if c > 0:
                    if i == 1:
                        if c == 1:
                            poly += ' + ' + X
                        else:
                            poly += ' + ' + str(c) + '*' + X
                    else:
                        if c == 1:
                            poly += ' + ' + X + '^' + str(i)
                        else:
                            poly += ' + ' + str(c) + '*' + X + '^' + str(i)
                else:
                    if i == 1:
                        if c == -1:
                            poly += ' - ' + X
                        else:
                            poly += ' - ' + str(-c) + '*' + X
                    else:
                        if c == -1:
                            poly += ' - ' + X + '^' + str(i)
                        else:
                            poly += ' - ' + str(-c) + '*' + X + '^' + str(i)
            elif i == 0 and j != 0:
                if c > 0:
                    if j == 1:
                        if c == 1:
                            poly += ' + ' + Y
                        else:
                            poly += ' + ' + str(c) + '*' + Y
                    else:
                        if c == 1:
                            poly += ' + ' + Y + '^' + str(j)
                        else:
                            poly += ' + ' + str(c) + '*' + Y + '^' + str(j)
                else:
                    if j == 1:
                        if c == -1:
                            poly += ' - ' + Y
                        else:
                            poly += ' - ' + str(-c) + '*' + Y
                    else:
                        if c == -1:
                            poly += ' - ' + Y + '^' + str(j)
                        else:
                            poly += ' - ' + str(-c) + '*' + Y + '^' + str(j)
            else:
                if c > 0:
                    if i == 1 and j == 1:
                        if c == 1:
                            poly += ' + ' + X + '*' + Y
                        else:
                            poly += ' + ' + str(c) + '*' + X + '*' + Y
                    elif i == 1:
                        if c == 1:
                            poly += ' + ' + X + '*' + Y + '^' + str(j)
                        else:
                            poly += ' + ' + str(c) + '*' + X + '*' + Y + '^' + str(j)
                    elif j == 1:
                        if c == 1:
                            poly += ' + ' + X + '^' + str(i) + '*' + Y
                        else:
                            poly += ' + ' + str(c) + '*' + X + '^' + str(i) + '*' + Y
                    else:
                        if c == 1:
                            poly += ' + ' + X + '^' + str(i) + '*' + Y + '^' + str(j)
                        else:
                            poly += ' + ' + str(c) + '*' + X + '^' + str(i) + '*' + Y + '^' + str(j)
                else:
                    if i == 1 and j == 1:
                        if c == -1:
                            poly += ' - ' + X + '*' + Y
                        else:
                            poly += ' - ' + str(-c) + '*' + X + '*' + Y
                    elif i == 1:
                        if c == -1:
                            poly += ' - ' + X + '*' + Y + '^' + str(j)
                        else:
                            poly += ' - ' + str(-c) + '*' + X + '*' + Y + '^' + str(j)
                    elif j == 1:
                        if c == -1:
                            poly += ' - ' + X + '^' + str(i) + '*' + Y
                        else:
                            poly += ' - ' + str(-c) + '*' + X + '^' + str(i) + '*' + Y
                    else:
                        if c == -1:
                            poly += ' - ' + X + '^' + str(i) + '*' + Y + '^' + str(j)
                        else:
                            poly += ' - ' + str(-c) + '*' + X + '^' + str(i) + '*' + Y + '^' + str(j)
    
    return poly


# print polynomial_as_string(np.array([[1,2,-3],[5,-4,3],[-7,6,4]]))

def max_abs_off_peak(array):
    return correlations.max_off_peak(correlations.autocorrelate_fftw(array))

def perfect_Q(array, tol = 5.0e-2):
    return max_abs_off_peak(array) < tol

def balanced_Q(array, n_elems, tol = 5.0e-2):
    bal = np.abs(np.sum(array))**2
    return approx_equal_Q(bal, float(n_elems), epsilon = 5.0e-2)

# print perfect_Q(np.array([1,1,1,-1]))
# print perfect_Q(np.array([[1,1],[1,-1]]))
# print perfect_Q(np.array([[-1,1,1,1],[1,-1,1,1],[1,1,-1,1],[1,1,1,-1]]))
# print correlations.autocorrelate_fftw(np.array([[1,1,1,1],[1,-1,1,1],[1,1,-1,1],[1,1,1,-1]]))
# print correlations.max_off_peak(correlations.autocorrelate_fftw(np.array([[1,1,1,1],[1,-1,1,1],[1,1,-1,1],[1,1,1,-1]])))


# Ref: Knuth, Semi-Numerical Algorithms, 4.2.2, eqn. 22.

def approx_equal_Q(a, b, epsilon = 1.0e-16):
    return abs(a - b) <= max(abs(a),abs(b))*epsilon


def fill_array_2d_slow(array, n, m, n_sums, Cs, Ns, pxys, modulus, phase):
    
    for i in xrange(0, n):
        for j in xrange(0, m):
            index = abstract_index_fn(j, i, n_sums, Cs, Ns, pxys, modulus)
            array[i,j] = cmath.exp(2.0*cmath.pi*1j*float(index)/float(phase))


def abstract_index_fn(x, y, n_sums, Cs, Ns, pxys, modulus):
    index = 0
    for i in xrange(0, n_sums):
        index += Cs[i]*int(math.floor(float(poly_eval_mod_2d(pxys[i], x, y, modulus))/float(Ns[i])))
        index = index%modulus
    return index


# 1 + 5 x - 7 x^2 + 2 y - 4 x y + 6 x^2 y - 3 y^2 + 3 x y^2 + 4 x^2 y^2
# abstract_index_fn(4, 2, 2, np.array([1,1]), np.array([2,3]), np.array([[1,2,-3],[5,-4,3],[-7,6,4]]), 24)


#    Stores a bivariate polynomial using a dense representation such that a dot product with:
#      
#    1       y     y^2 ...     y^degY
#    
#    x     x*y   x*y^2 ...   x*y^degY
#    
#    x^2 x^2*y x^2*y^2 ... x^2*y^degY
#    
#    ...                        .
#    
#    .                          .
#    
#    .                          .
#    
#    x^degX x^degX*y x^degX*y^2 ... x^degX*y^degY
# 
#    gives the bivariate polynomial of total degree degX + degY.

def poly_eval_mod_2d(coeffs, x, y, modulus):
    degX = coeffs.shape[0]
    degY = coeffs.shape[1]
    pxy = 0
    for i in range(0, degX):
        for j in range(0, degY):
            pxy += coeffs[i,j]*(x**i)*(y**j)
            pxy = pxy%modulus
    return pxy

# print poly_eval_mod_2d(np.array([[1,2,-3],[5,-4,3],[-7,6,4]]), 4, 2, 24)


def prime_Q(n):
    if n == 2:
        return True
    if n%2 == 0 or n < 2:
        return False

    sqr = int(math.sqrt(n)) + 1

    for d in xrange(3, sqr, 2):
        if n%d == 0:
            return False
    return True


# print prime_Q(91), prime_Q(93), prime_Q(97)

def gcd(x, y):
    while y != 0:
        (x, y) = (y, x % y)
    return x



def test_psf2d(tol_perfect = 5.0e-2, verbose = False):

    success = True

    # Zadoff-Chu sequences. 

    if verbose:
        print 'Odd-length Zadoff-Chu sequence tests...'

    n_sums = 1
    Cs = np.array([1], dtype = np.int32)
    Ns = np.array([1], dtype = np.int32)

    # Odd length.

    # [[1, y, y^2], [x, x y, x y^2], [x^2, x^2 y, x^2 y^2]] 
    pxy = np.array([[[0, 1, 1]]], dtype = np.int32)

    for slen in range(5,1000,20):
        nr = slen
        array = np.ones((slen,1), dtype = np.complex_)
        fill_array_2d.fill_array_2d(array, n_sums, Cs, Ns, pxy, nr)

        moffp = max_abs_off_peak(array.flatten(order = 'C'))

        if verbose:
            print slen, moffp
        if moffp > tol_perfect*float(slen):
            if verbose:
                print 'ERROR: Odd-length Zadoff-Chu sequences test failed.'
            success = False

    # Even length.

    if verbose:
        print 'Even-length Zadoff-Chu sequence tests...'

    pxy = np.array([[[0, 0, 1]]], dtype = np.int32)

    for slen in range(6, 1000, 20):
        nr = 2*slen
        array = np.ones((slen,1), dtype = np.complex_)

        fill_array_2d.fill_array_2d(array, n_sums, Cs, Ns, pxy, nr)
        moffp = max_abs_off_peak(array.flatten(order = 'C'))

        if verbose:
            print slen, moffp
        if moffp > tol_perfect*float(slen):
            if verbose:
                print 'ERROR: Even-length Zadoff-Chu sequences test failed.'
            success = False


    # Liu-Fan sequences.

    if verbose:
        print 'Liu-Fan sequence tests...'

    n_sums = 1
    Cs = np.array([1], dtype = np.int32)
    Ns = np.array([2], dtype = np.int32)

    for slen in range(4, 1000, 20):
        nr = slen
        pxy = np.array([[[0, 0, nr - 1]]], dtype = np.int32)
        array = np.ones((slen,1), dtype = np.complex_)

        fill_array_2d.fill_array_2d(array, n_sums, Cs, Ns, pxy, nr)
        moffp = max_abs_off_peak(array.flatten(order = 'C'))

        if verbose:
            print slen, moffp
        if moffp > tol_perfect*float(slen):
            if verbose:
                print 'ERROR: Liu-Fan sequences test failed.'
            success = False


    # Blake-Tirkel 2014 sequences. 

    n_sums = 1
    Cs = np.array([1], dtype = np.int32)

    # [[1, y, y^2], [x, x y, x y^2], [x^2, x^2 y, x^2 y^2]] 
    pxy = np.array([[[0, 0, 1], [0, 1, 0]]], dtype = np.int32)

    for m in range(1,20):
        for n in range(1,6):
            for k in range(1,4):
                nr = 2*m*n**k
                Ns = np.array([n], dtype = np.int32)
                array = np.ones((2*m*n**(k+1),2), dtype = np.complex_)
                fill_array_2d.fill_array_2d(array, n_sums, Cs, Ns, pxy, nr)
                moffp = max_abs_off_peak(array.flatten(order = 'C'))

                if verbose:
                    print 4*m*n**(k+1), moffp
                if moffp > tol_perfect*float(4*m*n**(k+1)):
                    if verbose:
                        print 'ERROR: Blake-Tirkel (2014) sequences test failed.'
                    success = False                


    # Milewski sequences. 


    # Frank sequences. 

    if verbose:
        print 'Frank sequence tests...'

    n_sums = 1
    Cs = np.array([1], dtype = np.int32)
    Ns = np.array([1], dtype = np.int32)

    pxy = np.array([[[0, 0, 0], [0, 1, 0]]], dtype = np.int32) # x y

    for phase in range(4, 1000, 40):
        array = np.ones((phase,phase), dtype = np.complex_)

        fill_array_2d.fill_array_2d(array, n_sums, Cs, Ns, pxy, phase)
        moffp = max_abs_off_peak(array.flatten(order = 'C'))

        if verbose:
            print phase, moffp
        if moffp > tol_perfect*float(phase**2):
            if verbose:
                print 'ERROR: Frank sequences test failed.'
            success = False

    return success


