
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
import psf


small_scale = 20
medium_scale = 200
large_scale = 1000

scale = medium_scale


psf.test_psf2d(tol_perfect = 5.0e-2, verbose = True)


for k in range(1000):

  # Some specific 24-phase searches.

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[2]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[3]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[4]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[5]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[3*4]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[3*8]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[6*8]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[9*8]], n_sums = 1, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)
  
  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[2,3]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[4,3]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[2,5]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[4,5]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[3,5]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[3,6]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[2,8]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[3,8]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[6,8]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[9,8]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[9,16]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[6*8,6*8]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)

  psf.psf2d(phases = [24], size_ranges = [[1728, 1728]], denominators = [[6*16,6*16]], n_sums = 2, degX = 2, degY = 2, n_trials = 150*scale, diagonal_only = True)
  
# # 4-phase

  phase = 4
  mod = 4
  trials = 1 # scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)


# # 6-phase

  psf.psf2d(phases = [6], moduli = [3], size_ranges = [[8*9, 8*9]], denominators = [[2]], n_sums = 1,
          degX = 2, degY = 2, n_trials = 500000, diagonal_only = True)

  psf.psf2d(phases = [6], moduli = [3], size_ranges = [[8*9, 8*9]], denominators = [[3]], n_sums = 1, \
          degX = 2, degY = 2, n_trials = 500000, diagonal_only = True)

  psf.psf2d(phases = [6], moduli = [3], size_ranges = [[8*9, 8*9]], denominators = [[2,3]], n_sums = 2, \
          degX = 2, degY = 2, n_trials = 1000000, diagonal_only = True)

  psf.psf2d(phases = [6], moduli = [3], size_ranges = [[8*27, 8*27]], denominators = [[2]], n_sums = 1, \
          degX = 2, degY = 2, n_trials = 1000000, diagonal_only = True)

  psf.psf2d(phases = [6], moduli = [3], size_ranges = [[8*27, 8*27]], denominators = [[3]], n_sums = 1, \
          degX = 2, degY = 2, n_trials = 1000000, diagonal_only = True)


  psf.psf2d(phases = [6], moduli = [3], size_ranges = [[8*27, 8*27]], denominators = [[2,3]], n_sums = 2, \
          degX = 2, degY = 2, n_trials = 1000000, diagonal_only = True)

  psf.psf2d(phases = [6], size_ranges = [[2*36, 2*36]], n_sums = 1, degX = 2, degY = 2, n_trials = 100000)

  psf.psf2d(phases = [6], size_ranges = [[2*36, 2*36]], n_sums = 2, degX = 2, degY = 2, n_trials = 100000)

  psf.psf2d(phases = [6], size_ranges = [[3*36, 3*36]], n_sums = 1, degX = 2, degY = 2, n_trials = 100000)

  psf.psf2d(phases = [6], size_ranges = [[3*36, 3*36]], n_sums = 2, degX = 2, degY = 2, n_trials = 100000)

  psf.psf2d(phases = [6], size_ranges = [[4*36, 4*36]], n_sums = 1, degX = 2, degY = 2, n_trials = 10000)

  psf.psf2d(phases = [6], size_ranges = [[4*36, 4*36]], n_sums = 2, degX = 2, degY = 2, n_trials = 10000)


  phase = 6
  mod = 6
  trials = 15*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)


# # 8-phase

  phase = 8
  mod = 8
  trials = 10*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

# # 10-phase

  phase = 10
  mod = 10
  trials = 75*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

# # 12-phase

  phase = 12
  mod = 15
  trials = 50*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

# # 15-phase

  phase = 15
  mod = 15
  trials = 30*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

# # 16-phase

  phase = 16
  mod = 16
  trials = 25*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, 8*phase**2]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**3, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

# # 18-phase

  phase = 18
  mod = 18
  trials = 20*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, 9*phase**2]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**3, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

# # 20-phase

  phase = 20
  mod = 20
  trials = 15*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, 10*phase**2]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**3, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

# # 24-phase

  phase = 24
  mod = 24
  trials = 12*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, 12*phase**2]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**3, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

# # 30-phase

  phase = 30
  mod = 30
  trials = 5*scale
  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**2, 15*phase**2]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**3, phase**3]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)

  for ns in range(1,4):
    for deg in range(1,4):
      psf.psf2d(phases = [phase], moduli = [mod], size_ranges = [[phase**4, phase**4]], \
                  n_sums = ns, degX = deg, degY = deg, n_trials = trials, square_only = True)


  psf.psf2d(phases = [30], moduli = [2], size_ranges = [[30**3, 30**3]], n_sums = 3, degX = 2, degY = 2, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)

  psf.psf2d(phases = [30], moduli = [30], size_ranges = [[30**3, 30**3]], n_sums = 3, degX = 2, degY = 2, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)

  psf.psf2d(phases = [30], moduli = [2], size_ranges = [[30**3, 30**3]], n_sums = 3, degX = 3, degY = 3, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)

  psf.psf2d(phases = [30], moduli = [30], size_ranges = [[30**3, 30**3]], n_sums = 3, degX = 3, degY = 3, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)


  psf.psf2d(phases = [30], moduli = [2], size_ranges = [[30**4, 30**4]], n_sums = 3, degX = 2, degY = 2, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)

  psf.psf2d(phases = [30], moduli = [30], size_ranges = [[30**4, 30**4]], n_sums = 3, degX = 3, degY = 3, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)

  psf.psf2d(phases = [30], moduli = [2], size_ranges = [[30**4, 30**4]], n_sums = 3, degX = 2, degY = 2, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)

  psf.psf2d(phases = [30], moduli = [30], size_ranges = [[30**4, 30**4]], n_sums = 3, degX = 3, degY = 3, \
          denominators = [[2, 3, 5]], n_trials = 25*scale, square_only = True)



