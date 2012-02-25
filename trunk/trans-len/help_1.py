#!/usr/bin/env python

# Copyright (C) 2012 Tianyang Li
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License

"""
functions used in probability calculations
"""

from nzmath.combinatorial import stirling2
from scipy import comb, factorial

def composition(r, n, k):
    """
    number of integer solutions
    $x_1 + x_2 + ... + x_n = r$
    $0 < x_i \leq k$
    """
    comp = 0
    for m in range(n + 1):
        if r - 1 - m * k < n - 1:
            break
        if m % 2 == 0:
            comp += (comb(n, m, exact=True) * comb(r - 1 - m * k, n - 1, exact=True))
        else:
            comp -= (comb(n, m, exact=True) * comb(r - 1 - m * k, n - 1, exact=True))
    return comp

def read_split_contig(N, r, k):
    """
    number of ways to get 
    $Y_n - Y_1 = r$
    where there are $n$ distinct starting positions
    satisfying $0 < Y_{i + 1} - Y_i \leq k$
    """
    splits = 0
    if r == 0:
        splits = 1
    for n in range(2, min(r + 2, N + 1)):
        splits += (factorial(n, exact=True) * stirling2(N, n) * composition(r, n - 1, k))
    return splits

