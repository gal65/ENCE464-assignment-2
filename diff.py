#!/usr/bin/python3

import sys

with open(sys.argv[1]) as ref:
    arr = [[float(x) for x in line.split()] for line in ref]

arr1 = [[float(x) for x in line.split()] for line in sys.stdin]

for refl, testl in zip(arr, arr1):
    for refd, testd in zip(refl, testl):
        print("{: .5f} ".format(refd - testd), end="")
    print()

