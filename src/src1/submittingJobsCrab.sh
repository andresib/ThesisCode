#!/bin/bash
crab -submit 1-400 -c $1
crab -submit 401-800 -c $1
crab -submit 801-1200 -c $1
crab -submit 1201-1600 -c $1
crab -submit 1601-2000 -c $1
crab -submit 2001-2400 -c $1
crab -submit 2401-2800 -c $1
crab -submit 2801-3200 -c $1
crab -submit 3201-3600 -c $1
crab -submit 3601-4000 -c $1
crab -submit 4001-4400 -c $1
crab -submit 4401-4800 -c $1
crab -submit 4801-5100 -c $1
exit 0

