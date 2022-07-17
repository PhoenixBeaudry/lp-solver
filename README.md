# Linear Program Solver - Phoenix Beaudry - V00869986

## Running Instructions
To run this program use:
```console
user@linux:~$ python3 lpsolver.py < input.txt 
```
with a well formatted LP in input.txt

## Solver Architecture
This program uses a Linear Algebraic Simplex Method with Primal Dual methods for solving initially infeasible LPs. 

It can be run using Blands rule for pivoting or largest coefficient, it uses largest coefficient by default.

It relies on the Python Numpy package for Linear Algebra operations such as solving linear equations and matrix inversions.

## Novel Features
Instead of using the dictionary based Simplex method it uses a linear algebraic form of the Simplex Method.

It can be run using Blands rule for pivoting or largest coefficient, it uses largest coefficient by default.

Additionally instead of solving an auxillary $\Omega$ variable proxy problem it uses a Primal-Dual method for solving initially infeasible linear programs.

