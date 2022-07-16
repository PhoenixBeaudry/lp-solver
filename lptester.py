import sys
import os
import lpsolver
tests = ["optimal", "infeasible", "unbounded"]

testfiles = "/tests/" + tests[0]

for filename in os.listdir(os.getcwd() + testfiles):
    with open(os.getcwd() + testfiles + "/" + filename) as f:
        lines = f.readlines()
        # Read in LP from file
        lp = []
        for line in lines:
            coeff_line = []
            for coeff in line.split():
                coeff_line.append(float(coeff))
            lp.append(coeff_line)
        print("Solving: ", filename)
        lpsolver.solve_lp(lp)

testfiles = "/tests/" + tests[1]
for filename in os.listdir(os.getcwd() + testfiles):
    with open(os.getcwd() + testfiles + "/" + filename) as f:
        lines = f.readlines()
        # Read in LP from file
        lp = []
        for line in lines:
            coeff_line = []
            for coeff in line.split():
                coeff_line.append(float(coeff))
            lp.append(coeff_line)
        print("Solving: ", filename)
        lpsolver.solve_lp(lp)

testfiles = "/tests/" + tests[2]
for filename in os.listdir(os.getcwd() + testfiles):
    with open(os.getcwd() + testfiles + "/" + filename) as f:
        lines = f.readlines()
        # Read in LP from file
        lp = []
        for line in lines:
            coeff_line = []
            for coeff in line.split():
                coeff_line.append(float(coeff))
            lp.append(coeff_line)
        print("Solving: ", filename)
        lpsolver.solve_lp(lp)
    