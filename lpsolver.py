import numpy
import sys
debug = False


# Takes in default standard form LP as a list and outputs the same LP as its corresponding dictionary matrix
def standard_form_to_dictionary(lp):
    non_basic = [[0] + lp[0]]
    basics = []
    for coeffs in lp[1:]:
        basics.append([coeffs[-1]] + [-1*i for i in coeffs[:-1]])
    dictionary = numpy.array(non_basic + basics)
    return dictionary





def main():
    # Read in LP from stdin
    lp = []
    for line in sys.stdin:
        coeff_line = []
        for coeff in line.split():
            coeff_line.append(float(coeff))
        lp.append(coeff_line)

    # Convert LP into dictionary matrix form
    dictionary = standard_form_to_dictionary(lp)

    if (debug):
        print(dictionary)


if __name__ == "__main__":
    main()

