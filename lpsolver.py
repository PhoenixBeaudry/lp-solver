import numpy
import sys
debug = True


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

    dual = -1 * dictionary.transpose()

    if (debug):
        print("Primal Dictionary:")
        print(dictionary)
        print("Dual Dictionary:")
        print(dual)

    # Check primal-dual feasibility
    if (numpy.all(dictionary[:,0] >= 0)):
        print("Dictionary is Primal Feasible")
    if (numpy.all(dictionary[0] <= 0)):
        print("Dictionary is Dual Feasible")




if __name__ == "__main__":
    main()

