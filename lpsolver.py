import numpy
import sys
debug = True


# Takes in default standard form LP as a list and outputs the same LP as its corresponding A, b, and c vectors
def standard_form_to_eq(lp):
    A = []
    b = []
    for coeffs in lp[1:]:
        A.append(coeffs[:-1])
        b.append(coeffs[-1])
    A = numpy.array(A)
    c = lp[0] + [0]*len(b)
    #Add the identity basis (slack variables)
    A = numpy.concatenate((A,numpy.identity(A.shape[0])), axis=1)
    return (A,b,c)

def determine_pivot_index(dictionary):
    # Largest Coefficient
    return numpy.argmax(dictionary[0])


def primal_simplex(A,b,c,B,N):
    return

def dual_simplex():
    return


def main():
    # Read in LP from stdin
    lp = []
    for line in sys.stdin:
        coeff_line = []
        for coeff in line.split():
            coeff_line.append(float(coeff))
        lp.append(coeff_line)

    # Convert LP into dictionary matrix form
    (A,b,c) = standard_form_to_eq(lp)

    if(debug):
        print(A)
        print(b)
        print(c)




if __name__ == "__main__":
    main()

