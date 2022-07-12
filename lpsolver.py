from operator import matmul
import numpy
import sys
debug = False


# Takes in default standard form LP as a list and outputs the same LP as its corresponding A, b, and c vectors
def standard_form_to_eq(lp):
    A = []
    b = []
    for coeffs in lp[1:]:
        A.append(coeffs[:-1])
        b.append(coeffs[-1])
    A = numpy.array(A)
    N = [i for i in range(len(lp[0]))]
    B = [i for i in range(len(N),len(N)+len(b))]
    c = lp[0] + [0]*len(b)
    #Add the identity basis (slack variables)
    A = numpy.concatenate((A,numpy.identity(A.shape[0])), axis=1)
    return (A,b,c,B,N)

def determine_pivot_index(dictionary):
    # Largest Coefficient
    return numpy.argmax(dictionary[0])


def primal_simplex(A,b,c,B,N):
    # First let's make our basis matrix A_b (columns of A corresponding to our basis variables)
    A_B = A[:, B[0]]
    for basis_index in B[1:]:
        A_column = A[:, basis_index]
        A_B = numpy.column_stack((A_B, A_column))

    # Non basis matrix
    A_N = A[:, N[0]]
    for non_basis_index in N[1:]:
        A_column = A[:, basis_index]
        A_N = numpy.column_stack((A_N, A_column))
    
    # Take its inverse (TODO CHANGE LATER TO NOT USE INVERSE)
    A_B_i = numpy.linalg.inv(A_B)

    # Compute initial value of x
    x_b = numpy.matmul(A_B_i, b)
    x_n = numpy.zeros(len(N))

    while True:
        # Compute our required variables based on B and N

        # Basis Matrix
        A_B = A[:, B[0]]
        for basis_index in B[1:]:
            A_column = A[:, basis_index]
            A_B = numpy.column_stack((A_B, A_column))

        # Non basis matrix
        A_N = A[:, N[0]]
        for non_basis_index in N[1:]:
            A_column = A[:, basis_index]
            A_N = numpy.column_stack((A_N, A_column))

        #Objective Vectors
        c_B = numpy.array(c)[B]
        c_N = numpy.array(c)[N]
        z_B = numpy.zeros(len(x_b))
        z_N = numpy.matmul(numpy.matmul(A_B_i, A_N).transpose(), c_B) - c_N

        #Check optimality
        if(numpy.all(z_N >= 0)):
            print("OPTIMAL!")
            print(z_N)
            #TODO Compute Eta*
            return

        # Choose pivot variable:
        

        
    
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
    (A,b,c,B,N) = standard_form_to_eq(lp)

    if(debug):
        print(A)
        print(b)
        print(c)
        print(B)
        print(N)
        

    #Initial Basis and N
    primal_simplex(A,b,c,B,N)


if __name__ == "__main__":
    main()

