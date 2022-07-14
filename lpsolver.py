from operator import matmul
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
    N = [i for i in range(len(lp[0]))]
    B = [i for i in range(len(N),len(N)+len(b))]
    c = lp[0] + [0]*len(b)
    #Add the identity basis (slack variables)
    A = numpy.concatenate((A,numpy.identity(A.shape[0])), axis=1)
    b = numpy.array(b)
    c = numpy.array(c)
    return (A,b,c,B,N)

def determine_pivot_index(dictionary):
    # Largest Coefficient
    return numpy.argmax(dictionary[0])


def primal_simplex(A,b,c,B,N):
    if debug: print("Starting PrimalSimplex")
    iteration = 1
    while True:
        # Compute our required variables based on B and N
        if debug: print("Iteration: ", iteration)

        if debug: print("Current Basis: ", B)
        if debug: print("Current Non Basis: ", N)

        # First let's make our basis matrix A_b (columns of A corresponding to our basis variables)
        A_B = A[:, B[0]]
        for basis_index in B[1:]:
            A_column = A[:, basis_index]
            A_B = numpy.column_stack((A_B, A_column))

        # Non basis matrix
        A_N = A[:, N[0]]
        for non_basis_index in N[1:]:
            A_column = A[:, non_basis_index]
            A_N = numpy.column_stack((A_N, A_column))
        
        # Take its inverse (TODO CHANGE LATER TO NOT USE INVERSE)
        A_B_i = numpy.linalg.inv(A_B)

        # Compute value of x
        x_B = numpy.dot(A_B_i, b)
        x_N = numpy.zeros(len(N))

        #Objective Vectors
        c_B = numpy.array(c)[B]
        c_N = numpy.array(c)[N]

        z_B = numpy.zeros(len(x_B))
        z_N = numpy.dot(numpy.dot(A_B_i, A_N).transpose(), c_B) - c_N

        #Check optimality
        if(numpy.all(z_N >= 0)):
            final = numpy.dot(numpy.dot(c_B, A_B_i), b)
            if debug: print("OPTIMAL!")
            if debug: print(final)
            return (B,N)

        # Choose pivot variable:
        # For Blands rule: Sort N lowest to highest, then test in order to see if coefficient is right sign
        z_N_index = 0
        entering_index = 0
        for index in N:
            entering_coefficient = z_N[z_N_index]
            if (entering_coefficient < 0):
                entering_index = index
                break
            z_N_index += 1

        # Choose Leaving Variable 
        delta_x_B = numpy.dot(A_B_i, A[:, entering_index])
        delta_x_N = numpy.zeros(len(N))

        # Check if delta_x_B <= 0, if it is, then the LP is unbounded.
        if(numpy.all(delta_x_B <= 0)):
            print("UNBOUNDED LP")
            return

        # Calculate the ratio of x_B / delta_x_B
        ratios = numpy.divide(x_B, delta_x_B)
        leaving_index = numpy.where(ratios > 0, ratios, numpy.inf).argmin() # smallest nonzero value

        leaving_index = B[leaving_index]
        # Updating B and N
        if debug: print("Entering Variable: ", entering_index)
        if debug: print("Leaving Variable: ", leaving_index)

        B.remove(leaving_index)
        B.append(entering_index)
        B.sort()

        N.remove(entering_index)
        N.append(leaving_index)
        N.sort()
        
        iteration += 1

    return

def dual_simplex(A,b,c,B,N):
    if debug: print("Starting DualSimplex")
    iteration = 1
    while True:
        # Compute our required variables based on B and N
        if debug: print("Iteration: ", iteration)
        if debug: print("Current Basis: ", B)
        if debug: print("Current Non Basis: ", N)

        # First let's make our basis matrix A_b (columns of A corresponding to our basis variables)
        A_B = A[:, B[0]]
        for basis_index in B[1:]:
            A_column = A[:, basis_index]
            A_B = numpy.column_stack((A_B, A_column))

        # Non basis matrix
        A_N = A[:, N[0]]
        for non_basis_index in N[1:]:
            A_column = A[:, non_basis_index]
            A_N = numpy.column_stack((A_N, A_column))
        
        # Take its inverse (TODO CHANGE LATER TO NOT USE INVERSE)
        A_B_i = numpy.linalg.inv(A_B)

        # Compute value of x
        x_B = numpy.dot(A_B_i, b)
        x_N = numpy.zeros(len(N))

        #Objective Vectors
        c_B = numpy.array(c)[B]
        c_N = numpy.array(c)[N]

        z_B = numpy.zeros(len(x_B))
        z_N = numpy.dot(numpy.dot(A_B_i, A_N).transpose(), c_B) - c_N

        if(numpy.all(z_N < 0)):
            print("NOT DUAL FEASIBLE")
            return

        #Check optimality
        if(numpy.all(x_B >= 0)):
            final = numpy.dot(numpy.dot(c_B, A_B_i), b)
            if debug: print("OPTIMAL!")
            if debug: print(final)
            return (B,N)

        # Choose pivot variable:
        # For Blands rule: Sort B lowest to highest, then test in order to see if coefficient is right sign
        x_B_index = 0
        leaving_index = 0
        for index in B:
            leaving_coefficient = x_B[x_B_index]
            if (leaving_coefficient < 0):
                leaving_index = index
                break
            x_B_index += 1

        # Choose Entering Variable
        # TODO MUST FIND U VECTOR
        u = numpy.zeros(len(z_B))
        index = 0
        for basic_index in B:
            if(index==basic_index):
                u[index] = 1
            else:
                u[index] = 0
        delta_z_N = -numpy.dot(numpy.dot(A_N.transpose(), numpy.linalg.inv(A_B.transpose())), u)
        delta_z_B = numpy.zeros(len(B))
        print(u)
        print(delta_z_N)
        return
        # Check if delta_z_B <= 0, if it is, then the LP is unbounded.
        '''
        if(numpy.all(delta_z_N >= 0)):
            print("UNBOUNDED LP")
            return
        '''
        

        # Calculate the ratio of x_B / delta_x_B
        ratios = numpy.divide(z_N, delta_z_N)
        entering_index = numpy.where(ratios > 0, ratios, numpy.inf).argmin() # smallest nonzero value

        entering_index = N[entering_index]
        # Updating B and N
        if debug: print("Entering Variable: ", entering_index)
        if debug: print("Leaving Variable: ", leaving_index)

        B.remove(leaving_index)
        B.append(entering_index)
        B.sort()

        N.remove(entering_index)
        N.append(leaving_index)
        N.sort()
        
        iteration += 1
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

    # LP is primal feasible
    if(numpy.all(b >= 0)):
        primal_simplex(A,b,c,B,N)
    elif(numpy.all(c <= 0)):
        dual_simplex(A,b,c,B,N)
    else:
        (B,N) = dual_simplex(A,b,numpy.zeros(len(c)),B,N)
        primal_simplex(A,b,c,B,N)


if __name__ == "__main__":
    main()

