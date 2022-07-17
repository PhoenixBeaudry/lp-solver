import sys
import numpy
from fractions import Fraction
numpy.seterr(divide='ignore', invalid='ignore')
debug = False
debug_iterations = 364
total_iterations = 1
precision = 10

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
    c = (lp[0] + [0]*len(b))
    #Add the identity basis (slack variables)
    A = numpy.concatenate((A,numpy.identity(A.shape[0])), axis=1).astype(float)
    b = numpy.array(b).astype(float)
    c = numpy.array(c).astype(float)
    return (A,b,c,B,N)

def primal_simplex(A,b,c,B,N):
    global total_iterations
    if debug: print("Starting PrimalSimplex")
    while True:
        if debug and total_iterations == debug_iterations: return
        # Compute our required variables based on B and N
        if debug: print("Primal Iteration: ", total_iterations)

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
        
        # Take its inverse
        A_B_i = numpy.around(numpy.linalg.inv(A_B), precision)

        # Compute value of x
        x_B = numpy.around(numpy.dot(A_B_i, b), precision)
        x_N = numpy.zeros(len(N))

        
        #Objective Vectors
        c_B = numpy.array(c)[B]
        c_N = numpy.array(c)[N]

        z_B = numpy.zeros(len(x_B))
        z_N = numpy.around(numpy.dot(numpy.dot(A_B_i, A_N).transpose(), c_B) - c_N, precision)

        if debug: print("Current z_N")
        if debug: print(z_N)


        #Check optimality
        if(numpy.all(z_N >= 0)):
            # Find final value of variables 0...n
            final_variables = []
            for i in range(0, len(N)):
                try:
                    value = x_N[N.index(i)]
                except:
                    value = x_B[B.index(i)]
                final_variables.append(value)
            # Find final optimal value
            final = numpy.dot(numpy.dot(c_B, A_B_i), b)
            # Output
            out = ["optimal", final, final_variables]
            return (B,N,out)


        # Choose pivot variable:
        # Largest coefficient
        
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
        delta_x_B = numpy.around(numpy.dot(A_B_i, A[:, entering_index]), precision)
        delta_x_N = numpy.zeros(len(N))

        # Check if delta_x_B <= 0, if it is, then the LP is unbounded.
        if(numpy.all(delta_x_B <= 0)):
            out = ["unbounded"]
            return (B,N,out)

        if debug: print("Current x_B: ")
        if debug: print(x_B)
        if debug: print("Current delta_x_B: ")
        if debug: print(delta_x_B)
        
        delta_x_B = numpy.ma.array(delta_x_B, mask=(delta_x_B <= 0))
        if debug: print("Current masked delta_x_B: ")
        if debug: print(delta_x_B)
        # Calculate the ratio of x_B / delta_x_B
        ratios = numpy.around(numpy.divide(x_B, delta_x_B), precision)
        ratios =  numpy.ma.array(numpy.nan_to_num(ratios, nan=numpy.nan, posinf=numpy.nan, neginf=numpy.nan))

        if debug: print("Current Basis: ")
        if debug: print(B)
        if debug: print("Current Non Basis: ")
        if debug: print(N)
        if debug: print("Ratios to choose leaving variable: ")
        if debug: print(ratios)
        
        leaving_index = numpy.ma.argmin(ratios)

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
        
        total_iterations += 1

    return

def dual_simplex(A,b,c,B,N):
    global total_iterations
    if debug: print("Starting DualSimplex")
    while True:
        if debug and total_iterations == debug_iterations: return
        # Compute our required variables based on B and N
        if debug: print("Dual Iteration: ", total_iterations)

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
        
        # Take its inverse
        A_B_i = numpy.around(numpy.linalg.inv(A_B), precision)

        # Compute value of x
        x_B = numpy.around(numpy.dot(A_B_i, b), precision)
        x_N = numpy.zeros(len(N))

        #Objective Vectors
        c_B = numpy.array(c)[B]
        c_N = numpy.array(c)[N]

        z_B = numpy.zeros(len(x_B))
        z_N = numpy.around(numpy.dot(numpy.dot(A_B_i, A_N).transpose(), c_B) - c_N, precision)

        #Check optimality
        if(numpy.all(x_B >= 0)):
            # Find final value of variables 0...n
            final_variables = []
            for i in range(0, len(N)):
                try:
                    value = x_N[N.index(i)]
                except:
                    value = x_B[B.index(i)]
                final_variables.append(value)
            
            final = numpy.dot(numpy.dot(c_B, A_B_i), b)


            out = ["optimal", final, final_variables]
            return (B,N,out)

        # Choose leaving pivot variable:
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
        # Find U vector
        u = numpy.zeros(len(z_B))
        index = 0
        for basic_index in B:
            if(basic_index==leaving_index):
                u[index] = 1
            else:
                u[index] = 0
            index += 1

        # Compute delta_z_N and delta_z_B
        v = numpy.around(numpy.linalg.solve(A_B.transpose(), u), precision)
        delta_z_N = -numpy.around(numpy.dot(A_N.transpose(), v), precision)
        delta_z_B = numpy.zeros(len(B))

        if debug: print("Current z_N: ")
        if debug: print(z_N)
        if debug: print("Current delta_z_N: ")
        if debug: print(delta_z_N)

        # Check if delta_z_B <= 0, if it is, then the LP is unbounded.
        if(numpy.all(delta_z_N <= 0)):
            out = ["infeasible"]
            return (B,N,out)

        delta_z_N = numpy.ma.array(delta_z_N, mask=(delta_z_N <= 0))

        # Calculate the ratio of z_N and delta_z_N and mask invalid values
        ratios = numpy.around(numpy.divide(z_N, delta_z_N), precision)
        ratios =  numpy.ma.array(numpy.nan_to_num(ratios, nan=numpy.nan, posinf=numpy.nan, neginf=numpy.nan))
        
        if debug: print("Current Non Basis: ")
        if debug: print(N)
        if debug: print("Current Basis: ")
        if debug: print(B)
        if debug: print("Ratios to choose entering variable: ")
        if debug: print(ratios)

        entering_index = numpy.ma.argmin(ratios)
        

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
        
        total_iterations += 1
    return


def solve_lp(lp):
    # Convert LP into dictionary matrix form
    (A,b,c,B,N) = standard_form_to_eq(lp)

    if debug: print("Current LP: ")
    if debug: print("A: ")
    if debug: print(A)
    if debug: print("b: ")
    if debug: print(b)
    if debug: print("c: ")
    if debug: print(c)

    # LP is primal feasible
    if(numpy.all(b >= 0)):
        if debug: print("LP is Primal Feasible")
        (B,N,out) = primal_simplex(A,b,c,B,N)
        print(out[0])
        try:
            print('{:.7g}'.format(out[1]))
            output_variables = ' '.join(['{:.7g}'.format(variable) for variable in out[2]])
            print(output_variables)
        except:
            return
        return
    elif(numpy.all(c <= 0)):
        if debug: print("LP is Dual Feasible")
        if(numpy.all(c == 0)):
            # Weird degenerate 0 case run with -1
            (B,N,out) = dual_simplex(A,b,numpy.full(len(c), -1),B,N)
        else:
            (B,N,out) = dual_simplex(A,b,c,B,N)
        print(out[0])
        try:
            print('{:.7g}'.format(out[1]))
            output_variables = ' '.join(['{:.7g}'.format(variable) for variable in out[2]])
            print(output_variables)
        except:
            return
        return
    else:
        if debug: print("LP is neither Primal nor Dual Feasible")
        (B,N,out) = dual_simplex(A,b,numpy.full(len(c), -1),B,N)
        if(out[0] == "infeasible" or out[0] == "unbounded"):
            print(out[0])
            return
        (B,N,out) = primal_simplex(A,b,c,B,N)
        print(out[0])
        try:
            print('{:.7g}'.format(out[1]))
            output_variables = ' '.join(['{:.7g}'.format(variable) for variable in out[2]])
            print(output_variables)
        except:
            return
        return

def main():
    # Read in LP from stdin
    lp = []
    for line in sys.stdin:
        coeff_line = []
        for coeff in line.split():
            coeff_line.append(float(coeff))
        lp.append(coeff_line)
    solve_lp(lp)
    return


if __name__ == "__main__":
    main()

