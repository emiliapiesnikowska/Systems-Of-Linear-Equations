import math
import time
from matrix_vector_methods import *

from matplotlib import pyplot


def Matrix(N, a1, a2, a3):
    A = []
    for i in range(N):
        row = []
        for j in range(N):
            if i == j:
                row.append(a1)
            elif i - 1 <= j <= i + 1:
                row.append(int(a2))
            elif i - 2 <= j <= i + 2:
                row.append(int(a3))
            else:
                row.append(int(0))
        A.append(row)
    return A


def Matrix_A(N, a1, a2, a3):
    return Matrix(N, a1, a2, a3)


def Vector_b(N, f):
    b = []
    for i in range(N):
        b.append(math.sin(i * (f + 1)))
    return b


def jacobi(A, b, N):
    # init
    matrix_a = copy_matrix(A)
    vector_b = copy_vector(b)
    iterations = 0
    prev_x = fill_vector(N, 0)
    for i in range(N):
        prev_x.append(1)

    x = fill_vector(N, 0)

    # perform
    while True:
        for i in range(N):
            s = 0
            for j in range(N):
                if i != j:
                    s += matrix_a[i][j] * prev_x[j]
            x[i] = (vector_b[i] - s)/matrix_a[i][i]
        iterations += 1
        prev_x = copy_vector(x)
        res = substract_two_vectors(dot_product(matrix_a, x, N), vector_b)

        if normRes(res, N) <= pow(10, -9):
            break

    # results
    print('number of iterations: ', iterations)


def gauss_seidel(A, b, N):
    matrix_a = copy_matrix(A)
    vector_b = copy_vector(b)
    iterations = 0
    prev_x = fill_vector(N, 1)
    x = fill_vector(N, 0)

    while True:
        for i in range(N):
            s = 0
            for j in range(i):
                s += matrix_a[i][j] * x[j]
            for j in range(N):
                if j >= i+1:
                    s += matrix_a[i][j] * prev_x[j]
            x[i] = (vector_b[i] - s) / matrix_a[i][i]

        iterations += 1
        prev_x = copy_vector(x)
        res = substract_two_vectors(dot_product(matrix_a, x, N), vector_b)
        if normRes(res, N) <= pow(10, -9):
            break

    print('number of iterations: ', iterations)


def LU(A, b, N):
    matrix_A = copy_matrix(A)
    matrix_L = diagonal_to_square_matrix(vec_ones(N))
    matrix_U = matrix_zeros(N, N)

    vector_b = copy_vector(b)
    vector_x = fill_vector(N, 1)

    for j in range(N):
        for i in range(j + 1):
            matrix_U[i][j] += matrix_A[i][j]
            for k in range(i):
                matrix_U[i][j] -= matrix_L[i][k] * matrix_U[k][j]

        for i in range(j + 1, N):
            for k in range(j):
                matrix_L[i][j] -= matrix_L[i][k] * matrix_U[k][j]
            matrix_L[i][j] += matrix_A[i][j]
            matrix_L[i][j] /= matrix_U[j][j]

    vector_y = fill_vector(N, 0)
    for i in range(N):
        s = 0
        for j in range(i):
            s += matrix_L[i][j] * vector_y[j]

        vector_y[i] = (vector_b[i] - s) / matrix_L[i][i]

    for i in range(N - 1, -1, -1):
        s = 0
        for j in range(i + 1, N):
            s += matrix_U[i][j] * vector_x[j]
        vector_x[i] = (vector_y[i] - s) / matrix_U[i][i]

    res = substract_two_vectors(dot_product(matrix_A, vector_x,N), vector_b)
    print("Residuum norm:", normRes(res, N))


if __name__ == '__main__':

    # moj numer indeksu to 180112
    # zatem e = 1 => a1 = 5+1
    # a2 = a3 = -1
    # N = 9cd, gdzie c=1, d=2
    # b jest wektorem o długości N, którego n−ty element ma wartość sin(n ·(f + 1)), gdzie f = 0
    # sin(n)
    e = 1
    f = 0
    N = 912
    a1 = 5+e
    a2 = -1
    a3 = -1

    # ZADANIE A
    A = Matrix_A(N, a1, a2, a3)
    b = Vector_b(N, f)

    # ZADANIE B
    print("Jacobi")
    time1 = time.time()
    jacobi(A, b, N)
    time2 = time.time()
    print('time:', time2 - time1)

    print("Gauss-Seidel")
    time1 = time.time()
    gauss_seidel(A, b, N)
    time2 = time.time()
    print('time:', time2 - time1)
    print()
    # ZADANIE C
    print("Zadanie C")
    a1C = 3
    A = Matrix_A(N, a1C, a2, a3)
    #  print("Jacobi"); time1 = time.time(); jacobi(A, b, N); time2 = time.time(); print('time:', time2 - time1)
    #  print("Gauss-Seidel"); time1 = time.time(); gauss_seidel(A, b, N); time2 = time.time(); print('time:', time2 - time1)

    # ZADANIE D
    print("LU")
    time1 = time.time()
    LU(A, b, N)
    time2 = time.time()
    print('time:', time2 - time1)


    #ZADANIE E

    N2 = [100, 500, 1000, 2000]
    time_jacobi = []
    time_gs = []
    time_lu = []
    for n in N2:
        print("Size:", n)
        matrix_A = Matrix_A(n, a1, a2, a3)
        vector_b = Vector_b(n, f)

        time1 = time.time()
        jacobi(matrix_A, vector_b, n)
        time_jacobi.append(time.time() - time1)

        time1 = time.time()
        gauss_seidel(matrix_A, vector_b, n)
        time_gs.append(time.time() - time1)

        time1 = time.time()
        LU(matrix_A, vector_b, n)
        time_lu.append(time.time() - time1)

    pyplot.plot(N2, time_gs, label="Gauss-Seidel", color="blue")
    pyplot.plot(N2, time_lu, label="LU", color="green")
    pyplot.plot(N2, time_jacobi, label="Jacobi", color="red")
    pyplot.legend()
    pyplot.grid(True)
    pyplot.ylabel('Time(s)')
    pyplot.xlabel('Number of unknowns')
    pyplot.title('Time dependence on the number of unknowns')
    pyplot.show()



