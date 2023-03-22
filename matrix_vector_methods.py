


def matrix_zeros(x, y):
    matrix = []
    for _ in range(y):
        row = []
        for _ in range(x):
            row.append(int(0))
        matrix.append(row)
    return matrix


def diagonal_to_square_matrix(vec):
    m = matrix_zeros(len(vec), len(vec))
    for i in range(len(vec)):
        m[i][i] = vec[i]
    return m


def vec_ones(N):
    v = []
    for _ in range(N):
        v.append(1.0)
    return v


def copy_matrix(_matrix):
    copy = []
    for row in _matrix:
        new_row = []
        for i in row:
            new_row.append(i)
        copy.append(new_row)
    return copy


def copy_vector(vector):
    vec = []
    for i in vector:
        vec.append(i)
    return vec


def fill_vector(N, z):
    vec = []
    for _ in range(N):
        vec.append(z)
    return vec


def substract_two_vectors(A, b):
    diff = copy_vector(A)
    for i in range(len(diff)):
        diff[i] -= b[i]
    return diff


def dot_product(A, _b, N):
    matrix_a = copy_matrix(A)
    vector_b = copy_vector(_b)
    n = len(matrix_a[0])
    c = fill_vector(N, 0)

    for i in range(N):
        for l in range(n):
            c[i] += matrix_a[i][l] * vector_b[l]
    return c


def normRes(res, N):
    norm = 0
    for i in range(N):
        norm += pow(res[i], 2)
    return pow(norm, 0.5)

