
import numpy
import numpy as np
import math

S = 0
eps_cur = 0
a = 0
b = 3
c = 0
d = 1

def m1 (y):
    return np.sin(np.pi * y) ** 2
def m2 (x):
    return np.cosh(x * x - 3 * x) - 1
def m3 (y):
    return 0
def m4 (x):
    return 0

def delta_u (x, y):
    return np.cosh(x - y)
def _delta_u (j, i):
    return delta_u(a + j * h2, c + i * k2)

def initial_approximation():

    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                V[i][j] = m2(a + j * h2)
            if i == m:
                V[i][j] = m4(a + j * h2)
            if j == 0:
                V[i][j] = m1(c + i * k2)
            if j == n:
                V[i][j] = m3(c + i * k2)

    ##вектор F
    for i in range(1, m):
        for j in range(1, n):
            sum = 0
            if j == 1:
                sum += (-h) * m1(c + i * k2)
            else:
                if j == (n - 1):
                    sum += (-h) * m3(c + i * k2)
            if i == 1:
                sum += (-k) * m2(a + j * h2)
            else:
                if i == (m - 1):
                    sum += (-k) * m4(a + j * h2)
            F[i][j] = -sum - _delta_u(j, i)


def Matrix_multiplication(X):
    Res = np.zeros((m + 1,  n + 1))
    for i in range(1, m):
        for j in range(1, n):
            if (i != 1 and i != m - 1):
                if (j != 1 and j != n - 1):
                    mult = -(k * X[i - 1][j] + h * X[i][j - 1] + a2 * X[i][j] + h * X[i][j + 1] + k * X[i + 1][j])
                else:
                    if (j == 1):
                        mult = -(k * (X[i - 1][j] + X[i + 1][j]) + a2 * X[i][j] + h * X[i][j + 1])
                    else:
                        if (j == n - 1):
                            mult = -(k * (X[i - 1][j] + X[i + 1][j]) + h * X[i][j - 1] + a2 * X[i][j])
            else:
                if (i == 1):
                    if (j == 1):
                        mult = -(a2 * X[i][j] + h * X[i][j + 1] + k * X[i + 1][j])
                    else:
                        if (j != n - 1):
                            mult = -(h * X[i][j - 1] + a2 * X[i][j] + h * X[i][j + 1] + k * X[i + 1][j])
                        else:
                            if (j == n - 1):
                                mult = -(h * X[i][j - 1] + a2 * X[i][j] + k * X[i + 1][j])

                else:
                    if (i == m - 1):
                        if (j == 1):
                            mult = -(k * X[i - 1][j] + a2 * X[i][j] + h * X[i][j + 1])
                        else:
                            if (j != n - 1):
                                mult = -(k * X[i - 1][j] + h * X[i][j - 1] + a2 * X[i][j] + h * X[i][j + 1])
                            else:
                                if (j == n - 1):
                                    mult = -(k * X[i - 1][j] + h * X[i][j - 1] + a2 * X[i][j])

            Res[i][j] = mult
    return Res
def Discrepancy(Ax, b):
    Res = np.zeros((m + 1, n + 1))
    for i in range(1, m):
        for j in range(1, n):
            Res[i][j] = Ax[i][j] - b[i][j]
    return Res
def Ao(r, h):
    s = 0
    b = 0
    Ah = Matrix_multiplication(h)
    for i in range(1, m):
        for j in range(1, n):
            s += r[i][j] * h[i][j]
            b += Ah[i][j] * h[i][j]
    A = - s / b
    return A
def Bo(h, r):
    l = 0
    s = 0
    Ah = Matrix_multiplication(h)
    for i in range(1, m):
        for j in range(1, n):
            s += Ah[i][j] * r[i][j]
            l += Ah[i][j] * h[i][j]
    B = s / l
    return B


def Conjugate_gradient( N, M, N_max, EPS):
    global  n, m, S, h2, k2, h, k, a2, V, U, F, Nmax, eps, eps_max, eps_tesk, max_r, max_ro, xm, ym
    n = N; m = M; Nmax = N_max; eps = EPS
    h2 = (b - a) / n
    k2 = (d - c) / m
    h = (-1) / (h2 * h2)
    k = (-1) / (k2 * k2)
    a2 = (-2) * (h + k)
    V = np.zeros((m + 1, n + 1))
    U = np.zeros((m + 1, n + 1))
    F = np.zeros((m + 1, n + 1))
    initial_approximation()
    H = np.zeros((m + 1, n + 1))
    h_old = np.zeros((m + 1, n + 1))
    r = Discrepancy(Matrix_multiplication(V), F)
    max_ro = 0
    for i in range(1, m):
        for j in range(1, n):
            if np.abs(max_ro) < np.abs(r[i][j]):  # норма невязки на начальном приближении
                max_ro = np.abs(r[i][j])
    eps_max = 0
    for i in range(1, m):
        for j in range(1, n):
            H[i][j] = - r[i][j]
    A = Ao(r, H)
    for i in range(1, m):
        for j in range(1, n):
            v_old = V[i][j]
            v_new = V[i][j] + A * H[i][j]
            eps_cur = np.abs(v_old - v_new)
            if eps_cur > eps_max:
                eps_max = eps_cur
            V[i][j] = v_new
    S = 1
    while (S < Nmax - 1):
        eps_max = 0
        h_old = H
        r = Discrepancy(Matrix_multiplication(V), F)
        B = Bo(h_old, r)
        for i in range(1, m):
            for j in range(1, n):
                H[i][j] = - r[i][j] + B * h_old[i][j]
        A = Ao(r, H)
        for i in range(1, m):
            for j in range(1, n):
                v_old = V[i][j]
                v_new = V[i][j] + A * H[i][j]
                eps_cur = np.abs(v_old - v_new)
                if eps_cur > eps_max:
                    eps_max = eps_cur
                V[i][j] = v_new
        if eps_max <= eps:
            break
        S = S + 1
    r = Discrepancy(Matrix_multiplication(V), F)
    max_r = 0
    xm = 0
    ym = 0
    for i in range(1, m):
        for j in range(1, n):
            if np.abs(max_r) < np.abs(r[i][j]):  # норма невязки на начальном приближении
                max_r = np.abs(r[i][j])
                xm = a + j * h2
                ym = c + i * k2
    print(S)
    return V