import math

import scipy

from root_approximator.root_approximator import secantMethod, separateRoots


def legendre_polynoms(degree):
    polynoms = ['1', 'x']
    for n in range(2, degree + 1):
        c1 = (2 * n - 1) / float(n)
        c2 = (n - 1) / float(n)
        prev1 = polynoms[n - 1]
        prev2 = polynoms[n - 2]
        curr = '(%s*%s*x-%s*%s)' % (c1, prev1, c2, prev2)
        polynoms.append(curr)
    return polynoms, polynoms[degree]


def get_roots(polynom):
    equation = lambda x: eval(polynom)
    intervals = separateRoots(equation, -1., 1., 1000)[0]
    roots = []
    for interval in intervals:
        root = secantMethod(equation, interval[0], interval[1], 0.000000000000001)
        roots.append(root)
    return roots


def get_coefficient(nodes, nodes_count, k, polynomSquare):
    node = nodes[k - 1]
    numerator = 2 * (1 - (node ** 2))
    denominator = (nodes_count ** 2) * polynomSquare(node)

    return numerator / denominator


def get_coefficients(polynoms, nodes, nodes_count):
    polynom = lambda x: eval(polynoms[nodes_count - 1])
    polynom_square = lambda x: polynom(x) * polynom(x)
    coefficients = []
    for k in range(1, nodes_count + 1):
        coef = get_coefficient(nodes, nodes_count, k, polynom_square)
        coefficients.append(coef)
    return coefficients


def gauss_formula(func, nodes_count, a=-1., b=1.):
    polynoms, polynom_degree = legendre_polynoms(nodes_count)
    nodes = get_roots(polynom_degree)
    coefficients = get_coefficients(polynoms, nodes, nodes_count)
    integral = 0
    for k in range(1, nodes_count + 1):
        ck = coefficients[k - 1]
        ak = 0.5 * (b - a) * ck
        tk = nodes[k - 1]
        xk = 0.5 * (b - a) * tk + 0.5 * (a + b)
        integral = integral + ak * func(xk)
    return nodes, coefficients, integral


if __name__ == '__main__':
    print('Задача вычисления интегралов при помощи КФ Гаусса')
    print()
    print('Вариант 6')
    print('[A, B] = [0, 1]')
    print('f(x) = x * ln(1 + x)')
    print()

    print("--------------------------------------------")
    print('Узлы и коэффициенты КФ Гаусса при N = 1...8')
    print("--------------------------------------------")
    print()

    for n in range(1, 9):
        nodes, coefficients, integral = gauss_formula(lambda x: 1, n)
        print('N =', n)
        for i in range(0, len(nodes)):
            print(f"Узел: {nodes[i]}, Коэффициент: {coefficients[i]}")
        print()

    print("--------------------------------------------")
    print('Проверка на многочленах')
    print('[A, B] = [-1, 1]')
    print("--------------------------------------------")
    print()

    for n in [3, 4, 5]:
        degree = 2 * n - 1
        polynom_str = 'x^%s' % degree
        polynom = lambda x: x ** degree

        _, _, gauss_value = gauss_formula(polynom, n)
        exact_value = scipy.integrate.quad(polynom, -1, 1)[0]
        gauss_error = abs(exact_value - gauss_value)

        print('N:', n)
        print(f"Polynom: {polynom_str}")
        print(f"Exact value: {exact_value}")
        print(f"Gauss value: {gauss_value}")
        print(f"Abs error: {gauss_error}")
        print()

    print("--------------------------------------------")
    print('Вычисление интеграла функции с помощью КФ Гаусса')
    print('f(x) = x * ln(1 + x)')
    print("--------------------------------------------")
    print()

    func = lambda x: x * math.log(1 + x)

    print('Введите [A, B]')
    a = float(input('Введите A: '))
    b = float(input('Введите B: '))
    print()

    to_continue = True
    while to_continue:
        nodes_list = input("Введите N_1, N_2, N_3 узлов для вычисления интеграла с помощью КФ Гаусса: ").split(" ")
        for nodes_count in nodes_list:
            exact_value = scipy.integrate.quad(func, a, b)[0]
            nodes, coefficients, gauss_value = gauss_formula(func, int(nodes_count), a, b)
            gauss_error = abs(exact_value - gauss_value)

            print('N:', nodes_count)
            print(f"Nodes: {nodes}")
            print(f"Coefficients: {coefficients}, sum: {sum(coefficients)}")
            print(f"Exact value: {exact_value}")
            print(f"Gauss value: {gauss_value}")
            print(f"Abs error: {gauss_error}")
            print()
        to_continue_letter = input("Ввести новые значения <a, b>? y/n: ")
        if to_continue_letter == 'n':
            to_continue = False