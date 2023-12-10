import math

import numpy
import scipy

from root_approximator.root_approximator import separateRoots, bisectionMethod


def simpson(func, a, b):
    param = (b - a) / 6.0
    sum1 = func(a)
    sum2 = 4 * func((a + b) / 2.0)
    sum3 = func(b)

    return param * (sum1 + sum2 + sum3)


def calculate_moments(weight, nodes_count, a, b):
    moments = []
    for k in range(0, nodes_count):
        k_func = lambda x: weight(x) * (x ** k)
        moment = scipy.integrate.quad(k_func, a, b)[0]
        moments.append(moment)
    print(f"Моменты: {moments}")
    return moments


def calculate_coefficients(values, moments):
    rows = []
    for k in range(0, len(values)):
        line = [x ** k for x in values]
        rows.append(line)
    matrix = numpy.array(rows)
    vector = numpy.array(moments)
    coefficients = list(numpy.linalg.solve(matrix, vector))
    print(f"Коэффициенты КФ: {coefficients}")
    return coefficients


def orthogonal_polynom(weight, nodes_count, a, b):
    moments = calculate_moments(weight, nodes_count * 2, a, b)
    rows = []
    for k in range(0, nodes_count):
        line = []
        for i in range(k, nodes_count + k):
            line.append(moments[i])
        rows.append(line)
    matrix = numpy.array(rows)
    part = moments[nodes_count:(2 * nodes_count)]
    vector = numpy.array([-t for t in part])
    coefficients = numpy.linalg.solve(matrix, vector)
    coefficients = list(reversed(coefficients))
    return coefficients


def calculate_nodes(weight, nodes_count, a, b):
    orthogonal_coefficients = orthogonal_polynom(weight, nodes_count, a, b)
    equation_str = f"x ** {nodes_count}"
    for i in range(0, nodes_count):
        coef = orthogonal_coefficients[i]
        degree = nodes_count - i - 1
        equation_str = f"{equation_str} + {coef} * x ** {degree}"
    print(f"Ортогональный многочлен: {equation_str}")
    equation_lambda = lambda x: eval(equation_str)
    intervals = separateRoots(equation_lambda, a, b, 1000)[0]
    roots = []
    for interval in intervals:
        root = bisectionMethod(equation_lambda, interval[0], interval[1], 0.000000001)
        roots.append(root)
    print(f"Корни ортогонального многочлена: {roots}")
    return roots


def calculate_integral(func, nodes, weight, a, b):
    nodes_count = len(nodes)
    moments = calculate_moments(weight, nodes_count, a, b)
    coefficients = calculate_coefficients(nodes, moments)
    integral = 0
    for i in range(0, nodes_count):
        integral += coefficients[i] * func(nodes[i])
    return integral


def create_nodes(nodes_count, a, b):
    nodes = []
    step = (b - a) / nodes_count
    for i in range(0, nodes_count):
        nodes.append(a + i * step)
    return nodes


if __name__ == "__main__":
    print("Задача приближенного вычисления интеграла при помощи ИКФ и КФ НАСТ")
    print("Вариант 6")
    print("f(x): sin(x), p(x): -x * ln(x)")
    print("[a, b] = [0, 1]")

    f_x = lambda x: math.sin(x)
    p_x = lambda x: -x * math.log(x)
    full_func = lambda x: f_x(x) * p_x(x)

    a = float(input('Введите A: '))
    b = float(input('Введите B: '))

    nodes_count = int(input('Введите число узлов: '))

    IQF_polynom = lambda x: x ** (nodes_count - 1)
    IQF_polynom_weight = lambda x: 1
    IQF_polynom_exact = scipy.integrate.quad(IQF_polynom, a, b)[0]

    QFNAST_polynom = lambda x: x ** (2 * nodes_count - 1)
    QFNAST_polynom_weight = lambda x: 1
    QFNAST_polynom_exact = scipy.integrate.quad(QFNAST_polynom, a, b)[0]

    exact_value = scipy.integrate.quad(full_func, a, b)[0]
    print("-------------------------------")
    print(f"Точное значение: {exact_value}")
    print("-------------------------------")

    IQF_nodes = create_nodes(nodes_count, a, b)
    IQF_integral = calculate_integral(f_x, IQF_nodes, p_x, a, b)
    IQF_error = abs(exact_value - IQF_integral)
    print(f"Узлы IQF: {IQF_nodes}")
    print(f"Приближенное значение: {IQF_integral}")
    print(f"Абсолютная погрешность: {IQF_error}")

    print(f"\n------ Проверка на многочлене: x ** {(nodes_count - 1)} ------")
    IQF_polynom_integral = calculate_integral(IQF_polynom, IQF_nodes, IQF_polynom_weight, a, b)
    print(f"Точное значение: {IQF_polynom_exact}")
    print(f"Приближенное значение: {IQF_polynom_integral}")
    print(f"Абсолютная погрешность: {abs(IQF_polynom_exact - IQF_polynom_integral)}")
    print("-------------------------------")

    QFNAST_nodes = calculate_nodes(p_x, nodes_count, a, b)
    QFNAST_integral = calculate_integral(f_x, QFNAST_nodes, p_x, a, b)
    QFNAST_error = abs(exact_value - QFNAST_integral)
    print(f"Узлы КФ НАСТ: {QFNAST_nodes}")
    print(f"Приближенное значение: {QFNAST_integral}")
    print(f"Абсолютная погрешность: {QFNAST_error}")

    print(f"\n------ Проверка на многочлене: x ** {(2 * nodes_count - 1)} ------")
    QFNAST_polynom_nodes = calculate_nodes(QFNAST_polynom_weight, nodes_count, a, b)
    QFNAST_polynom_integral = calculate_integral(QFNAST_polynom, QFNAST_polynom_nodes, QFNAST_polynom_weight, a, b)
    print(f"Точное значение: {QFNAST_polynom_exact}")
    print(f"Приближенное значение: {QFNAST_polynom_integral}")
    print(f"Абсолютная погрешность: {abs(QFNAST_polynom_exact - QFNAST_polynom_integral)}")
    print("-------------------------------")
