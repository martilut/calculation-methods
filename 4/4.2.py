import math

import scipy


class MethodData:
    def __init__(self, name, exact_value, value, theoretical_error):
        self.name = name
        self.exact_value = exact_value
        self.value = value
        self.abs_error = abs(exact_value - value)
        self.rel_error = abs(self.abs_error / exact_value) if exact_value != 0 else 0
        self.theoretical_error = theoretical_error


class CalculationData:
    def calculate_data(self, func, func_int, a, b, step, parts):
        exact_integral = scipy.integrate.quad(func_int, a, b)
        exact_value = exact_integral[0]
        left_rect = leftRectangle(func, a, step, parts)
        left_rect_theor_err = ((b - a) ** 2) / (2 * parts) * abs(math.exp(b))
        left_rect_data = MethodData("Left rectangle", exact_value, left_rect, left_rect_theor_err)
        right_rect = rightRectangle(func, a, step, parts)
        right_rect_theor_err = ((b - a) ** 2) / (2 * parts) * abs(math.exp(b))
        right_rect_data = MethodData("Right rectangle", exact_value, right_rect, right_rect_theor_err)
        center_rect = centerRectangle(func, a, step, parts)
        center_rect_theor_err = ((b - a) ** 3) / (24 * parts ** 2) * abs(math.exp(b))
        center_rect_data = MethodData("Center rectangle", exact_value, center_rect, center_rect_theor_err)
        trapez = trapezoid(func, a, b, step, parts)
        trapez_theor_err = ((b - a) ** 3) / (12 * parts ** 2) * abs(math.exp(b))
        trapez_data = MethodData("Trapezoid", exact_value, trapez, trapez_theor_err)
        simp = simpson(func, a, b, step, parts)
        simp_theor_err = ((b - a) ** 5) / (2880 * parts ** 4) * abs(math.exp(b))
        simp_data = MethodData("Simpson", exact_value, simp, simp_theor_err)
        return [left_rect_data, right_rect_data, center_rect_data, trapez_data, simp_data]


def print_calculation_data(stats, theor_err=False):
    for data in stats:
        print(f"---{data.name}---")
        print("Exact value: " + str(data.exact_value))
        print("Value: " + str(data.value))
        print("Absolute error: " + str(data.abs_error))
        print("Relative error: " + str(data.rel_error))
        if theor_err:
            print("Theoretical error: " + str(data.theoretical_error))
        print()


def leftRectangle(func, a, step, parts):
    res = 0
    for j in range(0, parts):
        arg = a + j * step
        res += func(arg)
    return step * res


def rightRectangle(func, a, step, parts):
    res = 0
    for j in range(0, parts):
        arg = a + (j + 1) * step
        res += func(arg)
    return step * res


def centerRectangle(func, a, step, parts):
    res = 0
    for j in range(0, parts):
        arg = a + (j + 0.5) * step
        res += func(arg)
    return step * res


def trapezoid(func, a, b, step, parts):
    res = func(a) + func(b)
    for j in range(1, parts):
        arg = a + j * step
        res += 2 * func(arg)
    return 0.5 * step * res


def simpson(func, a, b, step, parts):
    res = func(a) + func(b)
    for j in range(1, parts):
        arg = a + j * step
        res += 2 * func(arg)
    for j in range(0, parts):
        arg = a + (j + 0.5) * step
        res += 4 * func(arg)
    return (step / 6.0) * res


def refineStats(start_data, new_data, L):
    stats = []
    for i in range(0, len(start_data)):
        start_current = start_data[i]
        new_current = new_data[i]
        if start_current.name == 'Среднего пр-ка' or start_current.name == 'Трапеции':
            d = 1
        elif start_current.name == 'Симпсона':
            d = 3
        else:
            d = 0
        value_refined = ((L ** (d + 1)) * new_current.value - start_current.value) / ((L ** (d + 1)) - 1)
        value_exact = start_current.exact_value
        refinedStats = MethodData(start_current.name, value_exact, value_refined, None)
        stats.append(refinedStats)
    return stats


if __name__ == '__main__':
    print('Задача приближенного вычисления интеграла по составным квадратурным формулам')
    print('f(x) = exp(x)')
    print()

    f = lambda x: math.exp(x)
    p = lambda x: 1
    func_int = lambda x: f(x) * p(x)

    polynom_0 = lambda x: 100
    polynom_1 = lambda x: x ** 1
    polynom_2 = lambda x: x ** 2
    polynom_3 = lambda x: x ** 3

    a = float(input('Введите A: '))
    b = float(input('Введите B: '))

    print()
    parts = int(input('Введите число разбиений [A, B]: '))

    start_steps = (b - a) / float(parts)
    start_stats = CalculationData().calculate_data(f, func_int, a, b, start_steps, parts)

    print()
    print('m = ', parts)
    print('h =', start_steps)
    print_calculation_data(start_stats, True)
    print("------------------------------------------------")
    print("Проверка на полиномах")
    print("------------------------------------------------")

    print('f(x) = 100')
    polynom_0_data = CalculationData().calculate_data(polynom_0, polynom_0, a, b, start_steps, parts)
    print_calculation_data(polynom_0_data)
    print("------------------------------------------------")

    print('f(x) = x')
    polynom_1_data = CalculationData().calculate_data(polynom_1, polynom_1, a, b, start_steps, parts)
    print_calculation_data(polynom_1_data)
    print("------------------------------------------------")

    print('f(x) = x^2')
    polynom_2_data = CalculationData().calculate_data(polynom_2, polynom_2, a, b, start_steps, parts)
    print_calculation_data(polynom_2_data)
    print("------------------------------------------------")

    print('f(x) = x^3')
    polynom_3_data = CalculationData().calculate_data(polynom_3, polynom_3, a, b, start_steps, parts)
    print_calculation_data(polynom_3_data)

    print("------------------------------------------------")
    L = int(input('Введите L: '))
    print("------------------------------------------------")

    parts *= L
    steps_new = (b - a) / float(parts)
    data_new = CalculationData().calculate_data(f, func_int, a, b, steps_new, parts)
    print()
    print('h1 =', start_steps)
    print_calculation_data(start_stats)
    print("------------------------------------------------")

    print()
    print('h2 =', steps_new)
    print_calculation_data(data_new)
