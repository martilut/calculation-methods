import math
import numpy
import scipy


def simpson(a, b):
    param = (b - a) / 6.0
    sum1 = func(a)
    sum2 = 4 * func((a + b) / 2.0)
    sum3 = func(b)

    return param * (sum1 + sum2 + sum3)


def formula38(a, b):
    param = (b - a)
    h = (b - a) / 3.0

    sum1 = 1.0 / 8.0 * func(a)
    sum2 = 3.0 / 8.0 * func(a + h)
    sum3 = 3.0 / 8.0 * func(a + 2 * h)
    sum4 = 1.0 / 8.0 * func(b)

    return param * (sum1 + sum2 + sum3 + sum4)


def res(func, a, b):
    exact = scipy.integrate.quad(func, a, b)
    exact_value = exact[0]
    print(f"\nТочное значение: {exact_value}")
    leftRect = (b - a) * func(a)
    leftRectError = abs(exact_value - leftRect)
    print(f"КФ левого прямоугольника: {leftRect}, ошибка: {leftRectError}")
    right_rect = (b - a) * func(b)
    right_rect_err = abs(exact_value - right_rect)
    print(f"КФ правого прямоугольника: {right_rect}, ошибка: {right_rect_err}")
    center = (b - a) / 2.0
    center_rect = (b - a) * func(center)
    center_rect_err = abs(exact_value - center_rect)
    print(f"КФ среднего прямоугольника: {center_rect}, ошибка: {center_rect_err}")
    valSum = func(a) + func(b)
    trapezoid = (b - a) / 2.0 * valSum
    trapezoid_err = abs(exact_value - trapezoid)
    print(f"КФ трапеции: {trapezoid}, ошибка: {trapezoid_err}")
    simpson = simpson(a, b)
    simpson_err = abs(exact_value - simpson)
    print(f"КФ Симпсона: {simpson}, ошибка: {simpson_err}")
    formula38 = formula38(a, b)
    formula38_err = abs(exact_value - formula38)
    print(f"КФ 3/8: {formula38}, ошибка: {formula38_err}")


if __name__ == '__main__':
    print('Задача приближенного вычисления определенного интеграла')
    print('f(x) = e^(x) * sqrt(1 - x)')

    func = lambda x: math.exp(x) * math.sqrt(1 - x)
    fx = lambda x: math.exp(x)
    px = lambda x: math.sqrt(1 - x)

    print('Введите [A, B]')
    a = float(input('Введите A: '))
    b = float(input('Введите B: '))
