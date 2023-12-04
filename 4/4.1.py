import math
import scipy


def simpson(func, a, b):
    param = (b - a) / 6.0
    sum1 = func(a)
    sum2 = 4 * func((a + b) / 2.0)
    sum3 = func(b)

    return param * (sum1 + sum2 + sum3)


def formula38(func, a, b):
    param = (b - a)
    h = (b - a) / 3.0

    sum1 = 1.0 / 8.0 * func(a)
    sum2 = 3.0 / 8.0 * func(a + h)
    sum3 = 3.0 / 8.0 * func(a + 2 * h)
    sum4 = 1.0 / 8.0 * func(b)

    return param * (sum1 + sum2 + sum3 + sum4)


def left_rect(func, a, b, exact_value):
    leftRect = (b - a) * func(a)
    leftRectError = abs(exact_value - leftRect)
    return leftRect, leftRectError


def right_rect(func, a, b, exact_value):
    right_rect = (b - a) * func(b)
    right_rect_err = abs(exact_value - right_rect)
    return right_rect, right_rect_err


def center(func, a, b, exact_value):
    center = (b - a) / 2.0
    center_rect = (b - a) * func(center)
    center_rect_err = abs(exact_value - center_rect)
    return center_rect, center_rect_err


def trapezoid(func, a, b, exact_value):
    valSum = func(a) + func(b)
    trapezoid = (b - a) / 2.0 * valSum
    trapezoid_err = abs(exact_value - trapezoid)
    return trapezoid, trapezoid_err


def calculate(func, a, b):
    exact = scipy.integrate.quad(func, a, b)
    exact_value = exact[0]
    print(f"\nТочное значение: {exact_value}")
    leftRect, leftRectError = left_rect(func, a, b, exact_value)
    print(f"КФ левого прямоугольника: {leftRect}, ошибка: {leftRectError}")
    right_rect_val, right_rect_err = right_rect(func, a, b, exact_value)
    print(f"КФ правого прямоугольника: {right_rect_val}, ошибка: {right_rect_err}")
    center_rect, center_rect_err = center(func, a, b, exact_value)
    print(f"КФ среднего прямоугольника: {center_rect}, ошибка: {center_rect_err}")
    trapezoid_val, trapezoid_err = trapezoid(func, a, b, exact_value)
    print(f"КФ трапеции: {trapezoid_val}, ошибка: {trapezoid_err}")
    simpson_val = simpson(func, a, b)
    simpson_err = abs(exact_value - simpson_val)
    print(f"КФ Симпсона: {simpson_val}, ошибка: {simpson_err}")
    formula38_val = formula38(func, a, b)
    formula38_err = abs(exact_value - formula38_val)
    print(f"КФ 3/8: {formula38_val}, ошибка: {formula38_err}")


if __name__ == '__main__':
    print('Задача приближенного вычисления определенного интеграла по квадратурным формулам')
    print('f(x) = e^(x)')
    e_x = lambda x: math.e ** x

    print('Введите [A, B]')
    a = float(input('Введите A: '))
    b = float(input('Введите B: '))
    calculate(e_x, a, b)

    polynom_0 = lambda x: 100
    polynom_1 = lambda x: x ** 1
    polynom_2 = lambda x: x ** 2
    polynom_3 = lambda x: x ** 3

    print("------------------------------------------------")
    print("Проверка на полиномах")
    print("------------------------------------------------")
    print('f(x) = 100')
    calculate(polynom_0, a, b)
    print("------------------------------------------------")
    print('f(x) = x')
    calculate(polynom_1, a, b)
    print("------------------------------------------------")
    print('f(x) = x^2')
    calculate(polynom_2, a, b)
    print("------------------------------------------------")
    print('f(x) = x^3')
    calculate(polynom_3, a, b)
