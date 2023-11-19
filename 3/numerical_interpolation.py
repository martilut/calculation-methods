import math
from prettytable import PrettyTable
from reverse_interpolation import print_table


def get_table(a, m, h, func):
    table = {}
    for j in range(0, m + 1):
        z = a + j * h
        table[z] = func(z)
    return table


def first_derivatives(y, m, h):
    derivatives = []
    for i in range(0, m + 1):
        if i < 2:
            derivatives.append((-3 * y[i] + 4 * y[i + 1] - y[i + 2]) / (2 * h))
            continue
        derivatives.append((3 * y[i] - 4 * y[i - 1] + y[i - 2]) / (2 * h))
    return derivatives


def second_derivatives(y, m, h):
    derivatives = [None]
    for i in range(1, m):
        derivatives.append((y[i + 1] - 2 * y[i] + y[i - 1]) / (h ** 2))
    derivatives.append(None)
    return derivatives


def build_table(pretty, x, m, func, first_der_func, second_der_func, first_derivatives, second_derivatives):
    for i in range(0, m + 1):
        first_abs = abs(first_derivatives[i] - first_der_func(x[i]))
        first_rel_abs = first_abs / abs(first_der_func(x[i]))
        second_abs = abs(second_derivatives[i] - second_der_func(x[i])) if second_derivatives[i] is not None else None
        second_rel_abs = second_abs / abs(second_der_func(x[i])) if second_abs is not None else None
        pretty.add_row(
            [
                {x[i]}, {func(x[i])},
                {first_derivatives[i]}, {first_abs}, {first_rel_abs},
                {second_derivatives[i]}, {second_abs}, {second_rel_abs}
            ]
        )


if __name__ == "__main__":
    print("Нахождение производных таблично-заданной функции по формулам численного дифференцирования")
    print("Вариант №6")
    print("f(x): e^(3x)\n")

    func = lambda y: math.e ** (3 * y)
    first_derivative = lambda y: 3 * (math.e ** (3 * y))
    second_derivative = lambda y: 9 * (math.e ** (3 * y))

    exit_program = False

    while not exit_program:
        args_count = int(input("Введите число значений в таблице (m+1): "))
        m = args_count - 1
        a = float(input("Введите a: "))
        h = float(input("Введите h > 0: "))

        table = get_table(a, m, h, func)
        print()
        print_table(table)
        print()
        first = first_derivatives(list(table.values()), m, h)
        second = second_derivatives(list(table.values()), m, h)

        pretty = PrettyTable()
        pretty.field_names = ["x", "f(x)",
                              "f'(x)ЧД", "Абс. погрешность f'(x)ЧД", "Отн. погрешность f'(x)ЧД",
                              "f''(x)ЧД", "Абс. погрешность f''(x)ЧД", "Отн. погрешность f''(x)ЧД"]
        build_table(pretty, list(table.keys()), m, func, first_derivative, second_derivative, first, second)
        print(pretty)
        print()

        exit_program = input('Завершить программу? (y/n): ').lower().strip() == 'y'
        print()
