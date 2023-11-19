from rootApproximator import separateRoots, bisectionMethod


def get_table(a, b, m, func):
    step = (b - a) / m
    table = {}
    for j in range(0, m + 1):
        z = a + j * step
        table[z] = func(z)
    return table


def sort_table(table, x, n):
    table_sorted = {k: v for k, v in sorted(table.items(), key=lambda item: abs(item[0] - x))}
    table_sorted_keys = list(table_sorted.keys())
    table_new = {}
    i = 0
    while i <= n:
        table_new[table_sorted_keys[i]] = table_sorted[table_sorted_keys[i]]
        i += 1
    return table_new


def print_table(table, mode=0):
    if mode == 0:
        print("    x_k   |   f(x_k)")
    else:
        print("  f(x_k)  |  f^(-1)(f(x_k))")
    for key in table.keys():
        if key >= 0:
            print(" " + "{:.6f}".format(key), end=" | ")
        else:
            print("{:.6f}".format(key), end=" | ")
        if table[key] >= 0:
            print(" " + "{:.6f}".format(table[key]))
        else:
            print("{:.6f}".format(table[key]))


def change_table(table):
    new_table = {}
    for key in table.keys():
        new_table[table[key]] = key
    return new_table


# table: dict {x0 : f(x0)}
def LagrangeInterpolation(table, x, n):
    args = list(table.keys())
    res = 0
    for k in range(0, n + 1):
        x_k = args[k]
        x_k_value = table[x_k]
        local_res = 1
        for i in range(0, n + 1):
            if k != i:
                local_res *= (x - args[i]) / (args[k] - args[i])
        res += x_k_value * local_res
    return res


def LagrangeInterpolationPolynom(table, n):
    args = list(table.keys())
    res_dict = {}
    for k in range(0, n + 1):
        local_data = []
        x_k = args[k]
        x_k_value = table[x_k]
        for i in range(0, n + 1):
            if k != i:
                local_data.append((args[i], args[k] - args[i]))
        res_dict[x_k_value] = local_data
    return res_dict


def local_data_mult(x, x_k_value, local_data):
    local_res = 1
    for data in local_data:
        local_res *= (x - data[0]) / data[1]
    return x_k_value * local_res


def polynom(res_dict, F):
    return lambda x: sum([local_data_mult(x, x_k_value, res_dict[x_k_value]) for x_k_value in res_dict.keys()]) - F


# get divided differences
def GetDiv(table, arguments):
    if len(arguments) == 2:
        return (table[arguments[1]] - table[arguments[0]]) / (arguments[1] - arguments[0])
    else:
        numOfArguments = len(arguments)
        return (GetDiv(table, arguments[1: numOfArguments]) - GetDiv(table, arguments[0: numOfArguments - 1])) / (
                arguments[numOfArguments - 1] - arguments[0])


# table: dict {x0 : f(x0)}
def NewtonInterpolation(table, x, n):
    arguments = list(table.keys())
    result = 0
    rightMultiplicant = 1

    pairIndex = 0
    result += table[arguments[0]]
    for arguementValuePair in table:
        if pairIndex == 0:
            pairIndex += 1
            continue

        rightMultiplicant *= (x - arguments[pairIndex - 1])
        leftMultiplicant = GetDiv(table, arguments[0: pairIndex + 1])
        result += leftMultiplicant * rightMultiplicant
        pairIndex += 1
    return result


def monotonic_function(func):
    args_count = int(input("Введите число значений в таблице (m+1): "))
    m = args_count - 1
    a = float(input("Введите A >= 0 (начало отрезка): "))
    b = float(input("Введите B > 0 (конец отрезка): "))
    while a < 0 or b <= 0 or a >= b:
        if 0 >= b > a:
            break
        print("Введенный промежуток не удовлетворяет условию строгой монотонности функции, введите другой")
        a = float(input("Введите A >= 0 (начало отрезка): "))
        b = float(input("Введите B > 0 (конец отрезка): "))

    table = get_table(a, b, m, func)
    print()
    print_table(table)
    print()

    new_table = change_table(table)
    print_table(new_table, 1)
    print()

    exit_program = False

    while not exit_program:
        F = float(input("Введите F (значение f(x)): "))
        n = int(input(f"Введите N (степень интерполяционного многочлена, не превышающая {m}): "))
        while n > m:
            n = int(input("Получено недопустимое значение n, введите другое: "))

        new_sorted_table = sort_table(new_table, F, n)
        print()
        print_table(new_sorted_table)
        print()

        X = LagrangeInterpolation(new_sorted_table, F, n)
        print(f"X = {X}")
        print(f"|f(X) - F|: {abs(func(X) - F)}")
        print()

        exit_program = input('Перейти ко второму способу решения? (y/n): ').lower().strip() == 'y'
        print()
    return


def non_monotonic_function(func):
    args_count = int(input("Введите число значений в таблице (m+1): "))
    m = args_count - 1
    a = float(input("Введите A (начало отрезка): "))
    b = float(input("Введите B (конец отрезка): "))

    table = get_table(a, b, m, func)
    print()
    print_table(table)
    print()

    exit_program = False

    while not exit_program:
        F = float(input("Введите F (значение f(x)): "))
        n = int(input(f"Введите N (степень интерполяционного многочлена, не превышающая {m}): "))
        while n > m:
            n = int(input("Получено недопустимое значение n, введите другое: "))
        h = int(input("Введите H (кол-во интервалов разделения корней): "))
        eps = float(input("Введите эпсилон (погрешность): "))

        lagrange = LagrangeInterpolationPolynom(table, n)
        poly = polynom(lagrange, F)

        foundIntervals, foundRoots = separateRoots(poly, a, b, h)
        roots = []
        for interval in foundIntervals:
            roots.append(bisectionMethod(poly, interval[0], interval[1], eps, info=False))

        print(f"\nКорни P_n(x) = F, найденные методом бисекций: {roots}")
        for i in range(len(roots)):
            root = roots[i]
            print(f"|f({root}) - F|: {abs(func(root) - F)}")

        exit_program = input('Завершить программу? (y/n): ').lower().strip() == 'y'
        print()


if __name__ == "__main__":
    print("Задача обратного интерполирования")
    print("Вариант №6")
    print("f(x): x^2 / (1 + x^2)\na=0, b=1, m=10, ε=10^(-12)\n")
    func = lambda y: y ** 2 / (1 + y ** 2)

    print("1) Первый способ решения. Функция строго монотонна (x >= 0 или x <= 0)\n")
    monotonic_function(func)
    print("2) Второй способ решения. Произвольная функция\n")
    non_monotonic_function(func)
