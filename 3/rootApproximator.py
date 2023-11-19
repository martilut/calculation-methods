from math import cos, sin


def separateRoots(func, a, b, n):
    step_counter = 0
    step = (b - a) / n
    foundIntervals = list()
    foundRoots = list()

    start = a
    end = start + step

    while end <= b:
        start_value = func(start)
        end_value = func(end)

        if start_value == 0:
            foundRoots.append(start_value)
        elif start_value * end_value < 0:
            step_counter += 1
            foundIntervals.append((start, end))
        else:
            pass

        start = end
        end = end + step

    return foundIntervals, foundRoots


def printApproximationResults(name, initial_approx, step_counter, approx, last_segment_length, residual):
    print(f'Method name: {name}\n'
          f'Initial approximation: {initial_approx}\n'
          f'Step count: {step_counter}\n'
          f'Approximation by method: {approx}\n'
          f'Last difference: {last_segment_length}\n'
          f'|f(x)-0|: {residual}')


def getStartingValue(a, b, func, second_derivative):
    step = (b - a) / 1000
    x = a

    while True:
        if func(x) * second_derivative(x) > 0:
            break
        x += step
        if x > b:
            raise Exception('No convergence')

    return x


def bisectionMethod(func, a, b, eps, info=True):
    step_counter = 0
    local_a = a
    local_b = b

    while (local_b - local_a) > 2 * eps:
        step_counter += 1
        c = (local_a + local_b) / 2
        if func(c) * func(local_a) <= 0:
            local_b = c
        else:
            local_a = c

    approx = (local_a + local_b) / 2
    initial_approx = (a + b) / 2
    residual = abs(func(approx))

    if info:
        printApproximationResults("Bisection method", initial_approx, step_counter, approx, local_b - local_a, residual)
    return approx


def newtonMethod(func, derivative, second_derivative, a, b, eps):
    step_counter = 0
    initial_approx = x_prev = getStartingValue(a, b, func, second_derivative)

    while True:
        step_counter += 1
        x_curr = x_prev - func(x_prev) / derivative(x_prev)
        if abs(x_curr - x_prev) <= eps:
            break
        x_prev = x_curr

    last_segment_length = x_prev - x_curr
    approx = x_curr
    residual = abs(func(approx))

    printApproximationResults('Newton method', initial_approx, step_counter, approx, last_segment_length, residual)


def modifiedNewtonMethod(func, derivative, second_derivative, a, b, eps):
    step_counter = 0
    initial_approx = x_prev = x0 = getStartingValue(a, b, func, second_derivative)

    while True:
        step_counter += 1
        x_curr = x_prev - func(x_prev) / derivative(x0)
        if abs(x_curr - x_prev) <= eps:
            break
        x_prev = x_curr

    last_segment_length = x_prev - x_curr
    approx = x_curr
    residual = abs(func(approx))

    printApproximationResults('Modified Newton method', initial_approx, step_counter, approx, last_segment_length,
                              residual)


def secantMethod(func, a, b, eps):
    step_counter = 0
    x_prev = a
    x_curr = b

    while True:
        step_counter += 1
        x_next = x_curr - func(x_curr) * (x_curr - x_prev) / (func(x_curr) - func(x_prev))

        if abs(x_next - x_curr) <= eps:
            break

        x_prev = x_curr
        x_curr = x_next

    last_segment_length = x_next - x_curr
    approx = x_next
    initial_approx = (a + b) / 2
    residual = abs(func(approx))

    printApproximationResults('Secant method', initial_approx, step_counter, approx, last_segment_length, residual)


if __name__ == "__main__":
    print('Численные методы решения нелинейных уравнений')
    print('f(x) = 8cos(x) - x - 6')
    print('[A ; B] = [-9 ; 1]')
    print('eps = 10^(-7)\n')

    function = lambda x: 8 * cos(x) - x - 6
    first_deriv = lambda x: -8 * sin(x) - 1
    second_deriv = lambda x: -8 * cos(x)

    a = float(input("Input A: "))
    b = float(input("Input B: "))
    n = int(input("Input N: "))
    eps = float(input("Input epsilon: "))

    print('\n---Root separation---')
    foundIntervals, foundRoots = separateRoots(function, a, b, n)
    print(f"\nRoots found during root separation: {len(foundRoots)}\n"
          f"List of roots: {foundRoots}\n")
    print(f"Sign changing intervals found during root separation: {len(foundIntervals)}\n"
          f"List of intervals:")
    for interval in foundIntervals:
        print(f" [{interval[0]}, {interval[1]}]")
    print('\n---Clarification of roots---\n')

    for interval in foundIntervals:
        print(f"Current interval: [{interval[0]}, {interval[1]}]\n")
        bisectionMethod(function, interval[0], interval[1], eps)
        print('\n')
        newtonMethod(function, first_deriv, second_deriv, interval[0], interval[1], eps)
        print('\n')
        modifiedNewtonMethod(function, first_deriv, second_deriv, interval[0], interval[1], eps)
        print('\n')
        secantMethod(function, interval[0], interval[1], eps)
        print('\n')