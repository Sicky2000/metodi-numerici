def newton(func, dfunc, x0, tol, max_iter):
    x1 = x0

    for i in range(max_iter):
        fx1 = func(x1)
        dfx1 = dfunc(x1)

        if dfx1 == 0:
            print("La derivata è zero.")
            return None

        x1 = x1 - fx1 / dfx1

        es = abs(x1 - x0)

        if es < tol:
            return x1, i + 1

        if x1 != 0 and (es / abs(x1)) < tol:
            return x1, i + 1

        x0 = x1

    return x1, max_iter