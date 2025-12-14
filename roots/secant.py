def secant(func, x0, x1, tol, max_iter):
    fx1 = func(x1)
    fx0 = func(x0)

    for i in range(max_iter):
        if abs(fx1 - fx0) < 1e-12:
            print("Errore: Il denominatore è zero")
            return None

        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)

        if abs(x2 - x1) < tol:
            print(f"Convergenza dopo {i+1} iterazioni")
            return x2, i+1

        x0 = x1
        fx0 = fx1

        x1 = x2
        fx1 = func(x2)

    print("Numero massimo di iterazioni raggiunto senza convergenza")
    return x1, max_iter