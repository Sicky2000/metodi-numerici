def false_pos(func, xl, xu, tol, max_iter):

    iter_count = 0
    xr = 0.0
    ea = 100.0
    fl = func(xl)
    fu = func(xu)
    iu = 0
    il = 0

    while iter_count < max_iter:
        xrold = xr
        xr = xu - fu * (xl - xu) / (fl - fu)
        fr = func(xr)
        iter_count += 1

        if xr != 0 and iter_count > 1:
            ea = abs((xr - xrold) / xr) * 100

        test_val = fl * fr

        if test_val < 0:
            xu = xr
            fu = func(xu)
            iu = 0
            il += 1

            if il >= 2:
                fl = fl / 2
        elif test_val > 0:
            xl = xr
            fl = func(xl)
            il = 0
            iu += 1
            if iu >= 2:
                fu = fu / 2
        else:
            ea = 0

        if ea < tol or iter_count >= max_iter:
            break

    return xr, iter_count