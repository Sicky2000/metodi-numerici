"""
Modulo per i Metodi di Bracketing (Metodi Chiusi).

Include algoritmi che richiedono un intervallo iniziale [a, b] in cui
la funzione cambia segno, garantendo la convergenza alla radice.

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""

def bisezione(f, a, b, tol=1e-6, max_iter=100):
    """
    Trova la radice di f(x) nell'intervallo [a, b] usando il metodo di Bisezione.

    Args:
        f (callable): La funzione di cui trovare lo zero.
        a (float): Estremo inferiore dell'intervallo.
        b (float): Estremo superiore dell'intervallo.
        tol (float): Tolleranza per l'errore relativo stimato.
        max_iter (int): Numero massimo di iterazioni.

    Returns:
        float: La radice approssimata.

    Raises:
        ValueError: Se f(a) e f(b) hanno lo stesso segno.
        RuntimeError: Se il metodo non converge entro max_iter.
    """

    # Valutiamo la funzione agli estremi
    fa = f(a)
    fb = f(b)

    if fa * fb >= 0:
        raise ValueError("La funzione deve avere segni opposti agli estremi a e b (Bracketing non valido).")

    xr = a  # Valore iniziale (o un valore a caso dentro l'intervallo)
    xr_old = a

    for i in range(max_iter):
        # Stima della radice (punto medio)
        xr = (a + b) / 2
        fxr = f(xr)

        # Calcolo dell'errore relativo (evitiamo divisione per zero)
        if xr != 0 and i > 0:
            ea = abs((xr - xr_old) / xr)
            if ea < tol:
                return xr

        # Se abbiamo trovato lo zero esatto
        if fxr == 0:
            return xr

        # Aggiornamento dell'intervallo (Regola dei segni)
        if fa * fxr < 0:
            b = xr
            fb = fxr  # Aggiorniamo fb (anche se in bisezione pura non serve strettamente, è buona prassi)
        else:
            a = xr
            fa = fxr  # Aggiorniamo fa (fondamentale perché 'a' è cambiato)

        xr_old = xr

    # Se usciamo dal ciclo senza return, significa che non abbiamo raggiunto la tolleranza
    raise RuntimeError(f"Il metodo di bisezione non ha convertito dopo {max_iter} iterazioni.")


def falsa_posizione(f, a, b, tol=1e-6, max_iter=100):
    """
    Trova la radice di f(x) in [a, b] usando il metodo di Falsa Posizione (Variante Illinois).

    La variante Illinois gestisce il caso in cui uno degli estremi rimane stagnante,
    dimezzando il valore della funzione in quell'estremo per accelerare la convergenza.

    Args:
        f (callable): La funzione di cui trovare lo zero.
        a (float): Estremo inferiore dell'intervallo.
        b (float): Estremo superiore dell'intervallo.
        tol (float): Tolleranza per l'errore relativo.
        max_iter (int): Numero massimo di iterazioni.

    Returns:
        float: La radice approssimata.

    Raises:
        ValueError: Se f(a) e f(b) hanno lo stesso segno.
        RuntimeError: Se il metodo non converge entro max_iter.
    """

    fa = f(a)
    fb = f(b)

    if fa * fb >= 0:
        raise ValueError("La funzione deve avere segni opposti agli estremi a e b.")

    xr = a
    xr_old = a

    # Contatori per la stagnazione (Illinois algorithm)
    ia = 0  # Contatore per estremo a
    ib = 0  # Contatore per estremo b

    for i in range(max_iter):
        # Formula della Falsa Posizione (intersezione della secante)
        # xr = b - fb * (a - b) / (fa - fb)
        xr = b - (fb * (a - b)) / (fa - fb)
        fxr = f(xr)

        # Calcolo errore relativo (dalla seconda iterazione)
        if xr != 0 and i > 0:
            ea = abs((xr - xr_old) / xr)
            if ea < tol:
                return xr

        # Check esatto
        if fxr == 0:
            return xr

        # Logica di aggiornamento e Variante Illinois
        test = fa * fxr

        if test < 0:
            # La radice è tra a e xr
            b = xr
            fb = fxr
            ib = 0  # b si è mosso, resetto contatore
            ia += 1  # a è stato fermo
            if ia >= 2:
                fa /= 2  # Penalizzo fa (Illinois)
        else:
            # La radice è tra xr e b
            a = xr
            fa = fxr
            ia = 0  # a si è mosso, resetto contatore
            ib += 1  # b è stato fermo
            if ib >= 2:
                fb /= 2  # Penalizzo fb (Illinois)

        xr_old = xr

    raise RuntimeError(f"Il metodo di Falsa Posizione non ha convertito dopo {max_iter} iterazioni.")