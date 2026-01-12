"""
Modulo per i Metodi Aperti.

Include algoritmi che richiedono uno o più punti iniziali e non necessitano
di un intervallo di bracketing. Possono divergere ma sono generalmente più veloci.

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""


def fixed_point(g, x0, tol=1e-6, max_iter=100):
    """
    Trova la radice usando il Metodo del Punto Fisso (Iterazione funzionale).
    Risolve x = g(x).

    Attenzione: La funzione g(x) deve essere derivata da f(x)=0 tale che x = g(x).
    Il metodo converge solo se |g'(x)| < 1 nell'intorno della soluzione.

    Args:
        g (callable): La funzione di iterazione g(x).
        x0 (float): Stima iniziale.
        tol (float): Tolleranza per l'errore relativo.
        max_iter (int): Numero massimo di iterazioni.

    Returns:
        float: La radice approssimata.

    Raises:
        RuntimeError: Se il metodo diverge o non converge entro max_iter.
    """
    xr = x0
    xr_old = x0

    for i in range(max_iter):
        xr = g(xr_old)

        # Calcolo errore relativo
        if xr != 0:
            ea = abs((xr - xr_old) / xr)
            if ea < tol:
                return xr

        # Caso raro: convergenza esatta a 0
        if xr == xr_old:
            return xr

        xr_old = xr

    raise RuntimeError(
        f"Il metodo del Punto Fisso non ha convertito dopo {max_iter} iterazioni (possibile divergenza).")


def newton_raphson(f, df, x0, tol=1e-6, max_iter=100):
    """
    Trova la radice di f(x) usando il metodo di Newton-Raphson (Metodo delle tangenti).
    Richiede la derivata prima analitica f'(x).

    Convergenza quadratica (molto veloce) se x0 è vicino alla soluzione.

    Args:
        f (callable): La funzione f(x).
        df (callable): La derivata prima f'(x).
        x0 (float): Stima iniziale.
        tol (float): Tolleranza per l'errore relativo.
        max_iter (int): Numero massimo di iterazioni.

    Returns:
        float: La radice approssimata.

    Raises:
        ValueError: Se la derivata si annulla (divisione per zero).
        RuntimeError: Se il metodo non converge entro max_iter.
    """
    xr = x0

    for i in range(max_iter):
        fx = f(xr)
        dfx = df(xr)

        # Controllo derivata nulla (tangente orizzontale)
        if dfx == 0:
            raise ValueError(f"Derivata nulla in x={xr}. Il metodo fallisce (divisione per zero).")

        # Passo di Newton: x_new = x_old - f(x)/f'(x)
        xr_new = xr - (fx / dfx)

        # Calcolo errore relativo
        if xr_new != 0:
            ea = abs((xr_new - xr) / xr_new)
            if ea < tol:
                return xr_new

        # Caso convergenza esatta (f(x) = 0)
        if fx == 0:  # O un check tipo abs(fx) < epsilon
            return xr

        # Aggiornamento per la prossima iterazione
        xr = xr_new

    raise RuntimeError(f"Il metodo di Newton non ha convertito dopo {max_iter} iterazioni.")


def secanti(f, x0, x1, tol=1e-6, max_iter=100):
    """
    Trova la radice di f(x) usando il metodo delle Secanti.

    È simile al metodo di Newton, ma approssima la derivata usando
    una differenza finita tra due iterazioni precedenti (x0 e x1).
    Utile quando la derivata analitica f'(x) è difficile da calcolare.

    Args:
        f (callable): La funzione f(x).
        x0 (float): Prima stima iniziale.
        x1 (float): Seconda stima iniziale.
        tol (float): Tolleranza per l'errore relativo.
        max_iter (int): Numero massimo di iterazioni.

    Returns:
        float: La radice approssimata.

    Raises:
        ValueError: Se la differenza f(x1)-f(x0) è troppo piccola (divisione per zero).
        RuntimeError: Se il metodo non converge entro max_iter.
    """

    f0 = f(x0)
    f1 = f(x1)

    xr = x1  # Inizializzazione

    for i in range(max_iter):

        denom = f1 - f0

        # Controllo sicurezza numerica (Secante orizzontale)
        if abs(denom) < 1e-12:
            raise ValueError("Il denominatore è nullo o troppo piccolo (Secante orizzontale).")

        # Formula delle Secanti
        xr = x1 - (f1 * (x1 - x0) / denom)

        # Calcolo errore relativo
        if xr != 0:
            ea = abs((xr - x1) / xr)
            if ea < tol:
                return xr
        elif abs(xr - x1) < tol:  # Fallback su errore assoluto se xr ~= 0
            return xr

        # Shift delle variabili per la prossima iterazione
        x0 = x1
        f0 = f1

        x1 = xr
        f1 = f(xr)

    raise RuntimeError(f"Il metodo delle Secanti non ha convertito dopo {max_iter} iterazioni.")