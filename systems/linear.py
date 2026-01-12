"""
Modulo per la risoluzione di Sistemi Lineari (Metodi Diretti).

Include algoritmi per risolvere sistemi Ax = b.

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""

import numpy as np

def gauss_elimination(A, b, tol=1e-6):
    """
    Risolve il sistema lineare Ax = b usando l'eliminazione di Gauss
    con pivoting parziale scalato.

    Args:
        A (list or np.array): Matrice dei coefficienti (n x n).
        b (list or np.array): Vettore dei termini noti (n).
        tol (float): Tolleranza per determinare se la matrice è singolare.

    Returns:
        np.array: Il vettore soluzione x.

    Raises:
        ValueError: Se le dimensioni non coincidono.
        np.linalg.LinAlgError: Se la matrice è singolare (pivot vicino a 0).
    """
    # Copia e conversione in float per evitare modifiche agli input originali
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)

    if A.shape != (n, n):
        raise ValueError("La matrice A deve essere quadrata e compatibile con il vettore b.")

    # Vettore di scaling (massimo valore assoluto per riga)
    # Serve per il "Scaled Partial Pivoting"
    s = np.max(np.abs(A), axis=1)

    # --- Eliminazione in avanti ---
    for k in range(n - 1):
        # 1. Scelta del Pivot
        # Cerchiamo la riga p (tra k e n) che massimizza |A[i,k]| / s[i]
        # Usiamo argmax con slicing per evitare cicli for
        relative_peaks = np.abs(A[k:, k]) / s[k:]
        max_idx = np.argmax(relative_peaks) # Indice relativo alla slice
        p = max_idx + k # Indice assoluto

        # Controllo singolarità
        if abs(A[p, k] / s[p]) < tol:
            raise np.linalg.LinAlgError("Matrice singolare (pivot troppo piccolo).")

        # 2. Scambio righe (se necessario)
        if p != k:
            A[[k, p]] = A[[p, k]]
            b[[k, p]] = b[[p, k]]
            s[[k, p]] = s[[p, k]]

        # 3. Eliminazione
        # Per ogni riga sotto k, eliminiamo l'elemento nella colonna k
        for i in range(k + 1, n):
            factor = A[i, k] / A[k, k]
            # Vettorializzazione: aggiorniamo tutta la riga i in un colpo solo
            # A[i, k:] = A[i, k:] - factor * A[k, k:]
            # Nota: ottimizziamo partendo da k+1 perché la colonna k diventa 0
            A[i, k+1:] -= factor * A[k, k+1:]
            b[i] -= factor * b[k]

    # Controllo finale sull'ultimo elemento
    if abs(A[n - 1, n - 1]) < tol:
         raise np.linalg.LinAlgError("Matrice singolare (ultimo pivot troppo piccolo).")

    # --- Sostituzione all'indietro ---
    x = np.zeros(n)
    x[n - 1] = b[n - 1] / A[n - 1, n - 1]

    for i in range(n - 2, -1, -1):
        # Somma prodotto riga i * soluzioni già trovate
        # sum_ax = np.sum(A[i, i+1:] * x[i+1:])
        # Oppure prodotto scalare:
        sum_ax = np.dot(A[i, i+1:], x[i+1:])
        x[i] = (b[i] - sum_ax) / A[i, i]

    return x


def thomas(e, f, g, b):
    """
    Risolve un sistema tridiagonale Ax = b usando l'algoritmo di Thomas (TDMA).

    La matrice A è definita da tre vettori:
    - e: Diagonale inferiore (sub-diagonal). e[0] è ignorato.
    - f: Diagonale principale (main diagonal).
    - g: Diagonale superiore (super-diagonal). g[n-1] è ignorato.

    Args:
        e (list/array): Vettore diagonale inferiore (lunghezza n).
        f (list/array): Vettore diagonale principale (lunghezza n).
        g (list/array): Vettore diagonale superiore (lunghezza n).
        b (list/array): Vettore dei termini noti (lunghezza n).

    Returns:
        np.array: Il vettore soluzione x.

    Raises:
        ValueError: Se le dimensioni dei vettori non coincidono.
    """
    n = len(f)
    if len(e) != n or len(g) != n or len(b) != n:
        raise ValueError("I vettori e, f, g, b devono avere tutti la stessa lunghezza n.")

    # Copiamo per non modificare gli originali
    # Nota: usiamo float per evitare divisioni intere
    e_work = np.array(e, dtype=float)
    f_work = np.array(f, dtype=float)
    b_work = np.array(b, dtype=float)
    g_work = np.array(g, dtype=float)

    # --- Decomposizione (Eliminazione in avanti) ---
    for k in range(1, n):
        if f_work[k - 1] == 0:
            # Se capita spesso, servirebbe pivoting (ma Thomas classico non lo fa)
            raise ValueError(f"Pivot nullo in k={k - 1}. Thomas algorithm fallisce.")

        factor = e_work[k] / f_work[k - 1]
        # Aggiorniamo diagonale principale
        f_work[k] = f_work[k] - factor * g_work[k - 1]
        # Aggiorniamo termine noto
        b_work[k] = b_work[k] - factor * b_work[k - 1]

    # --- Sostituzione all'indietro ---
    x = np.zeros(n)

    # Calcolo ultima x
    if f_work[n - 1] == 0:
        raise ValueError("Pivot nullo nell'ultimo elemento. Sistema singolare.")

    x[n - 1] = b_work[n - 1] / f_work[n - 1]

    # Ciclo inverso
    for k in range(n - 2, -1, -1):
        x[k] = (b_work[k] - g_work[k] * x[k + 1]) / f_work[k]

    return x