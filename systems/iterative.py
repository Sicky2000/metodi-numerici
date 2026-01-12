"""
Modulo per la risoluzione di Sistemi Lineari (Metodi Iterativi).

Ideale per sistemi di grandi dimensioni o matrici sparse dove
l'eliminazione di Gauss sarebbe troppo costosa.

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
"""

import numpy as np


def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=100, omega=1.0):
    """
    Risolve il sistema Ax = b usando il metodo di Gauss-Seidel.
    Supporta il rilassamento (SOR - Successive Over-Relaxation) tramite il parametro omega.

    Args:
        A (np.array): Matrice dei coefficienti (n x n).
        b (np.array): Vettore dei termini noti (n).
        x0 (np.array, optional): Stima iniziale. Se None, usa vettore nullo.
        tol (float): Tolleranza per l'errore relativo (norma euclidea).
        max_iter (int): Numero massimo di iterazioni.
        omega (float): Fattore di rilassamento (1.0 = Gauss-Seidel standard).
                       0 < omega < 1: Sotto-rilassamento (per convergenza difficile)
                       1 < omega < 2: Sovra-rilassamento (per accelerare)

    Returns:
        np.array: Il vettore soluzione x.

    Raises:
        ValueError: Se la matrice ha elementi diagonali nulli.
        RuntimeError: Se il metodo non converge.
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)

    # Controllo diagonale dominante o zeri sulla diagonale
    diag = np.diag(A)
    if np.any(diag == 0):
        raise ValueError("Elemento diagonale nullo. Impossibile applicare Gauss-Seidel.")

    if x0 is None:
        x = np.zeros(n)
    else:
        x = np.array(x0, dtype=float)

    for k in range(max_iter):
        x_old = x.copy()

        # Iterazione sulle righe
        for i in range(n):
            # Calcolo sigma: somma di A[i,j] * x[j] per tutti i j != i
            # Nota: x contiene già i valori aggiornati per j < i (caratteristica di GS)
            # Ottimizzazione: prodotto scalare intera riga - elemento diagonale
            sigma = np.dot(A[i, :], x) - A[i, i] * x[i]

            # Calcolo nuovo valore (Formula di Gauss-Seidel)
            x_new = (b[i] - sigma) / A[i, i]

            # Applicazione del rilassamento (SOR)
            x[i] = omega * x_new + (1 - omega) * x_old[i]

        # Controllo convergenza (Norma dell'errore relativo)
        # Evitiamo divisione per zero se x è nullo
        norm_x = np.linalg.norm(x)
        if norm_x == 0:
            diff = np.linalg.norm(x - x_old)
        else:
            diff = np.linalg.norm(x - x_old) / norm_x

        if diff < tol:
            return x

    raise RuntimeError(f"Gauss-Seidel non ha convertito dopo {max_iter} iterazioni.")


def jacobi(A, b, x0=None, tol=1e-6, max_iter=100):
    """
    Risolve il sistema Ax = b usando il metodo di Jacobi.
    A differenza di Gauss-Seidel, aggiorna tutte le componenti simultaneamente.

    Args:
        A (np.array): Matrice dei coefficienti.
        b (np.array): Vettore dei termini noti.
        x0 (np.array, optional): Stima iniziale.
        tol (float): Tolleranza.
        max_iter (int): Max iterazioni.

    Returns:
        np.array: Soluzione x.
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)

    diag = np.diag(A)
    if np.any(diag == 0):
        raise ValueError("Elemento diagonale nullo.")

    # Matrice D (Diagonale) e R (Resto: L + U)
    # x_new = D^-1 * (b - R * x_old)
    R = A - np.diag(diag)  # Matrice A con diagonale azzerata

    if x0 is None:
        x = np.zeros(n)
    else:
        x = np.array(x0, dtype=float)

    for k in range(max_iter):
        # Calcolo vettoriale (molto veloce in NumPy)
        # b - (A senza diagonale) * x
        numerator = b - np.dot(R, x)
        x_new = numerator / diag

        # Errore relativo
        norm_x = np.linalg.norm(x_new)
        if norm_x == 0:
            diff = np.linalg.norm(x_new - x)
        else:
            diff = np.linalg.norm(x_new - x) / norm_x

        if diff < tol:
            return x_new

        x = x_new

    raise RuntimeError(f"Jacobi non ha convertito dopo {max_iter} iterazioni.")