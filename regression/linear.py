"""
Modulo per la Regressione Lineare.

Questo modulo implementa il metodo dei minimi quadrati per il calcolo
della regressione lineare semplice, includendo il calcolo degli errori standard
e del coefficiente di determinazione R^2.

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
Descrizione: Programma per la Regressione (Libreria matematica)
"""

def linear_regression(x, y):
    """
    Esegue la regressione lineare col metodo dei minimi quadrati.
    Restituisce i coefficienti della retta y = a1*x + a0 e le statistiche di errore.

    Args:
        x (list): Lista dei valori della variabile indipendente
        y (list): Lista dei valori della variabile dipendente

    Returns:
        tuple: (a1, a0, syx, r2)
            - a1: Pendenza (slope)
            - a0: Intercetta (intercept)
            - syx: Errore standard della stima
            - r2: Coefficiente di determinazione (0 <= r2 <= 1)
        """
    n = len(x)

    # Controllo di sicurezza: servono almeno 3 punti per calcolare syx (n-2)
    if n < 3:
        raise ValueError("Sono necessari almeno 3 punti dati.")

    # Inizializzazione Variabili
    sumx = 0 # Somma di tutte le x (Σx). Serve per calcolare la media xm
    sumy = 0# Somma di tutte le y (Σy). Serve per calcolare la media ym
    sumxy = 0 # Somma dei prodotti x*y (Σxy). Fondamentale per la covarianza
    sumx2 = 0 # Somma dei quadrati di x (Σx²). Fondamentale per la varianza di x
    st = 0 # (Sum Total) Dispersione totale dei dati rispetto alla MEDIA.
    sr = 0 # (Sum Residuals) Dispersione dei dati rispetto alla RETTA (errore).

    # Primo Ciclo (Accumulo Somme)
    for i in range(n):
        sumx += x[i]
        sumy += y[i]
        sumxy += (x[i] * y[i])
        sumx2 += (x[i] ** 2)

    # Calcolo Medie e Coefficienti
    xm = sumx / n
    ym = sumy / n

    # Calcolo a1 (Pendenza)
    numerator = (n * sumxy) - (sumx * sumy)
    denominator = (n * sumx2) - (sumx * sumx)

    if denominator == 0:
        raise ValueError("Impossibile calcolare: il denominatore è 0 (tutti gli x sono uguali?)")

    a1 = numerator / denominator

    # Calcolo a0 (Intercetta)
    a0 = ym - (a1 * xm)

    # Secondo Ciclo (Calcolo Errori)
    for i in range(n):
        # st: scarto quadratico totale rispetto alla media
        st = st + (y[i] - ym) ** 2
        # sr: scarto quadratico dei residui (errore della regressione)
        sr = sr + (y[i] - a1 * x[i] - a0) ** 2

    # Calcolo Statistiche Finali
    # syx: Errore standard della stima
    syx = (sr / (n - 2)) ** 0.5

    # r2: Coefficiente di determinazione
    # Se st è 0 (tutti i valori y sono uguali), r2 non è definito matematicamente, gestiamo il caso.
    if st == 0:
        r2 = 1.0  # Adattamento perfetto (linea orizzontale)
    else:
        r2 = (st - sr) / st

    return a1, a0, syx, r2