"""
Modulo per l'Integrazione con il Metodo dei Trapezi.

Questo modulo implementa la regola del trapezio per il calcolo approssimato
di integrali definiti, gestendo sia l'applicazione su segmento singolo
che la formula composta per intervalli multipli.

Autore:      Sicky2005
Corso:       Metodi Numerici per l'Ingegneria
Descrizione: Programma per il Metodo dei Trapezi
"""

def trap_rule(func, a, b, n):
    """
    Calcola l'integrale definito usando la regola del Trapezio composta.

    Args:
        func: La funzione da integrare (deve accettare x e ritornare y)
        a (float): Inizio dell'intervallo
        b (float): Fine dell'intervallo
        n (int): Numero di intervalli

    Returns:
        float: Valore approssimato dell'integrale
    """
    if n < 1:
        raise ValueError("Il numero di intervalli n deve essere >= 1")

    # Calcolo del passo (h)
    h = (b - a) / n

    # Valutazione estremi (f(a) + f(b))
    somma = func(a) + func(b)

    # Sommatoria dei punti interni moltiplicati per 2
    for i in range(1, n):
        x = a + i * h
        somma += 2 * func(x)

    # Calcolo finale
    return (h / 2) * somma