"""
    Autore: Sicky2005
    Corso: Metodi Numerici per l'Ingegneria
    Descrizione: Programma per il Metodo di Interpolazione di Lagrange
"""

def interpolazione_lagrange(x_nodi, y_nodi, x_target):
    """
    x_nodi: Array dei punti x conosciuti (nodi)
    y_nodi: Array dei valori y corrispondenti
    x_target: Il punto x in cui vogliamo trovare il valore y
    """
    # ricavo n dalla lunghezza del vettore
    n = len(x_nodi)

    # inizializzo la somma finale a 0
    risultato_interpolazione = 0.0

    # Ciclo esterno
    for i in range(n):

        # inizializzo il prodotto base Li(x) a 1 (elemento neutro della moltiplicazione)
        polinomio_base_li = 1.0

        # ciclo interno
        for j in range(n):
            if i != j:
                numeratore = x_target - x_nodi[j]
                denominatore = x_nodi[i] - x_nodi[j]
                polinomio_base_li = polinomio_base_li * (numeratore / denominatore)

        # ora moltiplico il polinomio base Li per il valore yi
        termine_corrente = y_nodi[i] * polinomio_base_li

        # aggiungo alla somma totale
        risultato_interpolazione += termine_corrente

        # mostro quanto vale Li e quanto contribuisce questo nodo al totale
        print(f"Nodo {i} (x={x_nodi[i]}): L_{i} = {polinomio_base_li:.5f} -> Contributo (+ y[{i}]*L_{i}) = {termine_corrente:.10f}")

    return risultato_interpolazione


