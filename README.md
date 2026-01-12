# Metodi Numerici per l'Ingegneria

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Una libreria Python completa, modulare ed efficiente che implementa i principali algoritmi numerici per la risoluzione di equazioni, sistemi lineari e non lineari, interpolazione e integrazione.

Sviluppato come progetto per il corso di **Metodi Numerici per l'Ingegneria**.

## 🚀 Caratteristiche

Il progetto è stato rifattorizzato in pacchetti logici per facilitare l'utilizzo e la manutenzione. Fa ampio uso di **NumPy** per garantire efficienza e vettorializzazione dei calcoli, sostituendo i classici cicli lenti di Python dove possibile.

### 📦 Moduli Disponibili

#### 1. `roots` (Ricerca di Radici)
Metodi per trovare gli zeri di funzioni non lineari $f(x) = 0$.
- **Metodi Chiusi (Bracketing):**
  - Bisezione (`bisezione`)
  - Falsa Posizione (`falsa_posizione`)
- **Metodi Aperti:**
  - Newton-Raphson (`newton_raphson`)
  - Secanti (`secanti`)
  - Punto Fisso (`fixed_point`)

#### 2. `systems` (Sistemi di Equazioni)
Risolutori per sistemi lineari $Ax = b$ e non lineari.
- **Metodi Diretti:**
  - Eliminazione di Gauss con Pivoting Parziale Scalato (`gauss_elimination`)
  - Algoritmo di Thomas per matrici tridiagonali (`thomas`)
- **Metodi Iterativi:**
  - Jacobi (`jacobi`)
  - Gauss-Seidel con Rilassamento SOR (`gauss_seidel`)
- **Sistemi Non Lineari:**
  - Metodo di Broyden (`broyden`) - Metodo Quasi-Newton

#### 3. `interpolation` (Interpolazione)
Costruzione di polinomi interpolanti.
- Polinomio di Lagrange (`lagrange`)
- Metodo di Newton alle Differenze Divise (`newton`)
- Generatore di nodi di **Chebyshev** (`chebyshev_nodes`) per minimizzare il fenomeno di Runge.

#### 4. `integration` (Integrazione Numerica)
Calcolo di integrali definiti.
- Regola dei Trapezi composta (`trapezoidal`)
- **Regola di Simpson Mista** (`simpson`): Algoritmo intelligente che combina Simpson 1/3 (per intervalli pari) e Simpson 3/8 (per gestire intervalli dispari) mantenendo un ordine di accuratezza $O(h^4)$.

---

## 🛠️ Installazione e Requisiti

Assicurati di avere Python installato. Le dipendenze principali sono **NumPy** e **Matplotlib**.

~~~bash
# Clona la repository
git clone https://github.com/tuo-username/metodi-numerici.git
cd metodi-numerici

# Installa le dipendenze
pip install -r requirements.txt
~~~

---

## 💻 Esempio di Utilizzo

Ecco come utilizzare i moduli della libreria nei tuoi script. Grazie alla struttura a pacchetti, l'importazione è semplice e intuitiva.

### 1. Ricerca di Radici (Newton)
~~~python
from roots import newton_raphson

def f(x):
    return x**3 - x - 1

def df(x):
    return 3*x**2 - 1

# Trova la radice partendo da x0 = 1.5
radice = newton_raphson(f, df, x0=1.5)
print(f"Radice trovata: {radice:.6f}")
~~~

### 2. Integrazione Numerica (Simpson Misto)
~~~python
import numpy as np
from integration import simpson

def g(x):
    return np.sin(x)

# Calcola l'integrale di sin(x) da 0 a pi greco
# Nota: Simpson gestisce automaticamente intervalli pari e dispari
area = simpson(g, a=0, b=np.pi, n=11)
print(f"Area calcolata: {area:.6f}")
~~~

### 3. Interpolazione e Nodi di Chebyshev
~~~python
from interpolation import newton, chebyshev_nodes

# Genera nodi ottimali per evitare oscillazioni
x_nodes = chebyshev_nodes(a=-5, b=5, n=10)
y_nodes = [val**2 for val in x_nodes] # Esempio y = x^2

target = 2.5
stima = newton(x_nodes, y_nodes, target)
print(f"Valore interpolato in x={target}: {stima:.4f}")
~~~

---

## 📂 Struttura del Progetto

Il repository è organizzato come segue:

~~~text
metodi-numerici/
├── integration/          # Metodi di integrazione (Trapezi, Simpson)
│   ├── __init__.py
│   ├── simpson.py
│   └── trapezoidal.py
├── interpolation/        # Metodi di interpolazione
│   ├── __init__.py
│   └── polynomial.py
├── roots/                # Ricerca zeri
│   ├── __init__.py
│   ├── bracketing.py
│   └── open_methods.py
├── systems/              # Sistemi lineari e non lineari
│   ├── __init__.py
│   ├── iterative.py
│   ├── linear.py
│   └── nonlinear.py
├── requirements.txt      # Dipendenze del progetto
├── LICENSE               # Licenza GPL-3.0
└── README.md             # Documentazione
~~~

---

## ✍️ Autore

**Sicky2005**:
Studente di Ingegneria della Trasformazione Digitale.

## 📄 Licenza

Questo progetto è distribuito sotto licenza **GPL-3.0**.
Vedi il file [LICENSE](LICENSE) per maggiori dettagli.