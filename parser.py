# Librerias Esenciales
import sys  # Manejo del sistema
from collections import defaultdict  # Diccionarios avanzados
from tabulate import tabulate  # Tablas elegantes

# Funciones Principales

def leer_gramatica():
    """Obtiene la gramatica desde la entrada estandar"""
    n_reglas = int(input())
    
    reglas = defaultdict(list)
    no_terminales = []
    
    for _ in range(n_reglas):
        linea = input().strip()
        if not linea:
            continue
            
        lado_izq, lado_der = linea.split("->")
        lado_izq = lado_izq.strip()
        
        if lado_izq not in no_terminales:
            no_terminales.append(lado_izq)
            
        alternativas = lado_der.strip().split()
        
        for alt in alternativas:
            for prod in alt.split('|'):
                simbolos = [c if c != 'e' else 'ε' for c in prod.strip()]
                if not simbolos:
                    simbolos = ['ε']
                reglas[lado_izq].append(simbolos)
                
    return reglas, no_terminales[0]

def obtener_terminales(reglas):
    """Extrae los terminales de las producciones"""
    terminales = set()
    
    for producciones in reglas.values():
        for prod in producciones:
            for simb in prod:
                if simb != 'ε' and (not simb.isupper()):
                    terminales.add(simb)
                    
    terminales.discard('$')
    return sorted(terminales)

# Nucleo LL(1)

class AnalizadorLL1:
    def __init__(self, reglas, simbolo_inicial):
        self.reglas = reglas
        self.inicio = simbolo_inicial
        self.no_terminales = set(reglas.keys())
        self.terminales = obtener_terminales(reglas)
        self.primeros = self._calcular_primeros()
        self.siguientes = self._calcular_siguientes()
        self.tabla, self.es_ll1 = self._construir_tabla()

    def _calcular_primeros(self):
        """Calcula conjuntos FIRST"""
        primeros = {nt: set() for nt in self.no_terminales}
        
        for t in self.terminales:
            primeros[t] = {t}
            
        cambio = True
        while cambio:
            cambio = False
            for nt in self.no_terminales:
                for prod in self.reglas[nt]:
                    i, epsilon = 0, True
                    while i < len(prod) and epsilon:
                        simb = prod[i]
                        f = primeros[simb] if simb in primeros else {simb}
                        
                        antes = len(primeros[nt])
                        primeros[nt].update(f - {'ε'})
                        
                        if len(primeros[nt]) != antes:
                            cambio = True
                            
                        if 'ε' not in f:
                            epsilon = False
                        i += 1
                        
                    if epsilon:
                        if 'ε' not in primeros[nt]:
                            primeros[nt].add('ε')
                            cambio = True
        return primeros

    def _calcular_siguientes(self):
        """Calcula conjuntos FOLLOW"""
        siguientes = {nt: set() for nt in self.no_terminales}
        siguientes[self.inicio].add('$')
        
        cambio = True
        while cambio:
            cambio = False
            for nt in self.no_terminales:
                for prod in self.reglas[nt]:
                    trailer = siguientes[nt].copy()
                    
                    for simb in reversed(prod):
                        if simb in self.no_terminales:
                            antes = len(siguientes[simb])
                            siguientes[simb].update(trailer)
                            
                            if len(siguientes[simb]) != antes: 
                                cambio = True
                            
                            if 'ε' in self.primeros[simb]:
                                trailer = trailer.union(self.primeros[simb] - {'ε'})
                            else:
                                trailer = self.primeros[simb] - {'ε'}
                        else:
                            trailer = {simb}
        return siguientes

    def _construir_tabla(self):
        """Construye tabla de analisis LL(1)"""
        tabla = {}
        es_valida = True
        
        for nt in self.no_terminales:
            for prod in self.reglas[nt]:
                primeros_prod = set()
                i, todo_epsilon = 0, True
                
                while i < len(prod) and todo_epsilon:
                    simb = prod[i]
                    f = self.primeros[simb] if simb in self.primeros else {simb}
                    primeros_prod |= (f - {'ε'})
                    
                    if 'ε' not in f: 
                        todo_epsilon = False
                    i += 1
                
                for t in primeros_prod:
                    if (nt, t) in tabla:
                        es_valida = False
                    tabla[(nt, t)] = prod
                
                if todo_epsilon:
                    for t in self.siguientes[nt]:
                        if (nt, t) in tabla:
                            es_valida = False
                        tabla[(nt, t)] = prod
        return tabla, es_valida

    def analizar(self, entrada):
        """Realiza analisis sintactico LL(1)"""
        if 'ε' in entrada:
            print("Error: 'ε' no permitido en entrada")
            return False

        tokens = list(entrada.strip()) + ['$']
        pila = ['$', self.inicio]
        i = 0
        secuencia = []
        traza = []

        print(f"\n{'Paso':<4} {'Pila':<20} {'Entrada':<20} {'Produccion'}")
        paso = 1

        while pila:
            pila_str = ' '.join(reversed(pila))
            entrada_str = ''.join(tokens[i:])
            tope = pila.pop()
            actual = tokens[i]
            prod_str = "-"

            if tope == actual:
                traza.append((paso, pila_str, entrada_str, prod_str))
                print(f"{paso:<4} {pila_str:<20} {entrada_str:<20} {prod_str}")
                i += 1
                if tope == '$':
                    break
                    
            elif tope in self.no_terminales:
                prod = self.tabla.get((tope, actual))
                
                if not prod:
                    print(f"Error: No hay produccion para ({tope}, {actual})")
                    print(f"Pila: {pila_str}  Entrada: {entrada_str}")
                    return False
                
                secuencia.append((tope, prod))
                prod_str = f"{tope} -> {' '.join(prod)}"
                traza.append((paso, pila_str, entrada_str, prod_str))
                print(f"{paso:<4} {pila_str:<20} {entrada_str:<20} {prod_str}")
                
                if prod != ['ε']:
                    for simb in reversed(prod):
                        pila.append(simb)
                        
            else:
                print(f"Error: '{tope}' no coincide con '{actual}'")
                print(f"Pila: {pila_str}  Entrada: {entrada_str}")
                return False
                
            paso += 1

        if (not pila and i == len(tokens)) or (not pila and i == len(tokens)-1 and tokens[i-1] == '$'):
            print("\nCadena valida")
            print("\nProducciones aplicadas:")
            for nt, prod in secuencia:
                print(f"{nt} -> {' '.join(prod)}")
            return True
        else:
            print("Error: analisis incompleto")
            return False

    # Metodos de visualizacion
    def mostrar_producciones(self):
        print("\n[PRODUCCIONES]")
        for nt in sorted(self.no_terminales):
            prods = []
            for prod in self.reglas[nt]:
                prods.append(' '.join(prod) if prod != ['ε'] else "ε")
            print(f"{nt} -> " + ' | '.join(prods))

    def mostrar_primeros_siguientes(self):
        print("\n[PRIMEROS]")
        for nt in sorted(self.no_terminales):
            print(f"FIRST({nt}): " + '{' + ', '.join(sorted(self.primeros[nt])) + '}')
            
        print("\n[SIGUIENTES]")
        for nt in sorted(self.no_terminales):
            print(f"FOLLOW({nt}): " + '{' + ', '.join(sorted(self.siguientes[nt])) + '}')

    def mostrar_tabla(self):
        print("\n[TABLA LL(1)]")
        columnas = self.terminales + ['$']
        filas = []
        
        for nt in sorted(self.no_terminales):
            fila = [nt]
            for t in columnas:
                entrada = self.tabla.get((nt, t))
                if entrada:
                    fila.append(' '.join(entrada))
                else:
                    fila.append('-')
            filas.append(fila)
            
        print(tabulate(filas, headers=['NT'] + columnas, tablefmt='fancy_grid'))

# Nucleo SLR(1)

def clausura(items, reglas):
    """Calcula la clausura LR(0)"""
    resultado = set(items)
    cambio = True
    
    while cambio:
        cambio = False
        nuevos_items = set()
        
        for (nt, rhs, punto) in resultado:
            if punto < len(rhs):
                B = rhs[punto]
                
                if B in reglas:
                    for prod in reglas[B]:
                        item = (B, tuple(prod), 0)
                        if item not in resultado:
                            nuevos_items.add(item)
                            
        if nuevos_items:
            resultado |= nuevos_items
            cambio = True
            
    return resultado

def transicion(items, X, reglas):
    """Calcula la transicion GOTO"""
    movidos = set()
    
    for (nt, rhs, punto) in items:
        if punto < len(rhs) and rhs[punto] == X:
            movidos.add((nt, rhs, punto+1))
            
    return clausura(movidos, reglas)

def construir_estados(reglas, simbolo_inicial):
    """Construye el automata LR(0)"""
    aumentado = f"{simbolo_inicial}'"
    while aumentado in reglas: 
        aumentado += "'"
        
    todas_reglas = {aumentado: [[simbolo_inicial]]}
    todas_reglas.update(reglas)
    
    inicial = clausura({(aumentado, tuple(todas_reglas[aumentado][0]), 0)}, todas_reglas)
    estados = [inicial]
    transiciones = {}
    
    simbolos = set(x for prods in todas_reglas.values() for prod in prods for x in prod) | set(todas_reglas.keys())
    
    while True:
        nuevos_estados = []
        
        for i, estado in enumerate(estados):
            for X in simbolos:
                siguiente = transicion(estado, X, todas_reglas)
                
                if siguiente and siguiente not in estados and siguiente not in nuevos_estados:
                    nuevos_estados.append(siguiente)
                    
                if siguiente:
                    idx = estados.index(siguiente) if siguiente in estados else len(estados) + nuevos_estados.index(siguiente)
                    transiciones[(i, X)] = idx
                    
        if not nuevos_estados:
            break
            
        estados += nuevos_estados
        
    return estados, transiciones, todas_reglas, aumentado

def calcular_primeros_slr(reglas):
    """Version SLR de FIRST"""
    primeros = {nt: set() for nt in reglas}
    
    for nt in reglas:
        for rhs in reglas[nt]:
            for simb in rhs:
                if simb not in reglas: 
                    primeros[simb] = {simb}
                    
    cambio = True
    while cambio:
        cambio = False
        
        for nt in reglas:
            for rhs in reglas[nt]:
                i, epsilon = 0, True
                
                while i < len(rhs) and epsilon:
                    simb = rhs[i]
                    f = primeros.get(simb, {simb})
                    
                    antes = len(primeros[nt])
                    primeros[nt].update(f - {'ε'})
                    
                    if len(primeros[nt]) != antes: 
                        cambio = True
                        
                    if 'ε' not in f: 
                        epsilon = False
                    i += 1
                    
                if epsilon:
                    if 'ε' not in primeros[nt]: 
                        primeros[nt].add('ε')
                        cambio = True
    return primeros

def calcular_siguientes_slr(reglas, primeros, simbolo_inicial):
    """Version SLR de FOLLOW"""
    siguientes = {nt: set() for nt in reglas}
    siguientes[simbolo_inicial].add('$')
    
    cambio = True
    while cambio:
        cambio = False
        
        for nt in reglas:
            for rhs in reglas[nt]:
                trailer = siguientes[nt].copy()
                
                for simb in reversed(rhs):
                    if simb in reglas:
                        antes = len(siguientes[simb])
                        siguientes[simb].update(trailer)
                        
                        if len(siguientes[simb]) != antes: 
                            cambio = True
                            
                        if 'ε' in primeros[simb]:
                            trailer = trailer.union(primeros[simb] - {'ε'})
                        else:
                            trailer = primeros[simb] - {'ε'}
                    else:
                        trailer = {simb}
    return siguientes

def construir_tabla_slr(reglas, simbolo_inicial):
    """Construye tablas ACTION/GOTO SLR"""
    estados, transiciones, todas_reglas, aumentado = construir_estados(reglas, simbolo_inicial)
    primeros = calcular_primeros_slr(todas_reglas)
    siguientes = calcular_siguientes_slr(todas_reglas, primeros, aumentado)
    
    ACTION = [{} for _ in range(len(estados))]
    GOTO = [{} for _ in range(len(estados))]
    
    mapa_prods = {}
    idx = 0
    lista_prods = []
    
    for nt in todas_reglas:
        for rhs in todas_reglas[nt]:
            mapa_prods[(nt, tuple(rhs))] = idx
            lista_prods.append((nt, rhs))
            idx += 1
            
    conflictos = []
    
    for i, I in enumerate(estados):
        for item in I:
            lhs, rhs, punto = item
            
            if punto < len(rhs):
                a = rhs[punto]
                
                if a not in todas_reglas:
                    j = transiciones.get((i, a))
                    
                    if j is not None:
                        if a in ACTION[i]:
                            conflictos.append((i, a, 'shift/shift'))
                            
                        ACTION[i][a] = ('s', j)
                        
            else:
                if lhs == aumentado:
                    ACTION[i]['$'] = ('acc',)
                else:
                    prod_idx = mapa_prods[(lhs, rhs)]
                    
                    for a in siguientes[lhs]:
                        if a in ACTION[i]:
                            conflictos.append((i, a, 'shift/reduce or reduce/reduce'))
                            
                        ACTION[i][a] = ('r', prod_idx)
                        
        for A in todas_reglas:
            if A in todas_reglas:
                j = transiciones.get((i, A))
                if j is not None:
                    GOTO[i][A] = j
                    
    es_slr = not conflictos
    
    return ACTION, GOTO, lista_prods, es_slr, estados, todas_reglas

def analizar_slr(entrada, ACTION, GOTO, lista_prods, simbolo_inicial):
    """Analizador SLR(1)"""
    tokens = list(entrada.strip()) + ['$']
    pila = [0]
    i = 0
    
    while True:
        estado = pila[-1]
        actual = tokens[i]
        
        accion = ACTION[estado].get(actual)
        if not accion:
            return False
            
        if accion[0] == 's':
            pila.append(actual)
            pila.append(accion[1])
            i += 1
                
        elif accion[0] == 'r':
            prod = lista_prods[accion[1]]
            lhs, rhs = prod
            
            if rhs != ['ε']:
                for _ in range(2*len(rhs)):
                    pila.pop()
                    
            estado = pila[-1]
            pila.append(lhs)
            pila.append(GOTO[estado][lhs])
                
        elif accion[0] == 'acc':
            return True
        else:
            return False

# Visualizacion SLR

def mostrar_tablas_slr(ACTION, GOTO, lista_prods, todas_reglas):
    """Muestra tablas SLR formateadas"""
    print("\n[TABLA ACTION]")
    simbolos = sorted(set(sym for a in ACTION for sym in a.keys() if sym != 'ε'))
    filas = []
    
    for i, a in enumerate(ACTION):
        fila = [i]
        for t in simbolos:
            act = a.get(t, '')
            
            if not act:
                fila.append('-')
            elif act[0] == 's':
                fila.append(f"s{act[1]}")
            elif act[0] == 'r':
                lhs, rhs = lista_prods[act[1]]
                rhs_str = ' '.join(rhs) if rhs != ['ε'] else 'ε'
                fila.append(f"r{act[1]}({lhs}->{rhs_str})")
            elif act[0] == 'acc':
                fila.append("acc")
            else:
                fila.append(str(act))
        filas.append(fila)
        
    print(tabulate(filas, headers=['st'] + simbolos, tablefmt='fancy_grid'))

    print("\n[TABLA GOTO]")
    nts = sorted(nt for nt in todas_reglas if nt.isupper())
    filas = []
    
    for i, g in enumerate(GOTO):
        fila = [i]
        for nt in nts:
            fila.append(str(g.get(nt,'-')))
        filas.append(fila)
        
    print(tabulate(filas, headers=['st']+nts, tablefmt='fancy_grid'))

    print("\n[PRODUCCIONES]")
    for i, (lhs, rhs) in enumerate(lista_prods):
        rhs_str = ' '.join(rhs) if rhs != ['ε'] else 'ε'
        print(f"r{i}: {lhs} -> {rhs_str}")

# Punto de Entrada

def principal():
    """Funcion principal del programa"""
    reglas, inicio = leer_gramatica()
    
    # Analisis LL(1)
    ll1 = AnalizadorLL1(reglas, inicio)
    
    # Analisis SLR(1)
    ACTION, GOTO, lista_prods, es_slr, estados, todas_reglas = construir_tabla_slr(reglas, inicio)
    es_ll1 = ll1.es_ll1

    # Determinar tipo de gramatica
    if es_ll1:
        print("Gramatica LL(1)")
        parser = "LL(1)"
    elif es_slr:
        print("Gramatica SLR(1)")
        parser = "SLR(1)"
        return
    else:
        print("Gramatica no LL(1) ni SLR(1)")
    
    # Menu interactivo
    while True:
        opcion = input("Seleccione parser (T: LL(1), B: SLR(1), Q: Salir): ").strip().upper()

        if opcion == 'Q':
            print("Saliendo...")
            break

        if opcion == 'T':
            if not es_ll1:
                print("Gramatica no LL(1). Elija otro parser.")
                continue
            else:
                print("Usando parser LL(1)")
                while True:
                    cadena = input()
                    if cadena == '':
                        break
                    resultado = ll1.analizar(cadena)
                    print("sí" if resultado else "no")

        elif opcion == 'B':
            if not es_slr:
                print("Gramatica no SLR(1). Elija otro parser.")
                continue
            else:
                print("Usando parser SLR(1)")
                while True:
                    cadena = input()
                    if cadena == '':
                        break
                    resultado = analizar_slr(cadena, ACTION, GOTO, lista_prods, inicio)
                    print("sí" if resultado else "no")
        else:
            print("Opcion invalida. Seleccione T, B o Q.")
    
if __name__ == "__main__":
    principal()
