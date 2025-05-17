# Importación de librerías necesarias
import sys  # Para manejo de sistema y entrada/salida
from collections import defaultdict  # Para usar diccionarios con valores por defecto
from tabulate import tabulate  # type: ignore # Para imprimir tablas formateadas

# --- Lectura de gramática ---
def input_grammar():
    # Lee el número de producciones de la gramática
    n = int(input())
    
    # Diccionario para almacenar las producciones (no-terminal -> lista de producciones)
    productions = defaultdict(list)
    
    # Lista para mantener el orden de los no-terminales
    non_terminals = []
    
    # Lee cada producción
    for _ in range(n):
        line = input().strip()
        if not line:  # Si la línea está vacía, continuar
            continue
            
        # Divide cada producción en parte izquierda (no-terminal) y derecha
        nt, rhs = line.split("->")
        nt = nt.strip()  # Elimina espacios en blanco
        
        # Agrega el no-terminal a la lista si no está
        if nt not in non_terminals:
            non_terminals.append(nt)
            
        # Procesa las alternativas de la parte derecha (separadas por espacios o |)
        alternatives = rhs.strip().split()
        
        for alt in alternatives:
            # Divide las producciones separadas por |
            for prod in alt.split('|'):
                # Reemplaza 'e' por ε (epsilon) y crea lista de símbolos
                symbols = [c if c != 'e' else 'ε' for c in prod.strip()]
                if not symbols:  # Si está vacío, es una producción epsilon
                    symbols = ['ε']
                # Agrega la producción al no-terminal correspondiente
                productions[nt].append(symbols)
                
    # Retorna las producciones y el primer no-terminal (símbolo inicial)
    return productions, non_terminals[0]

def get_terminals(productions):
    # Conjunto para almacenar terminales
    terms = set()
    
    # Recorre todas las producciones
    for rhss in productions.values():
        for rhs in rhss:
            for sym in rhs:
                # Un terminal es cualquier símbolo que no sea epsilon y no sea mayúscula
                if sym != 'ε' and (not sym.isupper()):
                    terms.add(sym)
                    
    # Elimina $ si está presente (símbolo de fin de cadena)
    terms.discard('$')
    
    # Retorna terminales ordenados
    return sorted(terms)

# --- Análisis LL(1) ---
class GrammarAnalyzerLL1:
    def __init__(self, productions, start_symbol):
        # Inicializa las estructuras básicas
        self.productions = productions  # Diccionario de producciones
        self.start = start_symbol  # Símbolo inicial
        self.non_terminals = set(productions.keys())  # Conjunto de no-terminales
        self.terminals = get_terminals(productions)  # Lista de terminales
        self.first = self._first()  # Calcula conjuntos FIRST
        self.follow = self._follow()  # Calcula conjuntos FOLLOW
        self.table, self.is_ll1 = self._ll1_table()  # Construye tabla LL(1)

    def _first(self):
        # Inicializa diccionario de conjuntos FIRST
        first = {nt: set() for nt in self.non_terminals}
        
        # FIRST de un terminal es el terminal mismo
        for t in self.terminals:
            first[t] = {t}
            
        changed = True
        while changed:
            changed = False
            for nt in self.non_terminals:
                for prod in self.productions[nt]:
                    i, add_epsilon = 0, True
                    # Calcula FIRST para cada símbolo en la producción
                    while i < len(prod) and add_epsilon:
                        sym = prod[i]
                        f = first[sym] if sym in first else {sym}
                        
                        # Guarda tamaño antes para detectar cambios
                        before = len(first[nt])
                        
                        # Agrega FIRST(sym) - {ε} a FIRST(nt)
                        first[nt].update(f - {'ε'})
                        
                        # Si hubo cambio, marca para continuar
                        if len(first[nt]) != before:
                            changed = True
                            
                        # Si ε no está en FIRST(sym), no continuar
                        if 'ε' not in f:
                            add_epsilon = False
                        i += 1
                        
                    # Si todos los símbolos pueden derivar ε, agregar ε
                    if add_epsilon:
                        if 'ε' not in first[nt]:
                            first[nt].add('ε')
                            changed = True
        return first

    def _follow(self):
        # Inicializa diccionario de conjuntos FOLLOW
        follow = {nt: set() for nt in self.non_terminals}
        
        # FOLLOW del símbolo inicial contiene $
        follow[self.start].add('$')
        
        changed = True
        while changed:
            changed = False
            for nt in self.non_terminals:
                for prod in self.productions[nt]:
                    trailer = follow[nt].copy()
                    
                    # Procesa la producción de derecha a izquierda
                    for sym in reversed(prod):
                        if sym in self.non_terminals:
                            before = len(follow[sym])
                            
                            # Agrega trailer al FOLLOW del no-terminal
                            follow[sym].update(trailer)
                            
                            if len(follow[sym]) != before: 
                                changed = True
                            
                            # Si el símbolo puede derivar ε, agrega FIRST(sym)-{ε} al trailer
                            if 'ε' in self.first[sym]:
                                trailer = trailer.union(self.first[sym] - {'ε'})
                            else:
                                trailer = self.first[sym] - {'ε'}
                        else:
                            # Para terminales, el trailer es el terminal mismo
                            trailer = {sym}
        return follow

    def _ll1_table(self):
        table = {}  # Tabla de análisis LL(1)
        is_ll1 = True  # Bandera para gramática LL(1)
        
        for nt in self.non_terminals:
            for prod in self.productions[nt]:
                first_set = set()
                i, all_epsilon = 0, True
                
                # Calcula FIRST de la producción
                while i < len(prod) and all_epsilon:
                    sym = prod[i]
                    f = self.first[sym] if sym in self.first else {sym}
                    first_set |= (f - {'ε'})
                    
                    # Si un símbolo no deriva ε, termina
                    if 'ε' not in f: 
                        all_epsilon = False
                    i += 1
                
                # Para cada terminal en FIRST, agrega a la tabla
                for t in first_set:
                    if (nt, t) in table:  # Conflicto: no es LL(1)
                        is_ll1 = False
                    table[(nt, t)] = prod
                
                # Si toda la producción deriva ε, usa FOLLOW
                if all_epsilon:
                    for t in self.follow[nt]:
                        if (nt, t) in table:  # Conflicto: no es LL(1)
                            is_ll1 = False
                        table[(nt, t)] = prod
        return table, is_ll1

    def parse(self, inp):
        # Verifica que la cadena no contenga ε explícito
        if 'ε' in inp:
            print("Error: la cadena de entrada no debe contener 'ε'. Use una cadena vacía si corresponde.")
            return False

        # Prepara la cadena de entrada (añade $ al final)
        tokens = list(inp.strip()) + ['$']
        
        # Inicializa pila con $ y símbolo inicial
        stack = ['$', self.start]
        i = 0  # Índice para recorrer la entrada
        
        # Estructuras para guardar información del análisis
        sequence = []  # Secuencia de producciones aplicadas
        trace = []     # Traza de pasos del análisis

        # Encabezado para la traza
        print(f"\n{'Paso':<4} {'Pila':<20} {'Entrada':<20} {'Producción'}")
        paso = 1  # Contador de pasos

        while stack:
            # Prepara strings para mostrar estado actual
            pila_str = ' '.join(reversed(stack))
            entrada_str = ''.join(tokens[i:])
            top = stack.pop()  # Saca elemento de la pila
            cur = tokens[i]    # Símbolo actual de entrada
            prod_str = "-"      # Producción aplicada (inicialmente ninguna)

            if top == cur:  # Coincidencia de terminales
                trace.append((paso, pila_str, entrada_str, prod_str))
                print(f"{paso:<4} {pila_str:<20} {entrada_str:<20} {prod_str}")
                i += 1
                if top == '$':  # Fin del análisis
                    break
                    
            elif top in self.non_terminals:  # No-terminal en la pila
                # Busca producción en la tabla
                prod = self.table.get((top, cur))
                
                if not prod:  # No hay producción válida
                    print(f"Error: No se encontró entrada en la tabla para ({top}, {cur})")
                    print(f"Pila: {pila_str}  Entrada: {entrada_str}")
                    return False
                
                # Guarda la producción aplicada
                sequence.append((top, prod))
                prod_str = f"{top} → {' '.join(prod)}"
                trace.append((paso, pila_str, entrada_str, prod_str))
                print(f"{paso:<4} {pila_str:<20} {entrada_str:<20} {prod_str}")
                
                # Si no es producción epsilon, pone símbolos en pila en orden inverso
                if prod != ['ε']:
                    for sym in reversed(prod):
                        stack.append(sym)
                        
            else:  # Error: terminal en pila no coincide con entrada
                print(f"Error: símbolo en pila '{top}' no coincide con entrada '{cur}'")
                print(f"Pila: {pila_str}  Entrada: {entrada_str}")
                return False
                
            paso += 1  # Incrementa contador de pasos

        # Verifica condiciones de aceptación
        pila_final = (stack == [])
        entrada_final = (i == len(tokens))
        
        if pila_final and i == len(tokens):
            print("\nLa cadena pertenece al lenguaje.")
            print("\nSecuencia de producción utilizada:")
            for nt, prod in sequence:
                print(f"{nt} → {' '.join(prod)}")
            return True
            
        elif pila_final and i == len(tokens)-1 and tokens[i-1] == '$':
            print("\nLa cadena pertenece al lenguaje.")
            print("\nSecuencia de producción utilizada:")
            for nt, prod in sequence:
                print(f"{nt} → {' '.join(prod)}")
            return True
            
        else:
            print("Error: la cadena no fue completamente consumida o la pila no está vacía.")
            return False

    # ---------- Métodos para imprimir tablas ----------
    def print_productions(self):
        print("\n[Producciones]")
        for nt in sorted(self.non_terminals):
            prods = []
            for prod in self.productions[nt]:
                # Formatea producción (epsilon si es vacía)
                prods.append(' '.join(prod) if prod != ['ε'] else "ε")
            print(f"{nt} → " + ' | '.join(prods))

    def print_first_follow(self):
        print("\n[FIRST]")
        for nt in sorted(self.non_terminals):
            print(f"FIRST({nt}): " + '{' + ', '.join(sorted(self.first[nt])) + '}')
            
        print("\n[FOLLOW]")
        for nt in sorted(self.non_terminals):
            print(f"FOLLOW({nt}): " + '{' + ', '.join(sorted(self.follow[nt])) + '}')

    def print_ll1_table(self):
        print("\n[TABLA LL(1)]")
        col_terms = self.terminals + ['$']  # Columnas: terminales + $
        rows = []
        
        for nt in sorted(self.non_terminals):
            row = [nt]
            for t in col_terms:
                entry = self.table.get((nt, t))
                if entry:
                    row.append(' '.join(entry))  # Producción correspondiente
                else:
                    row.append('-')  # Celda vacía
            rows.append(row)
            
        # Imprime tabla formateada
        print(tabulate(rows, headers=['NT'] + col_terms, tablefmt='fancy_grid'))

# --- Análisis SLR(1) ---
def closure(items, productions):
    # Calcula la clausura de un conjunto de items LR(0)
    result = set(items)
    changed = True
    
    while changed:
        changed = False
        new_items = set()
        
        for (nt, rhs, dot) in result:
            if dot < len(rhs):  # Si el punto no está al final
                B = rhs[dot]    # Símbolo después del punto
                
                if B in productions:  # Si es no-terminal
                    for prod in productions[B]:
                        # Agrega items con punto al inicio
                        item = (B, tuple(prod), 0)
                        if item not in result:
                            new_items.add(item)
                            
        if new_items:
            result |= new_items
            changed = True
            
    return result

def goto(items, X, productions):
    # Calcula la función GOTO para un símbolo X
    moved = set()
    
    for (nt, rhs, dot) in items:
        if dot < len(rhs) and rhs[dot] == X:
            # Mueve el punto después de X
            moved.add((nt, rhs, dot+1))
            
    return closure(moved, productions)  # Retorna clausura del resultado

def build_states(productions, start_symbol):
    # Aumenta la gramática con un nuevo símbolo inicial
    augmented = f"{start_symbol}'"
    while augmented in productions: 
        augmented += "'"
        
    # Crea gramática aumentada
    all_prods = {augmented: [[start_symbol]]}
    all_prods.update(productions)
    
    # Estado inicial: clausura del item inicial
    initial = closure({(augmented, tuple(all_prods[augmented][0]), 0)}, all_prods)
    states = [initial]  # Lista de estados
    transitions = {}    # Diccionario de transiciones
    
    # Todos los símbolos (terminales y no-terminales)
    symbols = set(x for prods in all_prods.values() for prod in prods for x in prod) | set(all_prods.keys())
    
    # Construcción de los estados
    while True:
        new_states = []
        
        for i, state in enumerate(states):
            for X in symbols:
                nxt = goto(state, X, all_prods)
                
                if nxt and nxt not in states and nxt not in new_states:
                    new_states.append(nxt)
                    
                if nxt:
                    # Guarda transición (estado, símbolo) -> nuevo estado
                    idx = states.index(nxt) if nxt in states else len(states) + new_states.index(nxt)
                    transitions[(i, X)] = idx
                    
        if not new_states:  # No hay nuevos estados
            break
            
        states += new_states  # Agrega nuevos estados
        
    return states, transitions, all_prods, augmented

def compute_first(productions):
    # Calcula conjuntos FIRST para SLR (similar a LL(1))
    first = {nt: set() for nt in productions}
    
    for nt in productions:
        for rhs in productions[nt]:
            for sym in rhs:
                if sym not in productions: 
                    first[sym] = {sym}
                    
    changed = True
    while changed:
        changed = False
        
        for nt in productions:
            for rhs in productions[nt]:
                i, add_epsilon = 0, True
                
                while i < len(rhs) and add_epsilon:
                    sym = rhs[i]
                    f = first.get(sym, {sym})
                    
                    before = len(first[nt])
                    first[nt].update(f - {'ε'})
                    
                    if len(first[nt]) != before: 
                        changed = True
                        
                    if 'ε' not in f: 
                        add_epsilon = False
                    i += 1
                    
                if add_epsilon:
                    if 'ε' not in first[nt]: 
                        first[nt].add('ε')
                        changed = True
    return first

def compute_follow(productions, first, start_symbol):
    # Calcula conjuntos FOLLOW para SLR (similar a LL(1))
    follow = {nt: set() for nt in productions}
    follow[start_symbol].add('$')
    
    changed = True
    while changed:
        changed = False
        
        for nt in productions:
            for rhs in productions[nt]:
                trailer = follow[nt].copy()
                
                for sym in reversed(rhs):
                    if sym in productions:
                        before = len(follow[sym])
                        follow[sym].update(trailer)
                        
                        if len(follow[sym]) != before: 
                            changed = True
                            
                        if 'ε' in first[sym]:
                            trailer = trailer.union(first[sym] - {'ε'})
                        else:
                            trailer = first[sym] - {'ε'}
                    else:
                        trailer = {sym}
    return follow

def build_slr_table(productions, start_symbol):
    # Construye la tabla de análisis SLR
    states, transitions, all_prods, augmented = build_states(productions, start_symbol)
    first = compute_first(all_prods)
    follow = compute_follow(all_prods, first, augmented)
    
    # Inicializa tablas ACTION y GOTO
    ACTION = [{} for _ in range(len(states))]
    GOTO = [{} for _ in range(len(states))]
    
    # Mapeo de producciones a índices
    prod_map = {}
    idx = 0
    prod_list = []
    
    for nt in all_prods:
        for rhs in all_prods[nt]:
            prod_map[(nt, tuple(rhs))] = idx
            prod_list.append((nt, rhs))
            idx += 1
            
    conflicts = []  # Para detectar conflictos
    
    for i, I in enumerate(states):
        for item in I:
            lhs, rhs, dot = item
            
            if dot < len(rhs):  # Item con punto no al final
                a = rhs[dot]    # Símbolo después del punto
                
                if a not in all_prods:  # Si es terminal
                    j = transitions.get((i, a))
                    
                    if j is not None:
                        if a in ACTION[i]:  # Conflicto shift-shift
                            conflicts.append((i, a, 'shift/shift'))
                            
                        # Acción de desplazamiento
                        ACTION[i][a] = ('s', j)
                        
            else:  # Item con punto al final
                if lhs == augmented:  # Item de aceptación
                    ACTION[i]['$'] = ('acc',)
                else:
                    prod_idx = prod_map[(lhs, rhs)]
                    
                    # Para cada terminal en FOLLOW(lhs), agregar reducción
                    for a in follow[lhs]:
                        if a in ACTION[i]:  # Conflicto shift-reduce o reduce-reduce
                            conflicts.append((i, a, 'shift/reduce or reduce/reduce'))
                            
                        ACTION[i][a] = ('r', prod_idx)
                        
        # Llenar tabla GOTO para no-terminales
        for A in all_prods:
            if A in all_prods:
                j = transitions.get((i, A))
                if j is not None:
                    GOTO[i][A] = j
                    
    # La gramática es SLR(1) si no hay conflictos
    is_slr = not conflicts
    
    return ACTION, GOTO, prod_list, is_slr, states, all_prods

def slr_parse(inp, ACTION, GOTO, prod_list, start_symbol):
    # Analizador SLR
    tokens = list(inp.strip()) + ['$']  # Prepara entrada
    stack = [0]  # Pila de estados
    i = 0        # Índice de entrada
    
    while True:
        state = stack[-1]  # Estado actual
        cur = tokens[i]    # Símbolo actual
        
        # Obtiene acción de la tabla
        action = ACTION[state].get(cur)
        if not action:  # Error: no hay acción definida
            return False
            
        if action[0] == 's':  # Desplazamiento
            stack.append(cur)     # Empuja símbolo
            stack.append(action[1])  # Empuja nuevo estado
            i += 1               # Avanza en la entrada
                
        elif action[0] == 'r':  # Reducción
            prod = prod_list[action[1]]  # Producción a reducir
            lhs, rhs = prod
            
            if rhs != ['ε']:  # Si no es producción epsilon
                # Desapila 2*|rhs| elementos (símbolos y estados)
                for _ in range(2*len(rhs)):
                    stack.pop()
                    
            state = stack[-1]  # Estado después de desapilar
            stack.append(lhs)   # Empuja cabeza de producción
            stack.append(GOTO[state][lhs])  # Empuja nuevo estado según GOTO
                
        elif action[0] == 'acc':  # Aceptación
            return True
        else:  # Error
            return False

# --- Impresión de tablas SLR ---
def print_slr_tables(ACTION, GOTO, prod_list, all_prods):
    print("\n[TABLA ACTION]")
    # Filtra ε y elimina duplicados en columnas
    all_syms = sorted(set(sym for a in ACTION for sym in a.keys() if sym != 'ε'))
    rows = []
    
    for i, a in enumerate(ACTION):
        row = [i]
        for t in all_syms:
            act = a.get(t, '')
            
            if not act:
                row.append('-')
            elif act[0] == 's':  # Desplazamiento
                row.append(f"s{act[1]}")
            elif act[0] == 'r':  # Reducción
                lhs, rhs = prod_list[act[1]]
                rhs_str = ' '.join(rhs) if rhs != ['ε'] else 'ε'
                row.append(f"r{act[1]}({lhs}→{rhs_str})")
            elif act[0] == 'acc':  # Aceptación
                row.append("acc")
            else:
                row.append(str(act))
        rows.append(row)
        
    print(tabulate(rows, headers=['st'] + all_syms, tablefmt='fancy_grid'))

    print("\n[TABLA GOTO]")
    nts = sorted(nt for nt in all_prods if nt.isupper())  # No-terminales
    rows = []
    
    for i, g in enumerate(GOTO):
        row = [i]
        for nt in nts:
            row.append(str(g.get(nt,'-')))  # Estado destino o -
        rows.append(row)
        
    print(tabulate(rows, headers=['st']+nts, tablefmt='fancy_grid'))

    print("\n[PRODUCCIONES] (índices para reduce)")
    for i, (lhs, rhs) in enumerate(prod_list):
        rhs_str = ' '.join(rhs) if rhs != ['ε'] else 'ε'
        print(f"r{i}: {lhs} → {rhs_str}")

def main():
    # Punto de entrada principal
    productions, start = input_grammar()  # Lee gramática
    
    # Análisis LL(1)
    ll1 = GrammarAnalyzerLL1(productions, start)
    
    # Análisis SLR(1)
    ACTION, GOTO, prod_list, is_slr, states, all_prods = build_slr_table(productions, start)
    is_ll1 = ll1.is_ll1

    # Determina tipo de gramática
    if is_ll1:
        print("Grammar is LL(1).")
        parser_type = "LL(1)"
    elif is_slr:
        print("Grammar is SLR(1).")
        parser_type = "SLR(1)"
        return
    else:
        print("Grammar is neither LL(1) nor SLR(1).")
    
    # Menú interactivo
    while True:
        choice = input("Select a parser (T: for LL(1), B: for SLR(1), Q: quit): ").strip().upper()

        if choice == 'Q':
            print("Assume Q is given.")
            break

        if choice == 'T':  # Analizador LL(1)
            if not is_ll1:
                print("Grammar is not LL(1). Please choose another parser.")
                continue
            else:
                print("Using LL(1) parser.")
                while True:
                    cadena = input()
                    if cadena == '':
                        break
                    result = ll1.parse(cadena)
                    print("yes" if result else "no")

        elif choice == 'B':  # Analizador SLR(1)
            if not is_slr:
                print("Grammar is not SLR(1). Please choose another parser.")
                continue
            else:
                print("Using SLR(1) parser.")
                while True:
                    cadena = input()
                    if cadena == '':
                        break
                    result = slr_parse(cadena, ACTION, GOTO, prod_list, start)
                    print("yes" if result else "no")
        else:
            print("Invalid option. Please select T, B or Q.")
    
if __name__ == "__main__":
    main()
