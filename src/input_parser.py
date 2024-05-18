"""
Parse file with input data into AMPL .dat format
"""
import csv
def data_to_dict(input_path: str) -> dict:
    Data = {}
    with open(input_path, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=';')
        for line in reader:
            # print(line)
            if len(line[0]) > 2:
                # N na pewno nie będzie liczbą trzycyfrową
                # dlatego zakładamy, że jeśli line[0] ma długość >2
                # to jest tytułem nowej grupy danych do programu
                data_type = line[0]
                Data[data_type] = []
            elif Data.get(data_type) is not None:
                # print(line)
                Data[data_type].append([int(el) for el in line if el != ''])
    return Data

def merged(lst: list[list]):
    temp_dict = {}
    retlist = []
    for l in lst:
        key, val = l[0], l[1:]
        if temp_dict.get(key) is None:
            temp_dict[key] = val
        else:
            temp_dict[key].extend(val)  
    # print(temp_dict)
    for key, val in temp_dict.items():
        temp = [key]
        temp.extend(val)
        retlist.append(temp)
    return retlist

def transform_data(data: dict) -> dict:
    temp_data = {}
    for key, val in data.items():
        temp_data[key] = merged(val)
    return temp_data
                
def generate_dat(input_path: str, output_path: str = None) -> dict: # -> tuple:
    
    Data_raw = data_to_dict(input_path) # wszystkie dane zostały wczytane do tego słownika
    Data = transform_data(Data_raw)

    # for key, val in Data.items():
    #     print(f'{key=}\n{val=}\n\n')
    
    with open(output_path, 'w') as output_file:    
        N = len(list(Data.values())[1]) # liczba pielęgniarek
        K = len(list(Data.values())[0]) # liczba dni do zaplanowania
        J = len(list(Data.values())[0][0]) - 1 # liczba zmian dziennie
        
        print('data;\n', file=output_file)
        
        nurses = ' '.join(['n'+f'{i}' for i in range(1, N+1)])
        print(f'set NURSES := {nurses};', file=output_file)
        
        shifts = ' '.join(['s'+f'{j}' for j in range(1, J+1)])
        print(f'set SHIFTS := {shifts};', file=output_file)
        
        days = ' '.join(['d'+f'{k}' for k in range(1, K+1)])
        print(f'set DAYS := {days};\n', file=output_file)
        
        print(f'param N := {N};', file=output_file)
        print(f'param J := {J};', file=output_file)
        print(f'param K := {K};\n', file=output_file)
        
        print(f'param demand (tr): {days} :=', file=output_file)
        for j in range(J):
            print(f's{j+1}', file=output_file, end = '')
            for day in Data['demand']:
                print(f" {day[j+1]}", file=output_file, end = '')
            print('', file=output_file) if j < J-1 else print(';', file=output_file)
        
        print(f'\nparam workhours :=', file=output_file)
        for i in range(N):
            print(f'n{i+1} ', file=output_file, end = '')
            workhours = Data['workhours'][i][1:]
            if workhours == []: # jeśli nie ma godzin do przepracowania
                l = ' 0'
            else: # jeśli ma ileś godzin do przepracowania
                workhours.sort()
                l = ' '+str(workhours[0])
            print(l, file=output_file, end='')
            print('', file=output_file) if i < N-1 else print(';', file=output_file)
            
        print(f'\nparam vacation (tr): {days} :=', file=output_file)
        for i in range(N):
            print(f'n{i+1} ', file=output_file, end = '')
            vac_days = Data['vacation'][i][1:]
            if vac_days == []: # jeśli 'i' nie wzięła żadnego dnia wolnego
                l = ' 0'*K
            else: # jeśli wzięła przynajmniej jeden dzień wolny
                # print(f'{vac_days=}')
                vac_days.sort()
                l = ' 0'*(vac_days[0] - 1) + ' 1'
                for i in range(1, len(vac_days)):
                    l += ' 0'*(vac_days[i] - vac_days[i-1] - 1) + ' 1'
                l += ' 0'*(K - vac_days[-1])
            print(l, file=output_file, end='')
            print('', file=output_file) if i < N-1 else print(';', file=output_file)
            
        print(f'\nparam pref_comp (tr): {nurses} :=', file=output_file)
        for i in range(N):
            print(f'n{i+1} ', file=output_file, end = '')
            pref_comp = Data['preferred_companions'][i][1:]
            if pref_comp == []: # jeśli 'i' nie preferuje żadnej koleżanki
                l = ' 0'*K
            else: # jeśli preferuje przynajmniej jedną
                pref_comp.sort()
                l = ' 0'*(pref_comp[0] - 1) + ' 1'
                for i in range(1, len(pref_comp)):
                    l += ' 0'*(pref_comp[i] - pref_comp[i-1] - 1) + ' 1'
                l += ' 0'*(N - pref_comp[-1])
            print(l, file=output_file, end='')
            print('', file=output_file) if i < N-1 else print(';', file=output_file)
            
        print(f'\nparam unpref_comp (tr): {nurses} :=', file=output_file)
        for i in range(N):
            print(f'n{i+1} ', file=output_file, end = '')
            unpref_comp = Data['unpreferred_companions'][i][1:]
            if unpref_comp == []: # jeśli 'i' nie preferuje żadnej koleżanki
                l = ' 0'*K
            else: # jeśli preferuje przynajmniej jedną
                unpref_comp.sort()
                l = ' 0'*(unpref_comp[0] - 1) + ' 1'
                for i in range(1, len(unpref_comp)):
                    l += ' 0'*(unpref_comp[i] - unpref_comp[i-1] - 1) + ' 1'
                l += ' 0'*(N - unpref_comp[-1])
            print(l, file=output_file, end='')
            print('', file=output_file) if i < N-1 else print(';', file=output_file)
            
        print(f'\nparam pref_shifts :=', file=output_file)
        for j in range(J):
            print(f'\n[*,*,s{j+1}]: {days} :=', file=output_file)
            for i in range(N):
                print(f'n{i+1} ', file=output_file, end = '')
                
                pref_shifts = Data['preferred_shifts'][i][1:]
                ks = [x for i, x in enumerate(pref_shifts) if i % 2 == 0]
                ss = [x for i, x in enumerate(pref_shifts) if i % 2 == 1]
                
                indices = [i for i, x in enumerate(ss) if x == j+1]
                # print(f'{ks=}\n{ss=}')
                
                if indices == []:
                    l = ' 0'*K
                else:
                    ks = [ks[x] for x in indices] # lista dni dot. zmiany j
                    ks.sort()
                    l = ' 0'*(ks[0] - 1) + ' 1'
                    for i in range(1, len(ks)):
                        l += ' 0'*(ks[i] - ks[i-1] - 1) + ' 1'
                    l += ' 0'*(K - ks[-1])
                    # print(f'{ks=}')
                    # for s in ss:
                    #     if s != j+1:
                    #         l = ' 0'*K
                    #     else:
                    #         l = ' 0'*(ks[0] - 1)
                    #         if len(ks) == 1:
                    #             l += ' 1' + ' 0'*(K-ks[0])
                    #         else:
                    #             # przypadek większej liczby deklaracji
                    #             # zostaje do implementacji!
                    #             pass
                            
                print(l, file=output_file, end='')
                print('', file=output_file) if i < N-1 or j < J-1 else print(';', file=output_file)
                
        print(f'\nparam unpref_shifts :=', file=output_file)     
        for j in range(J):
            print(f'\n[*,*,s{j+1}]: {days} :=', file=output_file)
            for i in range(N):
                print(f'n{i+1} ', file=output_file, end = '')
                
                unpref_shifts = Data['unpreferred_shifts'][i][1:]
                ks = [x for i, x in enumerate(unpref_shifts) if i % 2 == 0]
                ss = [x for i, x in enumerate(unpref_shifts) if i % 2 == 1]
                
                ks.sort()
                ss.sort()
                
                # print(f'{ks=}\n{ss=}')
                
                if ks == [] or ss == [] :
                    l = ' 0'*K
                else:
                    for s in ss:
                        if s != j+1:
                            l = ' 0'*K
                        else:
                            l = ' 0'*(ks[0] - 1)
                            if len(ks) == 1:
                                l += ' 1' + ' 0'*(K-ks[0])
                            else:
                                # przypadek większej liczby deklaracji
                                # zostaje do implementacji!
                                pass
                            
                print(l, file=output_file, end='')
                print('', file=output_file) if i < N-1 or j < J-1 else print(';', file=output_file)
    return N, J, K

def employer_demands(input_path: str, J: int, K: int) -> str:
    demands = data_to_dict(input_path)['demand']
    # print(demands)
    time = ';'+';'.join([f'{k+1};;' for k in range(K)])
    # print(time)
    
    temp = []
    for k in demands:
        pass
        # k oznacza jak zwykle jeden dzień
        temp.extend([f'{x}' for x in k[1:]])
    # print(temp)
    temp = ';'+';'.join(temp)
    return time, temp