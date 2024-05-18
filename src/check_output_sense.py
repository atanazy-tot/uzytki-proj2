def check(output_path: str) -> dict:
    
    Variables = {}
    
    with open(output_path, 'r') as f:
        lines = f.read().splitlines()
        for i, line in enumerate(lines):
            split_line = line.split()
            if len(split_line) > 1:
                if split_line[0] == '1':
                    start_index = i
                if len(split_line[1]) > 1:
                    if split_line[1][1] == 'z':
                        end_index = i
                        break
                    
        # tylko te linie naprawdę zawierają istotne dane, tzn. zmienne 'x'        
        lines = lines[start_index:end_index]
        
        for i, line in enumerate(lines):
            temp_line = line.split()
            variable = temp_line[1].replace("\'", "")[3:-2].split(',')
            variable = tuple([int(x[1:]) for x in variable])
            value = int(temp_line[2])
            
            Variables[variable] = value
            
    return Variables


def timetable(Variables: dict,
              N: int,
              J: int,
              K: int,
              timetable_path: str,
              employer_demands: str,
              variant: str = 'employer',
              employee_id: int | None = None):
    
    if variant not in ('employee', 'employer'):
        raise ValueError("Argument 'variant' must be one of 'employee', 'employer'")
    
    if variant == 'employer':
        temp = {}
        for x in Variables:
            if Variables[x] == 0:
                continue
            if temp.get(x[0]) is None:
                temp[x[0]] = [[x[1], x[2]]]
            else:
                temp[x[0]].append([x[1], x[2]])
        time = ';'+';'.join([f'{k};;' for k in range(1, K+1)])
        # print(time)
        shifts = [0 for t in range(K*J)]
        for i in range(1, N+1):
            # mamy w ręku konkretną pielęgniarkę
            work = []
            for k in list(range(1, K+1)):
                # mamy w ręku konkretny dzień
                for j in range(1, J+1):
                    # mamy w ręku konkretną zmianę
                    
                    if temp.get(i) is not None and [j, k] in temp[i]:
                        # pielęgniarka 'i' ma przyjść tego dnia na tę zmianę
                        work.append('1')
                        shifts[(k-1)*J + (j-1)] += 1
                    else:
                        work.append('0')
            work = f'Nurse_{i};'+';'.join(work)
            temp[i] = work
            # print(temp[i])
            
        with open(timetable_path, 'w') as f:
            print(time, file=f)
            print((';'+';'.join([f'{j+1}' for j in range(J)]))*K, file=f)
            for i in range(1, N+1):
                print(temp[i], file=f)
            print('Nurses_on_shift;'+';'.join([f'{k}' for k in shifts]), file=f)
            print('Demanded_nurses'+employer_demands, file=f)
                
    else:
        # wariant z grafikiem zawierającym tylko jednego pracownika
        pass
        
        
    
                
                