import pandas as pd
import numpy as np
import datetime

def get_n(uni_list):
    n = []
    for uniprot in uni_list:
        df_pdbminer = pd.read_json(f"1core/results/{uniprot}/{uniprot}_all.json")
        n.append(len(df_pdbminer))
    n = sum(n)
    return n

uni_list = ["Q05397", "P62942", "P28698",
            "O75385", "P40692", "Q9GZQ8",
            "Q99856", "P10415", "P62937",
            "Q07817", "Q12931", "P48052"]

total_n = get_n(uni_list)

df = pd.DataFrame(columns=['cores', 'total', "per_pdb"])

directory_list = ["1core", "2core", "3core", "4core", "5core","6core"]

def read_core_runs(df, directory_list, n):
    for directory in directory_list:
        with open(f"{directory}/log.txt", "r") as f: 
            lines = f.readlines()
            t = lines[-2][:-1].split(": ")[-1]
            (h, m, s) = t.split(':')
            st = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))
            per_structure = st/n
            list_of_values = [directory, st, per_structure]
            df.loc[directory] = np.array(list_of_values, dtype="object")
    return df

df1 = read_core_runs(df, directory_list, total_n)

df1.to_csv("time.csv")
