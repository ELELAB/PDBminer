import pandas as pd
import numpy as np
import datetime
import glob


def read_commandline(logs, uni_dict):
    df = pd.DataFrame(columns=['cores','protein', 'total', "per_pdb", "n_structures"])
    for log in logs:
        gene = log.split("/")[-1].split("_")[0]
        uniprot = uni_dict[gene]    
        df_pdbminer = pd.read_json(f"commandline/results/{uniprot}/{uniprot}_all.json")
        n = len(df_pdbminer)
        with open(log, "r") as f: 
            lines = f.readlines()
            t = lines[3][:-1].split(": ")[-1]
            (h, m, s) = t.split(':')
            st = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))
            per_structure = st/n
            list_of_values = ["commandline", gene, st, per_structure, n]
            df.loc[gene] = np.array(list_of_values, dtype="object")
    df.to_csv("commandline_per_structure.csv")
    return df

logs = glob.glob("commandline/*.txt")
uni_dict = {"FAK1": "Q05397", "FKBP1A": "P62942", "MZF1": "P28698", 
            "ULK1": "O75385", "MHL1": "P40692", "MAP2LC3B": "Q9GZQ8", 
            "ARID3A": "Q99856", "BCL2": "P10415", "PPIA": "P62937", 
            "BCL2L1": "Q07817", "KRAS": "P01116", "TRAP1": "Q12931"}

command_df = read_commandline(logs, uni_dict)
total_n = sum(command_df.n_structures)

df1 = pd.DataFrame(columns=['cores', 'total', "per_pdb"])

command_df['total']=pd.to_timedelta(command_df['total'])
total_command = command_df['total'].sum()
command_df['per_pdb']=pd.to_timedelta(command_df['per_pdb'])
average_per_pdb = command_df['per_pdb'].mean()

list_of_values = ["commandline", total_command, average_per_pdb]

df1.loc["commandline"] = np.array(list_of_values, dtype="object")

directory_list = ["1core", "2core", "3core", "4core"]

def read_core_runs(df1, directory_list, n):
    for directory in directory_list:
        with open(f"{directory}/log.txt", "r") as f: 
            lines = f.readlines()
            t = lines[-2][:-1].split(": ")[-1]
            (h, m, s) = t.split(':')
            st = datetime.timedelta(hours=int(h), minutes=int(m), seconds=float(s))
            per_structure = st/n
            list_of_values = [directory, st, per_structure]
            df1.loc[directory] = np.array(list_of_values, dtype="object")
    return df1

df2 = read_core_runs(df1, directory_list, total_n)

df2.to_csv("time.csv")
