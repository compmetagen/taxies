import csv

def read_protein_names(input_fn):
        
    protein_names = set()
    with open(input_fn, 'r') as input_handle:
        for row in csv.reader(input_handle, delimiter='\t'):
            if len(row) > 0:
                protein_names.add(row[0])
                
    return list(protein_names)
