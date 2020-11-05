import requests as r

endpoint = "http://biosig.unimelb.edu.au/mcsm/stability_prediction_list"

with open("2ocj.pdb", 'rb') as pdb_file:
    data_pdb = pdb_file.read()

with open("mutation_list", 'rb') as mut:
    data_mut = mut.read()

data_val = {
    "wild": data_pdb,
    "mutation_list": [data_mut]
}
headers = {
'content-type' : 'application/octet-stream'
}
files = {'file': open("2ocj.pdb", 'rb')}
values = {"mutation_list": data_mut}
x = r.post(endpoint, data=data_val)
print(x.text)
