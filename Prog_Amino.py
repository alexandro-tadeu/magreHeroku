import requests
from bs4 import BeautifulSoup

class Atom:
    def __init__(self, line):
        self.atom_name = line[12:16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21]
        self.res_number = int(line[22:26])
        self.coordinates = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54])
        )

def extrair_sequencia(pdb_lines, chain):
    sequencia = ""

    for line in pdb_lines:
        if line.startswith("ATOM") and line[21] == chain:
            atom = Atom(line)
            sequencia += atom.res_name

    return sequencia

def obter_codigo_pdb(codigo_pdb):
    url = f'https://files.rcsb.org/view/{codigo_pdb}.pdb'
    response = requests.get(url)

    if response.status_code == 200:
        return response.text.split('\n')
    else:
        print(f"Erro ao obter o código PDB {codigo_pdb}. Verifique o código e tente novamente.")
        return None

def main():
    codigo_pdb = input("Insira o código PDB: ").upper()
    pdb_lines = obter_codigo_pdb(codigo_pdb)

    if pdb_lines:
        # Substitua 'A' pela cadeia desejada
        sequencia_aminoacidos = extrair_sequencia(pdb_lines, 'A')

        print("Sequência de Aminoácidos:", sequencia_aminoacidos)

if __name__ == "__main__":
    main()
