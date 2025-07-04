# -*- coding: utf-8 -*-
import requests
import Prog_Funcoes1
import Prog_Chama_GraficoN
import Prog_Chama_Grafico

arquivo = "C:/magre_dinamico/"

# ============================================================
arquivo1 = arquivo[0:17] + "/arquivos/pdb/"
arquivo2 = arquivo[0:17] + "/arquivos/testes/"
CodigoPDB = input("Informe codigo PDB : ")


pdbr = 'https://files.rcsb.org/view/' + CodigoPDB + '.pdb'

print("Codigo Solicitado:", CodigoPDB)

response = requests.get(pdbr)
print(response)
entrada = saida = arquivo1+CodigoPDB + ".pdb"
# print("resultado",saida)

A = response.text
proteinArq = ""
pdbLines = ""
A3DLines = ""

lin = 81
fim = int(len(A)/lin)
I = 0
chain = " "
contpdb = 0
desliga = 0
proteinArq = CodigoPDB + ".pdb"
chamador = 1
resultado = []

# entradax = arquivo1 + CodigoPDB + ".pdb"
# arq_entrada1 = open(entradax,'r')
# textoxx = arq_entrada1.readlines()

if response.status_code == 200:
    arq_saida = open(saida, 'w')
    # print("aaaaaaa")
    for I in range(0, fim-1, 1):
        words = A[lin*I:lin*(I+1)]
        arq_saida.write(words)
    arq_saida.close()

with open(entrada) as arq_entrada:
    pdbLines = arq_entrada.readlines()

for line in pdbLines:
    words = line
    # print(words)
    if len(words) < 1:
        pass
    elif (words[0:4] == "ATOM" and words[21:22] != chain):

        chain = words[21:22]
        print(proteinArq, pdbLines, A3DLines, chain, CodigoPDB, chamador)

        resultado = Prog_Funcoes1.movimenta(
            proteinArq, pdbLines, A3DLines, chain, CodigoPDB, chamador)

        resultado = Prog_Funcoes1.avalia_ruido(resultado)
        # print(resultado)
        PastaArq = arquivo2
        contatohpFilename = PastaArq + CodigoPDB[0:4] + "_" + chain + ".csv"
        with open(contatohpFilename, "w") as contatohpFile:
            for line in resultado:
                print(line, file=contatohpFile)
        Prog_Chama_GraficoN.Grafico(CodigoPDB, chain)
        Prog_Chama_Grafico.Grafico(CodigoPDB, chain)
        #Prog_Alinhamento.dash(CodigoPDB)
