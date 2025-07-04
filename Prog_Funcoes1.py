# -*- coding: utf-8 -*-
import os
import Prog_Mod
from atom import Atom

import matplotlib.pyplot as plt
import freesasa
import numpy as np
import pandas as pd
import Prog_Predicao

# Diretórios seguros e dinâmicos
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ARQUIVOS_DIR = os.path.join(BASE_DIR, 'arquivos')
PDB_DIR = os.path.join(ARQUIVOS_DIR, 'pdb')
SAIDA_DIR = os.path.join(ARQUIVOS_DIR, 'saida')
ESPECIAIS_DIR = os.path.join(ARQUIVOS_DIR, 'especiais')


def plot_2d_space(X, y, label='Classes'):
    colors = ['#1F77B4', '#FF7F0E']
    markers = ['o', 's']
    for l, c, m in zip(np.unique(y), colors, markers):
        plt.scatter(
            X[y == l, 0],
            X[y == l, 1],
            c=c, label=l, marker=m
        )
    plt.title(label)
    plt.legend(loc='upper right')
    plt.show()


def avalia_ruido(resultado):
    resultadox = []

    # Se resultado muito pequeno, retorna cópia direta
    if len(resultado) < 3:
        return resultado.copy()

    # Percorre do segundo até o penúltimo elemento para evitar índices inválidos
    for i in range(1, len(resultado) - 1):
        resultado11 = resultado[i-1].split(";")
        resultado12 = resultado[i].split(";")
        resultado13 = resultado[i+1].split(";")
        if float(resultado12[3]) > 0.5:
            if float(resultado11[3]) < 0.5 and float(resultado13[3]) < 0.5:
                resultados = resultado12[0] + ";" + resultado12[1] + ";" + resultado12[2] + ";"
                resultados += str(float(resultado12[3]) / 2) + ";" + resultado12[4] + ";"
                resultadox.append(resultados)
            else:
                resultadox.append(resultado[i])
        else:
            resultadox.append(resultado[i])

    # Mantém o primeiro e o último elemento intactos, adicionando-os nas bordas
    resultadox.insert(0, resultado[0])
    resultadox.append(resultado[-1])

    return resultadox


def calculascore(entrada, saida):
    contador = 0
    with open(entrada, 'r') as arq_entrada, open(saida, 'w') as arq_saida:
        for linha in arq_entrada:
            contador += 1
            AGG3D = Prog_Mod.apuraAGG3D(linha)
            arq_saida.write(linha)
    print("Gravou", contador, "sequencias")


def apuraclasses(entrada, saida):
    janela = 7
    contador = 0
    conteudo1 = ""

    with open(entrada, 'r') as arq_entrada, open(saida, 'w') as arq_saida:
        texto = arq_entrada.readlines()
        for linha in texto:
            contador += 1
            print("Seq. =", contador, "Tamanho =", len(linha))
            tabela2 = Prog_Mod.apura_agregacao(linha, " ")
            for I1 in range(0, len(tabela2) - janela):
                for I2 in range(I1+1, I1 + janela+1):
                    codigo = linha[I2]
                    texto_desc = Prog_Mod.descritores(codigo, " ")
                    conteudo1 += texto_desc
                if tabela2[I1 + 3] == "S":
                    conteudo1 += "Sim\r\n"
                else:
                    conteudo1 += "Nao\r\n"
                if (conteudo1[1] != " "):
                    arq_saida.write(conteudo1)
                    conteudo1 = ""
        print("Gravou", contador, "sequencias")


def movimenta(proteinArq, pdbLines, A3DLines, chain, CodigoPDB, chamador):
    resultado = []
    grava = " "
    matriz = []
    residuos = []
    atomos = []
    resNum = []
    cadeia = []
    A3D = []
    atoms = []
    coordenadaCA = []
    residuoCA = []
    resnumberCA = []
    desliga = 0
    contpdb = conta3d = 0
    opcao = 2

    for line in pdbLines:
        if len(line) < 1:
            continue
        if (line[0:4] == "ENDM"):
            desliga = 1
        elif (line[0:4] == "ATOM") and (desliga == 0):
            if (line[13:16] == "CA ") and (line[16:17] == " " or line[16:17] == "A"):
                contpdb += 1
                newAtom = Atom()
                newAtom.resID = line[8:3]
                newAtom.atomName = line[13:16]
                newAtom.resName = line[17:20]
                newAtom.resNumber = line[23:26]
                newAtom.coordinates = [line[31:38], line[40:46], line[48:54]]
                coordenadas = (float(newAtom.coordinates[0]), float(newAtom.coordinates[1]), float(newAtom.coordinates[2]))
                matriz.append(coordenadas)
                residuos.append(newAtom.resName)
                cadeia.append(line[21])
                resNum.append(newAtom.resNumber)
                atomos.append(newAtom.atomName)
                newAtom.bFactor = line[54:59]
                newAtom.tag = line[60:65]
                atoms.append(newAtom)

    if chamador == 0:
        for line1 in A3DLines:
            words1 = line1.split(",")
            if len(words1) >= 5 and words1[0] == "folded":
                A3D.append(float(words1[4]))
                conta3d += 1

    for ida in range(len(matriz)):
        coordenadaCA.append(matriz[ida])
        residuoCA.append(residuos[ida])
        resnumberCA.append(resNum[ida])

    grava = BuscaEsfera(proteinArq, coordenadaCA, residuoCA, resnumberCA, atomos, A3D, cadeia, chamador)

    if chamador == 0:
        contatohpFilename = os.path.join(SAIDA_DIR, f"{CodigoPDB[0:4]}_{chain}.csv")
        with open(contatohpFilename, "w") as contatohpFile:
            for line in grava:
                print(line, file=contatohpFile)
        return grava
    else:
        print("\nCadeia", chain)
        fez = 0
        for i in range(len(grava)):
            tamanho = len(grava) / 4
            prob = float(Prog_Predicao.checa_agregacao(grava[i], tamanho))
            hidrof = Prog_Mod.hidrofobico(residuos[i])
            grava1 = chain + ";" + str(int(resnumberCA[i])) + ';' + Prog_Mod.conversao(residuos[i], opcao) + ";" + str(round(prob, 3))
            grava1 += ";" + hidrof
            if grava1[0] == cadeia[i]:
                fez = 1
                resultado.append(grava1)
            elif fez == 1:
                break
        return resultado


def BuscaEsfera(ProteinArq, Posicao, residuos, resnumber, atomos, A3DScore, chain, chamador):
    completo = ""
    Estrutura = os.path.join(PDB_DIR if chamador == 1 else ARQUIVOS_DIR, ProteinArq)
    config_path = os.path.join(ESPECIAIS_DIR, "naccess.config")
    classifier = freesasa.Classifier(config_path)
    structure = freesasa.Structure(Estrutura, classifier)
    result = freesasa.calc(structure, freesasa.Parameters({'algorithm': freesasa.LeeRichards, 'n-slices': 100}))

    matrizR = [freesasa.selectArea(("TOT,resi 1", f"R1,chain {chain[ida]} and resi {resnumber[ida]}"), structure, result)["R1"] for ida in range(len(Posicao))]

    contatoHP = []
    opcao = 1
    agg = 0

    for ida, _ in enumerate(Posicao):
        if chamador != 1:
            Agrega1 = Prog_Mod.Agrega(A3DScore[ida])
        restrad = Prog_Mod.conversao(residuos[ida], opcao)
        if restrad:
            Total = sum([matrizR[j] for j in range(len(Posicao)) if 0 < Prog_Mod.dist(Posicao[ida], Posicao[j]) <= 10])
            AA = matrizR[ida]
            RSA = Prog_Mod.discrRSA((AA / Total) * 100 if Total > 0 else 0)

            gravar = ""
            if chamador != 1:
                gravar = restrad + Prog_Mod.descritores(residuos[ida], completo, chamador) + ";" + Agrega1 + ";0;00;00"
            else:
                gravar = ";" + restrad + Prog_Mod.descritores(residuos[ida], completo, chamador) + "00000"

            jda1 = 0
            for jda, _ in enumerate(Posicao):
                dist = Prog_Mod.dist(Posicao[ida], Posicao[jda])
                if 0 < dist <= 10:
                    restrad1 = Prog_Mod.conversao(residuos[jda], opcao)
                    RSA_j = Prog_Mod.discrRSA((matrizR[jda] / Total) * 100 if Total > 0 else 0)
                    if RSA_j or ida == jda:
                        jda1 += 1
                        if jda1 <= 21:
                            dist_disc = Prog_Mod.discrDist(round(dist, 2))
                            sinal = ";" if chamador != 1 else ""
                            gravar += ";" + "" + Prog_Mod.descritores(residuos[jda], completo, chamador) + sinal + str(round(dist, 2)) + sinal + RSA_j + sinal + Prog_Mod.apuraAGG(residuos[jda], agg)
            if chamador != 1:
                gravar += "; " + str(jda1)
            contatoHP.append(gravar)

    return contatoHP
