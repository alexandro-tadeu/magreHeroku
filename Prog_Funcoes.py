# -*- coding: utf-8 -*-
import os
import Prog_Mod
from atom import Atom

import matplotlib.pyplot as plt
import math
import freesasa
import numpy as np
import pandas as pd


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


def calculascore(entrada, saida):

    janela = 7
    contador = 0
    contsim = 0
    contnao = 0
    codigo = " "
    conteudo1 = ""
    completo = " "
    tabela = " "
    tabela1 = " "
    tabela2 = " "
    AGG3D = 0
    agg3d = 0
    agg = 0
    tabela3 = []
    tabagg = " "
    arq_entrada = open(entrada, 'r')
    arq_saida = open(saida, 'w')

    texto = arq_entrada.readlines()

    for linha in texto:
        contador = contador + 1
        tabagg = linha
        # print(linha)
        AGG3D = Prog_Mod.apuraAGG3D(tabagg)
       # if (AGG3D > 0.0):
      #     print(contador,AGG3D)
        # if contador <= 2:

        arq_saida.write(linha)

    print("Gravou ", contador, " sequencias")
    arq_entrada.close()
    arq_saida.close()
    return


def apuraclasses(entrada, saida):
    janela = 7
    contador = 0
    contsim = 0
    contnao = 0
    codigo = " "
    conteudo1 = ""
    completo = " "
    tabela = " "
    tabela1 = " "
    tabela2 = " "
    arq_entrada = open(entrada, 'r')
    arq_saida = open(saida, 'w')

    texto = arq_entrada.readlines()

    for linha in texto:
        contador = contador + 1
        print("Seq. = ", contador, "Tamanho =  ", len(linha))
        tabela2 = Prog_Mod.apura_agregacao(linha, tabela)
        print(tabela2)
        #print("tamanho agregacao = ",len(tabela2))
        print(linha)
        for I1 in range(0, len(tabela2) - janela):
            for I2 in range(I1+1, I1 + janela+1):
                codigo = linha[I2]  # print (I1,I2)

                texto = Prog_Mod.descritores(codigo, completo)
                conteudo1 += texto
            # print(I1)
            if tabela2[I1 + 3] == "S":
                conteudo1 += "Sim\r\n"

                contsim = contsim + 1
            else:
                contnao = contnao + 1
                conteudo1 += "Nao\r\n"
                # print(conteudo1)
            if (conteudo1[1] != " "):
                arq_saida.write(conteudo1)
                conteudo1 = ""

        # arq_saida.write(linha)

    print("Gravou ", contador, " sequencias")
    arq_entrada.close()
    arq_saida.close()
    return


def movimenta(proteinArq, pdbLines, A3DLines, chain, CodigoPDB, chamador):
    print("aa")
    arquivo = os.getcwd()
    resultado = []
    arquivo1 = arquivo[0:17] + "/arquivos/pdb/"
    grava = " "
    grava1 = ""
    matriz = []  # matriz com as coordenadas
    residuos = []  # matriz
    atomos = []  # matriz com os atomos
    resNum = []
    A3D = []
    atoms = []
    coordenadaCA = []
    residuoCA = []
    resnumberCA = []
    desliga = 0
    contpdb = conta3d = 0
    opcao = 2
    teve = 0
    print("inicio")
    for line in pdbLines:
        words = line
        if len(words) < 1:
            pass
        else:
            # print(words)
            if (words[0:4] == "ENDM"):
                desliga = 1
            else:
                if (words[0:4] == "ATOM") and (desliga == 0):
                    if (words[21] == chain):
                        if (words[13:16] == "CA "):

                            if (words[16:17] == " " or words[16:17] == "A"):
                                # print("dauqui",words)
                                contpdb = contpdb+1
                                newAtom = Atom()
                                newAtom.resID = words[8:3]
                                newAtom.atomName = words[13:16]
                                newAtom.resName = words[17:20]
                                newAtom.resNumber = words[23:26]
                                # Foi feita uma mudança em relação a leitura de cada coluna
                                newAtom.coordinates = [
                                    words[31:38], words[40:46], words[48:54]]
                                coordenadas = (float(newAtom.coordinates[0]), float(
                                    newAtom.coordinates[1]), float(newAtom.coordinates[2]))
                                # carrega uma matriz com as coordenadas
                                matriz.append(coordenadas)
                                residuos.append(newAtom.resName)
                                resNum.append(newAtom.resNumber)
                                atomos.append(newAtom.atomName)
                                newAtom.bFactor = words[54:59]
                                newAtom.tag = words[60:65]
                                atoms.append(newAtom)
    # print("segundo")

    if chamador == 0:
        for line1 in A3DLines:
            words = line1
            # print(words)
            words1 = words.split(",")
            if len(words1) < 1:
                pass
            # elif (words[4:5] == chain):
            elif (words1[0] == "folded" and words1[1] == chain):
                # print(words1)
                #newAtom.A3D = float(words[26:32])
                newAtom.A3D = float(words1[4])
                conta3d = conta3d+1
                A3D.append((newAtom.A3D))

    # print(chain,contpdb,conta3d)
    #print("contadores= ",chain,contpdb,conta3d)
    for ida, tupla_a in enumerate(matriz):
        coordenadaCA.append(matriz[ida])
        residuoCA.append(residuos[ida])
        resnumberCA.append(resNum[ida])
    # print("res",residuos)
    if (contpdb == conta3d) or (chamador == 1):
        grava = BuscaEsfera(proteinArq, coordenadaCA, residuoCA,
                            resnumberCA, atomos, A3D, chain, chamador)

    if chamador == 0:
        # PastaArq=arquivo1
        PastaArq = arquivo[0:17] + "/arquivos/"
        contatohpFilename = PastaArq + "saida/" + \
            CodigoPDB[0:4] + "_" + chain + ".csv"
        with open(contatohpFilename, "w") as contatohpFile:
            for line in grava:
                print(line, file=contatohpFile)
        #print("GRAVACAO:" + " " + contatohpFilename)
        return grava
    else:
        print(" ")
        print("Cadeia", chain)
        n = 6
        m = len(grava)
        tabela = [[0] * m for i in range(n)]
        for i in range(0, len(grava), 1):
            # print("resultado",grava[i])
            prob = float(Prog_Mod.checa_agregacao(grava[i]))
            hidrof = Prog_Mod.hidrofobico(residuos[i])
            # print(AAA)
            grava1 = chain + ";" + str(int(resnumberCA[i])) + ';' + Prog_Mod.conversao(
                residuos[i], opcao) + ";" + str(round(prob, 2))
            grava1 += ";" + hidrof
            tabela[0][i] = int(resnumberCA[i])
            tabela[1][i] = float(prob)
            tabela[3][i] = 0.5
            tabela[2][i] = int(resnumberCA[i])

            tabela[5][i] = int(resnumberCA[i])
            resultado.append(grava1)
            teve = 1

        arquivo2 = arquivo[0:17] + "/arquivos/A3D_entrada/" + \
            CodigoPDB + "/" + "A3D.csv"

        TabelaA = []
        #print("vai comparar")
        # TabelaA.append(Prog_Mod.Compara(resultado,arquivo2,chain))
      # print(TabelaA)

        for i in range(0, len(grava), 1):
            tabela[4][i] = TabelaA[0][i]

        plt.plot(tabela[0], tabela[1])
        plt.plot(tabela[2], tabela[3])
        plt.xlabel('Resíduos(Pos)', fontweight='bold')
        plt.ylabel('Probabilidade', fontweight='bold')
        plt.title("Propensão a Agregação - PDB: " + CodigoPDB +
                  "/ Chain: " + chain, fontweight='bold')
        plt.savefig(arquivo[0:17] + "/arquivos/testes/" +
                    CodigoPDB + '_' + chain + '.png', format='png')
        plt.clf()

        plt.plot(tabela[5], tabela[4])
        plt.xlabel('Resíduos(Pos)', fontweight='bold')
        plt.ylabel('Erros', fontweight='bold')
        # plt.legend()
        if (teve == 1):
            # print(tabela[2][0])
            plt.axis([int(tabela[5][0]), int(
                tabela[5][0]) + len(grava), -1.2, 1.2])
        # plt.axis(str(int(resnumberCA[0])),len(grava), 0 , 1])
        plt.title("Demonstrativo de Erros - PDB: " + CodigoPDB +
                  "/ Chain: " + chain, fontweight='bold')
        plt.savefig(arquivo[0:17] + "/arquivos/testes/" +
                    CodigoPDB + '_' + chain + 'A.png', format='png')
        plt.clf()

        return resultado


#################################################################################
def BuscaEsfera(ProteinArq, Posicao, residuos, resnumber, atomos, A3DScore, chain, chamador):

    arquivo = os.getcwd()
    # print("BuscaEsfera")
    if chamador == 1:
        arquivo1 = arquivo[0:17] + "/arquivos/pdb/"
    else:
        arquivo1 = arquivo[0:17] + "/arquivos/"

    arquivo2 = arquivo[0:17] + "/arquivos/especiais/"
    # print(proteinarq)

    completo = ""
    Estrutura = arquivo1 + ProteinArq

    classifier = freesasa.Classifier(arquivo2 + "naccess.config")
    structure = freesasa.Structure(Estrutura, classifier)
    result = freesasa.calc(structure,
                           freesasa.Parameters({'algorithm': freesasa.LeeRichards,
                                                # 'probe-radius' : 10,
                                                'n-slices': 100
                                                }))

    AA = 0
    Dist = 0
    jda1 = 0
    key = ""
    opcao = 1
    Total = 0
    contatoHP = []
    agg = 0

    matrizR = []

    for ida, tupla_a in enumerate(Posicao):
        A = 'R1' + ',' + 'chain ' + chain + ' and resi ' + str(resnumber[ida])

        selections = freesasa.selectArea(('TOT,resi 1', A),
                                         structure, result)
        key = "R1"
        AA = selections[key]
        matrizR.append(selections[key])

    for ida, tupla_a in enumerate(Posicao):
        if chamador != 1:
            # print(ida)
            Agrega1 = Prog_Mod.Agrega(A3DScore[ida])
        restrad = Prog_Mod.conversao(residuos[ida], opcao)
        Total = 0
        gravar = " "
        for jdat, tupla_b in enumerate(Posicao):
            if (Prog_Mod.dist(Posicao[ida], Posicao[jdat]) >= 0 and Prog_Mod.dist(Posicao[ida], Posicao[jdat]) <= 10.00):
                AA = matrizR[jdat]
                Total = Total + AA
        AA = matrizR[ida]
        # print(round(AA,0))
        if (Total > 0):
            RSA = (AA/Total)*100
        else:
            RSA = 0
        # RSA=0 #especial
        RSA = Prog_Mod.discrRSA(RSA)
        # print(RSA)
        if chamador != 1:
            gravar = restrad + ";" + "0.0" + ";" + str(round(RSA, 0))
            gravar += ";" + Agrega1
        else:
            gravar += ";" + restrad + ";" + "0.0" + ";" + str(round(RSA, 0))
        gravar += Prog_Mod.descritores(residuos[ida], completo)
        # print(gravar)
        # print(gravar)
        # print("AA",residuos[ida],Prog_Mod.apuraAGG(residuos[ida],agg))
        gravar += ";" + str(float(Prog_Mod.apuraAGG(residuos[ida], agg)))
        # print(gravar)
        for jda, tupla_b in enumerate(Posicao):
            if (Prog_Mod.dist(Posicao[ida], Posicao[jda]) > 0 and Prog_Mod.dist(Posicao[ida], Posicao[jda]) <= 10.00):
                restrad1 = Prog_Mod.conversao(residuos[jda], opcao)
                AA = matrizR[jda]
                if (Total > 0):
                    RSA = (AA/Total)*100
                else:
                    RSA = 0
                if (RSA >= 4.5 or ida == jda):
                    RSA = Prog_Mod.discrRSA(RSA)
                    jda1 = jda1 + 1
                   # RSA=0 # especial
                    # print(str(round(RSA,0)
                    Dist = Prog_Mod.discrDist(
                        round(Prog_Mod.dist(Posicao[ida], Posicao[jda]), 2))
                    gravar += ";" + restrad1 + ";" + \
                        str(round(Dist, 2)) + ";" + str(round(RSA, 0))
                    # print(gravar)
                    gravar += Prog_Mod.descritores(residuos[jda], completo)
                    #gravar +=";" +  str(float(Prog_Mod.apuraAGG(residuos[jda],agg)))
                    gravar += ";" + "0.0"
                else:
                    pass
            else:
                pass
        #print("linha = ",gravar)
        if chamador != 1:
            gravar += "; " + str(jda1)  # so pra chegar o tamanho

        contatoHP.append(gravar)
        jda1 = 0

    return contatoHP
