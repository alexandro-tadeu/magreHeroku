# -*- coding: utf-8 -*-
import os
import Prog_Mod
from atom import Atom

import matplotlib.pyplot as plt
import freesasa
import numpy as np
import pandas as pd
import Prog_Predicao


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
   # print(resultado)
    resultado11 = ""
    resultado12 = ""
    resultado13 = ""
    resultados = ""
    resultadox = []

    for i in range(0, len(resultado)-1, 1):
        resultado11 = resultado[i-1].split(";")
        resultado12 = resultado[i].split(";")
        resultado13 = resultado[i+1].split(";")
        if float(resultado12[3]) > 0.5:
            if float(resultado11[3]) < 0.5 and float(resultado13[3]) < 0.5:
               # print(i,resultado12)
                resultados = resultado12[0]+";"
                resultados += resultado12[1]+";"
                resultados += resultado12[2]+";"
                resultados += str(float(resultado12[3])/2) + ";"
                resultados += resultado12[4]+";"
                resultadox.append(resultados)
            else:
                resultadox.append(resultado[i])
        else:
            resultadox.append(resultado[i])

    resultadox.append(resultado[i+1])

    return resultadox


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
        # print(tabela2)
        # print("tamanho agregacao = ",len(tabela2))
        # print(linha)
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

    arquivo = os.getcwd()
    resultado = []
    arquivo1 = arquivo[0:17] + "/arquivos/pdb/"
    grava = " "
    grava1 = ""
    matriz = []  # matriz com as coordenadas
    residuos = []  # matriz
    atomos = []  # matriz com os atomos
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
    teve = 0
    # print("inicio")
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
                    # if (words[21] == chain):
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
                            cadeia.append(words[21])
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
            # elif (words1[0] == "folded" and words1[1] == chain):
            elif (words1[0] == "folded"):
                # print(words1)
                # newAtom.A3D = float(words[26:32])
                newAtom.A3D = float(words1[4])
                conta3d = conta3d+1
                A3D.append((newAtom.A3D))

    # print(chain,contpdb,conta3d)
    # print("contadores= ",chain,contpdb,conta3d)
    for ida, tupla_a in enumerate(matriz):
        coordenadaCA.append(matriz[ida])
        residuoCA.append(residuos[ida])
        resnumberCA.append(resNum[ida])

    grava = BuscaEsfera(proteinArq, coordenadaCA, residuoCA,
                        resnumberCA, atomos, A3D, cadeia, chamador)

    if chamador == 0:
        PastaArq = arquivo[0:17] + "/arquivos/"
        contatohpFilename = PastaArq + "saida/" + \
            CodigoPDB[0:4] + "_" + chain + ".csv"
        with open(contatohpFilename, "w") as contatohpFile:
            for line in grava:
                print(line, file=contatohpFile)
        return grava
    else:
        print(" ")
        print("Cadeia", chain)

        n = 6
        m = len(grava)
        tabela = [[0] * m for i in range(n)]
        fez = 0
        for i in range(0, len(grava), 1):
            tamanho = len(grava)/4
            prob = float(Prog_Predicao.checa_agregacao(grava[i], tamanho))
            hidrof = Prog_Mod.hidrofobico(residuos[i])
            # print(AAA)
            grava1 = chain + ";" + str(int(resnumberCA[i])) + ';' + Prog_Mod.conversao(
                residuos[i], opcao) + ";" + str(round(prob, 3))
            grava1 += ";" + hidrof

            if (grava1[0] == cadeia[i]):
                fez = 1
                resultado.append(grava1)
            else:
                if fez == 1:
                    break

        return resultado

#################################################################################


def BuscaEsfera(ProteinArq, Posicao, residuos, resnumber, atomos, A3DScore, chain, chamador):
    # print(resnumber)
    arquivo = os.getcwd()
    # print("BuscaEsfera")
    # print(A3DScore)
    if chamador == 1:
        arquivo1 = arquivo[0:17] + "/arquivos/pdb/"
    else:
        arquivo1 = arquivo[0:17] + "/arquivos/"

    arquivo2 = arquivo[0:17] + "/arquivos/especiais/"
    # print(ProteinArq)

    completo = ""
    Estrutura = arquivo1 + ProteinArq
    # print(Estrutura)
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
        # print(ida,chain)
        A = 'R1' + ',' + 'chain ' + \
            chain[ida] + ' and resi ' + str(resnumber[ida])

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
        if restrad != "":
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
            gravar = " "
            if chamador != 1:
                gravar = restrad + Prog_Mod.descritores(
                    residuos[ida], completo, chamador) + ";" + Agrega1 + ";" + "0" + ";" + "00"
                gravar += ";" + "00"
            else:
                # gravar += ";" +  restrad + Prog_Mod.descritores(residuos[ida],completo) + ";" + "0" + ";" + "0" + ";" + "0"
                gravar = ";" + restrad + \
                    Prog_Mod.descritores(
                        residuos[ida], completo, chamador) + "00000"
                # gravar = ";" + "9" + restrad + Prog_Mod.descritores(residuos[ida],completo,chamador) +  "000"
                # ravar=";" + restrad
                # print("vai montar",gravar)
        #

            for jda, tupla_b in enumerate(Posicao):
                # if jda <= 50:
                if (Prog_Mod.dist(Posicao[ida], Posicao[jda]) > 0 and Prog_Mod.dist(Posicao[ida], Posicao[jda]) <= 10.00):
                    restrad1 = Prog_Mod.conversao(residuos[jda], opcao)
                    AA = matrizR[jda]
                    if chamador != 1:
                        sinal = ";"
                    else:
                        sinal = ""
                    if (Total > 0):
                        RSA = (AA/Total)*100
                    else:
                        RSA = 0
                    # if jda <= 50:
                    if (RSA > 0 or ida == jda):
                        RSA = Prog_Mod.discrRSA(RSA)
                        # print(RSA)
                        jda1 = jda1 + 1
                        if jda1 <= 21:
                            restrad1 = ""
                            Dist = Prog_Mod.discrDist(
                                round(Prog_Mod.dist(Posicao[ida], Posicao[jda]), 2))
                            # print(Prog_Mod.apuraAGG(residuos[jda],agg))
                        # gravar += ";" + restrad1 + ";" +  Prog_Mod.apuraAGG(residuos[jda],agg) + ";" + str(round(Dist,2)) + ";" + str(round(RSA,0))
                            gravar += ";" + restrad1 + Prog_Mod.descritores(residuos[jda], completo, chamador) + sinal + str(
                                round(Dist, 2)) + sinal + RSA + sinal + Prog_Mod.apuraAGG(residuos[jda], agg)
                            # gravar += Prog_Mod.descritores(residuos[jda],completo)
                            # gravar +=";"
                        else:
                            jda1 = 21
                    else:
                        pass
                    # else:
                        #    pass
                else:
                    pass
        # print("linha = ",gravar)
            if chamador != 1:
                gravar += "; " + str(jda1)  # so pra chegar o tamanho

        contatoHP.append(gravar)
        jda1 = 0

    return contatoHP
