# -*- coding: utf-8 -*-
import os
import joblib
import matplotlib.pyplot as plt
import numpy as np

def checa_agregacao(descritores, tamanho):
    arquivo = os.getcwd()
    modelo = "Mod_Oficial.sav"
    arquivo1 = arquivo[0:17] + "/ARQUIVOS/MODELOS/" + modelo

    tabela = descritores.split(";")
    m = 20
    n = 1
    tabela1 = [[0] * m for i in range(n)]
    tem = " "
    for i in range(0, m, 1):
        tabela1[0][i] = '0'
    cont = 0
    for i in range(0, len(tabela), 1):
        if i < 21:
            if i > 0:
                if tabela[i] != "":
                    tabela1[0][i-1] = str(tabela[i])
                    cont = cont + 1
                else:
                    tabela1[0][i-1] = '0'

    modelo = joblib.load(arquivo1)

    A = modelo.predict(tabela1)

    if (A > 0.5 and cont < 8):
        A = A - 0.2

    A = A + 0.01

    if (A > 1):
        A = 1

    return A


def Plota_Resultado(resultado, CodigoPDB, chain):
    arquivo = os.getcwd()

    n = 6
    m = len(resultado)

    tabela1 = [[0] * m for i in range(n)]
    tabela = tabelaa = ""
    i1 = 0
    for i in range(0, len(resultado), 1):
        tabela = resultado[i].split(";")
        tabela1[3][i] = 0.5
        tabela1[0][i] = int(tabela[1])
        tabela1[1][i] = float(tabela[3])
        tabela1[2][i] = int(tabela[1])
        if i > 36:
            if i < 100:
                i1 = i1 + 1
                tabela1[4][i] = int(tabelaa[1])
                tabela1[5][i] = float(tabelaa[3])
            else:
                tabela1[4][i] = 99
                tabela1[5][i] = 0.0
        else:
            tabela1[4][i] = 36
            tabela1[5][i] = 0.0

    plt.clf()

    plt.plot(tabela1[0], tabela1[1], color='blue')
    plt.plot(tabela1[0], tabela1[1], '.', color='Black')
    plt.plot(tabela1[2], tabela1[3], '--', color='green')

    plt.xlabel('Residues', fontweight='bold')
    plt.ylabel('Probability', fontweight='bold')
    plt.xlim(int(tabela1[0][0]), int(tabela1[0][i]))
    plt.ylim(0, 1)

    plt.yticks(np.arange(0.0, 1.1, 0.1))

    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6

    plt.title("Aggregation Propensity: " + CodigoPDB +
              "/" + chain, fontweight='bold')
    plt.savefig(arquivo[0:17] + "/ARQUIVOS/TESTES/" +
                CodigoPDB + '_' + chain + '.png', format='png')
    plt.figure(figsize=(10, 10))
    plt.clf()

    return
# MODULO PARA  CALCULAR O SCORE A3D
