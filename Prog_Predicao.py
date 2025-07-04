# -*- coding: utf-8 -*-
import os
import joblib
import matplotlib.pyplot as plt
import numpy as np

# --- Configurações de arquivos ---
MODEL_DIR = 'arquivos/modelos'
MODEL_FILENAME = 'Mod_Oficial.sav'
MODEL_PATH = os.path.join(MODEL_DIR, MODEL_FILENAME)

TESTES_DIR = 'arquivos/testes'

# Se preferir garantir caminho absoluto usando BASE_DIR:
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH_ABS = os.path.join(BASE_DIR, MODEL_PATH)
TESTES_DIR_ABS = os.path.join(BASE_DIR, TESTES_DIR)

def checa_agregacao(descritores, tamanho):
    if not os.path.isfile(MODEL_PATH_ABS):
        raise FileNotFoundError(f"Arquivo do modelo não encontrado: {MODEL_PATH_ABS}")

    tabela = descritores.split(";")
    m = 20
    n = 1
    tabela1 = [[0] * m for _ in range(n)]

    for i in range(m):
        tabela1[0][i] = '0'

    cont = 0
    for i in range(len(tabela)):
        if i < 21 and i > 0:
            if tabela[i] != "":
                tabela1[0][i-1] = str(tabela[i])
                cont += 1
            else:
                tabela1[0][i-1] = '0'

    modelo = joblib.load(MODEL_PATH_ABS)
    A = modelo.predict(tabela1)

    if (A > 0.5 and cont < 8):
        A = A - 0.2

    A = A + 0.01

    if (A > 1):
        A = 1

    return A


def Plota_Resultado(resultado, CodigoPDB, chain):
    n = 6
    m = len(resultado)

    tabela1 = [[0] * m for _ in range(n)]
    tabela = tabelaa = ""
    i1 = 0
    for i in range(m):
        tabela = resultado[i].split(";")
        tabela1[3][i] = 0.5
        tabela1[0][i] = int(tabela[1])
        tabela1[1][i] = float(tabela[3])
        tabela1[2][i] = int(tabela[1])
        if i > 36:
            if i < 100:
                i1 += 1
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
    plt.xlim(int(tabela1[0][0]), int(tabela1[0][m-1]))
    plt.ylim(0, 1)

    plt.yticks(np.arange(0.0, 1.1, 0.1))

    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6

    caminho_figura = os.path.join(TESTES_DIR_ABS, f"{CodigoPDB}_{chain}.png")
    plt.title(f"Aggregation Propensity: {CodigoPDB}/{chain}", fontweight='bold')
    plt.savefig(caminho_figura, format='png')
    plt.figure(figsize=(10, 10))
    plt.clf()

    return
