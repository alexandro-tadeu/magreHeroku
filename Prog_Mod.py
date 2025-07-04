# -*- coding: utf-8 -*-
import math
import os
from atom import Atom  # Importa a classe Atom de um módulo chamado "atom"

def apura3D(proteinArq, pdbLines, A3DLines, chain, CodigoPDB):
    # Função principal para analisar estruturas 3D de proteínas
    arquivo = os.getcwd()  # Obtém o diretório de trabalho atual
    resultado = []  # Lista para armazenar resultados
    arquivo1 = arquivo[0:17] + "/arquivos/pdb/"  # Gera um caminho de arquivo usando o diretório atual
    grava = ""
    grava1 = ""
    matriz = []  # Lista para armazenar coordenadas
    residuos = []  # Lista para armazenar nomes dos resíduos
    atomos = []  # Lista para armazenar nomes dos átomos
    resNum = []  # Lista para armazenar números dos resíduos
    A3D = []  # Lista para armazenar linhas do arquivo A3D
    atoms = []  # Lista para armazenar instâncias da classe Atom
    coordenadaCA = []  # Lista para armazenar coordenadas do átomo CA
    residuoCA = []  # Lista para armazenar nomes dos resíduos do átomo CA
    resnumberCA = []  # Lista para armazenar números dos resíduos do átomo CA
    analiseCA = []  # Lista para armazenar resultados da análise

    ind_i = []  # Lista para armazenar índices de início de sequências de peptídeos
    ind_j = []  # Lista para armazenar índices de fim de sequências de peptídeos
    desliga = 0  # Variável de controle para desligar a leitura
    contpdb = conta3d = 0  # Contadores
    opcao = 2  # Opção
    teve = 0  # Flag
    tabela = ""  # String para armazenar informações da tabela
    tabres = []  # Lista para armazenar linhas do arquivo "Residuos.csv"

    arquivo3 = arquivo[0:17] + "/arquivos/especiais/"
    entrad1 = arquivo3 + "Residuos.csv"

    # Lê as linhas do arquivo "Residuos.csv"
    with open(entrad1) as arq_entrada:
        resLines = arq_entrada.readlines()

    for line in resLines:
        tabres.append(line)

    # Loop através das linhas do arquivo PDB
    for line in pdbLines:
        words = line
        if len(words) < 1:
            pass
        else:
            if (words[0:4] == "ENDM"):
                desliga = 1
            else:
                if (words[0:4] == "ATOM") and (desliga == 0):
                    if (words[21] == chain):
                        if (words[13:16] == "CA "):
                            if (words[16:17] == " " or words[16:17] == "A"):
                                contpdb = contpdb + 1
                                newAtom = Atom()
                                newAtom.resID = words[8:13]
                                newAtom.atomName = words[13:16]
                                newAtom.resName = words[17:20]
                                newAtom.resNumber = words[23:26]
                                # Foi feita uma mudança em relação a leitura de cada coluna
                                newAtom.coordinates = [
                                    words[31:38], words[40:46], words[48:54]]
                                coordenadas = (float(newAtom.coordinates[0]), float(
                                    newAtom.coordinates[1]), float(newAtom.coordinates[2]))
                                # Carrega uma matriz com as coordenadas
                                matriz.append(coordenadas)
                                residuos.append(newAtom.resName)
                                resNum.append(newAtom.resNumber)
                                atomos.append(newAtom.atomName)
                                newAtom.bFactor = words[54:59]
                                newAtom.tag = words[60:65]
                                atoms.append(newAtom)

    # Preenche listas com informações relevantes dos átomos CA
    for ida, tupla_a in enumerate(matriz):
        resnumberCA.append(resNum[ida])
        restrad = conversao(residuos[ida], opcao)
        residuoCA.append(restrad)

    i1 = 0
    achou = 'N'
    # Identifica sequências de peptídeos
    for i in range(0, len(residuoCA), 1):
        tabela = ""
        for j in range(i, len(residuoCA), 1):
            tabela += residuoCA[j]
            if (j-i) >= 3:
                if peptideo(tabela, tabres) == "S":
                    ind_i.append(i)
                    ind_j.append(j)
                    i1 = i1 + 1
                    achou = 'S'
    if achou == 'N':
        tabela = ""
        for j in range(0, len(residuoCA), 1):
            tabela += residuoCA[j]
        if verresiduo(tabela, tabres) == "S":
            ind_i.append(0)
            ind_j.append(len(tabela)-1)
            i1 = i1 + 1
            achou = 'S'

    # Gera uma lista de resultados da análise
    for i in range(0, len(residuoCA), 1):
        if (i1 > 0 and i >= ind_i[0] and i <= ind_j[0]):
            analiseCA.append("1.0")
        else:
            if (i1 > 1 and i >= ind_i[1] and i <= ind_j[1]):
                analiseCA.append("1.0")
            else:
                if (i1 > 2 and i >= ind_i[2] and i <= ind_j[2]):
                    analiseCA.append("1.0")
                else:
                    if (i1 > 3 and i >= ind_i[3] and i <= ind_j[3]):
                        analiseCA.append("1.0")
                    else:
                        if (i1 > 4 and i >= ind_i[4] and i <= ind_j[4]):
                            analiseCA.append("1.0")
                        else:
                            if (i1 > 4 and i >= ind_i[4] and i <= ind_j[4]):
                                analiseCA.append("1.0")
                            else:
                                if (i1 > 5 and i >= ind_i[5] and i <= ind_j[5]):
                                    analiseCA.append("1.0")
                                else:
                                    if (i1 > 6 and i >= ind_i[6] and i <= ind_j[6]):
                                        analiseCA.append("1.0")
                                    else:
                                        if (i1 > 7 and i >= ind_i[7] and i <= ind_j[7]):
                                            analiseCA.append("1.0")
                                        else:
                                            analiseCA.append("0.0")

    print(ind_i)
    print(ind_j)

    # Itera sobre os elementos de residuoCA para gerar entradas para a lista resultado
    for i in range(0, len(residuoCA), 1):
        # Adiciona informações sobre a análise à string grava
        grava += "folded,"
        grava += chain + "," + residuoCA[i] + "," + resnumberCA[i] + "," + str(analiseCA[i])
        
        # Adiciona a string grava à lista resultado
        resultado.append(grava)
        
        # Zera a variável grava para a próxima iteração
        grava = ""

    # Retorna a lista resultado, que contém informações sobre a análise da estrutura 3D da proteína
    return resultado


def peptideo(campo, tabela):
    """
    Verifica se o campo está presente na tabela.
    
    Parâmetros:
    - campo (str): O campo a ser verificado.
    - tabela (list): A lista de strings a ser pesquisada.
    
    Retorna:
    - completo (str): Retorna 'S' se o campo está presente na tabela, caso contrário, retorna uma string vazia.
    """
    completo = ""  # Variável de controle
    
    campo1 = campo + "\n"  # Adiciona uma quebra de linha ao campo
    for I in range(1, len(tabela), 1):
        if tabela[I] == campo1:
            completo = "S"  # Se o campo estiver presente na tabela, atribui 'S' à variável completo
    
    return completo


def verresiduo(campo, tabela):
    """
    Verifica se o campo está presente na tabela utilizando a função find().
    
    Parâmetros:
    - campo (str): O campo a ser verificado.
    - tabela (list): A lista de strings a ser pesquisada.
    
    Retorna:
    - completo (str): Retorna 'S' se o campo está presente em algum lugar da tabela, caso contrário, retorna uma string vazia.
    """
    completo = ""  # Variável de controle
    
    campo1 = campo + "\n"  # Adiciona uma quebra de linha ao campo
    for I in range(0, len(tabela), 1):
        if tabela[I].find(campo) != -1:
            completo = "S"  # Se o campo estiver presente em algum lugar na tabela, atribui 'S' à variável completo
    
    return completo


def Compara(resultado, arquivo2, chain):
    """
    Compara os resultados obtidos com um arquivo de referência.

    Parâmetros:
    - resultado (list): Lista contendo os resultados a serem comparados.
    - arquivo2 (str): Caminho do arquivo de referência.
    - chain (str): Identificador da cadeia a ser considerada.

    Retorna:
    - TabelaA (list): Lista contendo resultados da comparação (-1, 0, 1).
    """

    TabelaA = []  # Lista para armazenar resultados da comparação
    ind = 0  # Variável de controle para índice
    VP = FN = VN = FP = 0  # Inicialização de contadores

    # Lê as linhas do arquivo de referência
    with open(arquivo2) as a3dFile:
        a3dLines = a3dFile.readlines()

    # Loop através das linhas do arquivo de referência
    for line in a3dLines:
        words = line
        words = words.split(",")  # Divide a linha em palavras usando ',' como delimitador
        if len(words) < 1:
            pass
        elif (words[0] == "folded" and words[1] == chain):
            ind = ind + 1
            if (ind <= len(resultado)):
                res = resultado[ind-1].split(";")  # Divide a linha de resultado em palavras usando ';' como delimitador

                # Comparação dos resultados
                if (float(res[3]) > 0.5):
                    if (float(words[4]) > 0):
                        VP = VP + 1
                        TabelaA.append(0)
                    else:
                        TabelaA.append(1)
                        FP = FP + 1
                else:
                    if (float(words[4]) > 0):
                        TabelaA.append(-1)
                        FN = FN + 1
                    else:
                        TabelaA.append(0)
                        VN = VN + 1

    # Cálculo de métricas
    acuracia = (VP + VN) / (VP + VN + FP + FN)

    if (VP + FN) > 0:
        rec = VP / (VP+FN)
    if (VP + FP) > 0:
        prec = VP / (VP + FP)

    # Impressão de métricas
    print("VP=", VP, "FN=", FN, "VN=", VN, "FP=", FP)
    print(" ")
    print("Acuracia =", acuracia)

    if (VP + FN) > 0:
        print("Recall = ", rec)
    if (VP + FP) > 0:
        print("Precisao = ", prec)
    print(" ")

    return TabelaA


def monta_label(numero):
    """
    Monta uma lista de rótulos (labels) com base em um número fornecido.

    Parâmetros:
    - numero (int): Número que determina o conjunto específico de rótulos.

    Retorna:
    - labels (list): Lista de rótulos gerados.
    """

    labels = []  # Lista para armazenar os rótulos

    # Loop para gerar rótulos com base no número fornecido
    for I in range(1, 40):
        labels.append("D" + str(I))
        labels.append("RSA" + str(I))

        if numero == 1:
            labels.append("H" + str(I))
        if numero == 2:
            labels.append("H" + str(I))
            labels.append("T" + str(I))
        if numero == 3:
            labels.append("H" + str(I))
            labels.append("T" + str(I))
            labels.append("C" + str(I))
        if numero == 4:
            labels.append("H" + str(I))
            labels.append("C" + str(I))
        if numero == 5:
            labels.append("T" + str(I))
            labels.append("C" + str(I))
        if numero == 6:
            labels.append("T" + str(I))
        if numero == 7:
            labels.append("C" + str(I))

        if I == 1:
            print(labels)

    return labels


def hidrofobico(argument):
    """
    Mapeia aminoácidos para categorias hidrofóbicas.

    Parâmetros:
    - argument (str): Código de três letras do aminoácido.

    Retorna:
    - str: Categoria hidrofóbica correspondente.
    """

    switcher = {
        'ALA': 'H',
        'ARG': 'P',
        'ASN': 'P',
        'ASP': 'P',
        'CYS': 'H',
        'GLN': 'N',
        'GLU': 'P',
        'GLY': 'N',
        'HIS': 'N',
        'ILE': 'H',
        'LEU': 'H',
        'LYS': 'P',
        'MET': 'H',
        'PHE': 'H',
        'PRO': 'P',
        'SER': 'N',
        'THR': 'N',
        'TRP': 'H',
        'TYR': 'H',
        'VAL': 'H',
        'HSD': 'N',
        'HSE': 'N',
        'HSP': 'N',
    }

    # Obtém a categoria hidrofóbica correspondente ao aminoácido
    return switcher.get(argument, "nothing")


# converte o codigo da proteina em letras ou numeros
def conversao(proteina, opcao):
    """
    Converte códigos de aminoácidos entre diferentes representações.

    Parâmetros:
    - proteina (str): Código de três letras do aminoácido.
    - opcao (int): Opção de conversão (1 para número, 2 para letra).

    Retorna:
    - str: Resultado da conversão (número ou letra).
    """

    completo = ""
    tabela = [["ALA", "A", "01"],
              ["ARG", "R", "02"],
              ["ASN", "N", "03"],
              ["ASP", "D", "04"],
              ["CYS", "C", "05"],
              ["GLN", "Q", "06"],
              ["GLU", "E", "07"],
              ["GLY", "G", "08"],
              ["HIS", "H", "09"],
              ["ILE", "I", "10"],
              ["LEU", "L", "11"],
              ["LYS", "K", "12"],
              ["MET", "M", "13"],
              ["PHE", "F", "14"],
              ["PRO", "P", "15"],
              ["SER", "S", "16"],
              ["THR", "T", "17"],
              ["TRP", "W", "18"],
              ["TYR", "Y", "19"],
              ["VAL", "V", "20"]]

    for I in range(0, 20):
        if tabela[I][0] == proteina:
            if opcao == 1:
                completo = tabela[I][2]
            else:
                completo = tabela[I][1]

    if completo == "":
        completo = "99"

    return completo
def conversao1(proteina, opcao):
    """
    Converte códigos de aminoácidos entre diferentes representações.

    Parâmetros:
    - proteina (str): Código de três letras do aminoácido.
    - opcao (int): Opção de conversão (1 para código de três letras, 2 para letra única).

    Retorna:
    - str: Resultado da conversão (código de três letras ou letra única).
    """

    completo = ""
    tabela = [["ALA", "A", "01"],
              ["ARG", "R", "02"],
              ["ASN", "N", "03"],
              ["ASP", "D", "04"],
              ["CYS", "C", "05"],
              ["GLN", "Q", "06"],
              ["GLU", "E", "07"],
              ["GLY", "G", "08"],
              ["HIS", "H", "09"],
              ["ILE", "I", "10"],
              ["LEU", "L", "11"],
              ["LYS", "K", "12"],
              ["MET", "M", "13"],
              ["PHE", "F", "14"],
              ["PRO", "P", "15"],
              ["SER", "S", "16"],
              ["THR", "T", "17"],
              ["TRP", "W", "18"],
              ["TYR", "Y", "19"],
              ["VAL", "V", "20"]]

    for I in range(0, 20):
        if tabela[I][0] == proteina:
            if opcao == 1:
                completo = tabela[I][0]
            else:
                completo = tabela[I][1]

    return completo

def conversao2(proteina, opcao):
    """
    Converte códigos de aminoácidos entre diferentes representações.

    Parâmetros:
    - proteina (str): Código de uma letra do aminoácido.
    - opcao (int): Opção de conversão (1 para código de três letras, 2 para código de número).

    Retorna:
    - str: Resultado da conversão (código de três letras ou número).
    """

    completo = ""
    tabela = [["ALA", "A", "01"],
              ["ARG", "R", "02"],
              ["ASN", "N", "03"],
              ["ASP", "D", "04"],
              ["CYS", "C", "05"],
              ["GLN", "Q", "06"],
              ["GLU", "E", "07"],
              ["GLY", "G", "08"],
              ["HIS", "H", "09"],
              ["ILE", "I", "10"],
              ["LEU", "L", "11"],
              ["LYS", "K", "12"],
              ["MET", "M", "13"],
              ["PHE", "F", "14"],
              ["PRO", "P", "15"],
              ["SER", "S", "16"],
              ["THR", "T", "17"],
              ["TRP", "W", "18"],
              ["TYR", "Y", "19"],
              ["VAL", "V", "20"]]

    for I in range(0, 20):
        if tabela[I][1] == proteina:
            if opcao == 1:
                completo = tabela[I][0]
            else:
                completo = tabela[I][2]

    return completo

def conversao3(proteina, opcao):
    """
    Converte códigos de aminoácidos entre diferentes representações.

    Parâmetros:
    - proteina (str): Número associado ao aminoácido.
    - opcao (int): Opção de conversão (1 para código de uma letra, 2 para código de três letras).

    Retorna:
    - str: Resultado da conversão (código de uma letra ou três letras do aminoácido).
    """

    completo = ""
    tabela = [["ALA", "A", "01"],
              ["ARG", "R", "02"],
              ["ASN", "N", "03"],
              ["ASP", "D", "04"],
              ["CYS", "C", "05"],
              ["GLN", "Q", "06"],
              ["GLU", "E", "07"],
              ["GLY", "G", "08"],
              ["HIS", "H", "09"],
              ["ILE", "I", "10"],
              ["LEU", "L", "11"],
              ["LYS", "K", "12"],
              ["MET", "M", "13"],
              ["PHE", "F", "14"],
              ["PRO", "P", "15"],
              ["SER", "S", "16"],
              ["THR", "T", "17"],
              ["TRP", "W", "18"],
              ["TYR", "Y", "19"],
              ["VAL", "V", "20"]]

    for I in range(0, 20):
        if tabela[I][2] == proteina:
            if opcao == 1:
                completo = tabela[I][0]
            else:
                completo = tabela[I][1]

    return completo

# Apura o AGG (tabela do Aggrescan)

def apuraAGG(codigo, agg):
    """
    Calcula a propensão para formação de agregados para um aminoácido.

    Parâmetros:
    - codigo (str): Código do aminoácido (uma letra).
    - agg (float): Valor inicial da propensão (inicializado como 0.0, mas não utilizado na função).

    Retorna:
    - float: Valor da propensão para formação de agregados do aminoácido.
    """

    agg = 0.0

    propensao = [["A", "ALA", "1.822", "12"],
                 ["R", "ARG", "1.754", "04"],
                 ["N", "ASN", "1.594", "03"],
                 ["D", "ASP", "1.380", "01"],
                 ["C", "CYS", "1.159", "13"],
                 ["Q", "GLN", "1.037", "05"],
                 ["E", "GLU", "0.910", "02"],
                 ["G", "GLY", "0.604", "08"],
                 ["H", "HIS", "-0.036", "06"],
                 ["I", "ILE", "-0.159", "20"],
                 ["L", "LEU", "-0.294", "16"],
                 ["P", "LYS", "-0.334", "07"],
                 ["M", "MET", "-0.535", "14"],
                 ["F", "PHE", "-0.931", "19"],
                 ["P", "PRO", "-1.033", "09"],
                 ["S", "SER", "-1.231", "10"],
                 ["T", "THR", "-1.240", "11"],
                 ["W", "TRP", "-1.302", "15"],
                 ["Y", "TYR", "-1.412", "16"],
                 ["V", "VAL", "-1.836", "18"]]

    for I in range(0, 20):
        if propensao[I][1] == codigo:
            agg = propensao[I][3]

    return agg

# Avalia o score do AGRESCAN 3d (arquivo)
def Agrega(Score):
    """
    Determina se uma pontuação indica a formação de agregados.

    Parâmetros:
    - Score (float): Pontuação a ser avaliada.

    Retorna:
    - str: '1' se a pontuação for maior que 0.0, '0' caso contrário.
    """

    agrega = " "

    if float(Score) > 0.0:
        agrega = '1'
    else:
        # if float(Score) == 0.0:
        #  agrega = '1'
        # else:
        agrega = '0'

    return agrega

# Discretiza o RSA

def discrRSA(AA):
    """
    Converte um valor de área acessível (AA) em uma categoria discreta de RSA.

    Parâmetros:
    - AA (float): Valor da área acessível.

    Retorna:
    - str: Categoria discreta de RSA conforme as faixas estabelecidas.
    """

    RSA = 0

    if AA < 1:
        RSA = "00"
    elif AA < 5:
        RSA = "01"
    elif AA < 10:
        RSA = "05"
    elif AA < 20:
        RSA = "10"
    elif AA < 30:
        RSA = "20"
    elif AA < 40:
        RSA = "30"
    elif AA < 50:
        RSA = "40"
    elif AA < 80:
        RSA = "60"
    elif AA < 100:
        RSA = "80"
    else:
        RSA = "99"

    return RSA

#   Discretiza a distancia
def discrDist(Dist):
    """
    Converte um valor de distância (Dist) em uma categoria discreta.

    Parâmetros:
    - Dist (float): Valor da distância.

    Retorna:
    - int: Categoria discreta de distância.
    """

    D = 0

    if Dist < 1:
        D = 0
    elif Dist < 2:
        D = 1
    elif Dist < 3:
        D = 2
    elif Dist < 4:
        D = 3
    elif Dist < 5:
        D = 4
    elif Dist < 6:
        D = 5
    elif Dist < 7:
        D = 6
    elif Dist < 8:
        D = 7
    elif Dist < 9:
        D = 8
    elif Dist < 10:
        D = 9
    elif Dist < 11:
        D = 10
    elif Dist < 12:
        D = 11
    elif Dist < 13:
        D = 12
    elif Dist < 14:
        D = 13
    elif Dist < 15:
        D = 14
    else:
        D = 15

    return D

def dist(res1, res2):
    """
    Calcula a distância Euclidiana tridimensional entre dois pontos.

    Parâmetros:
    - res1 (tuple): Coordenadas tridimensionais do primeiro ponto.
    - res2 (tuple): Coordenadas tridimensionais do segundo ponto.

    Retorna:
    - float: Distância Euclidiana entre os dois pontos.
    """

    dis = math.sqrt(((res1[0] - res2[0]) ** 2) +
                    ((res1[1] - res2[1]) ** 2) + ((res1[2] - res2[2]) ** 2))

    return dis

def descritores(codigo, completo, chamador):
    """
    Gera descritores com base em um código de aminoácido.

    Parâmetros:
    - codigo (str): Código do aminoácido.
    - completo (str): String que será preenchida com os descritores.
    - chamador (int): Parâmetro indicando quem chamou a função.

    Retorna:
    - str: String completa com os descritores.
    """

    tabela = [
        ["A", "ALA", "02", "1.410", "0.750", "1", "2"],
        ["R", "ARG", "01", "1.210", "0.850", "3", "3"],
        ["N", "ASN", "01", "0.730", "0.630", "2", "2"],
        ["D", "ASP", "01", "0.820", "0.550", "2", "1"],
        ["C", "CYS", "03", "0.850", "1.360", "1", "2"],
        ["Q", "GLN", "01", "1.260", "0.720", "3", "2"],
        ["E", "GLU", "01", "1.390", "0.650", "3", "1"],
        ["G", "GLY", "02", "0.440", "0.670", "1", "2"],
        ["H", "HIS", "01", "0.870", "0.990", "3", "3"],
        ["I", "ILE", "03", "1.040", "1.790", "3", "2"],
        ["L", "LEU", "03", "1.280", "1.150", "3", "2"],
        ["K", "LYS", "01", "1.170", "0.760", "3", "3"],
        ["M", "MET", "03", "1.260", "1.010", "3", "2"],
        ["F", "PHE", "03", "1.000", "1.400", "3", "2"],
        ["P", "PRO", "02", "0.440", "0.400", "2", "2"],
        ["S", "SER", "02", "0.760", "0.810", "1", "2"],
        ["T", "THR", "02", "0.780", "1.210", "2", "2"],
        ["W", "TRP", "03", "1.070", "1.230", "3", "2"],
        ["Y", "TYR", "02", "1.980", "1.370", "3", "2"],
        ["V", "VAL", "03", "0.910", "2.000", "2", "2"],
        ["X", "GLY", "02", "0.440", "0.670", "1", "2"]
    ]

    for I in range(0, 20):
        if tabela[I][1] == codigo:
            if chamador != 1:
                completo += ";" + tabela[I][2] + ";"
                completo += tabela[I][5] + ";"
                completo += tabela[I][6]
            else:
                completo += tabela[I][2]
                completo += tabela[I][5]
                completo += tabela[I][6]

    return completo

#   apura agregacao no modelo de treinamento
