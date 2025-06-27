import streamlit as st
import requests
import matplotlib.pyplot as plt
import os
import Prog_Funcoes1

# Variável de controle para verificar se o gráfico já foi gerado
graficos_gerados = set()

# Variável global para armazenar a cadeia
chain = " "

# Função para plotar o gráfico
def Plota_Resultado(resultado, codigo_pdb, chain):
    cor_linha = "black"
    cor_forte = "red"
    cor_fraca = "orange"
    limite_superior = 0.8
    limite_inferior = 0.2
    limiar = 0.5  # Altera de "Threshold" para "Limiar"

    tamanho = len(resultado)
    tabela1 = [[0] * tamanho for _ in range(2)]  # Inicializa tabela para resíduo e probabilidade

    for i in range(tamanho):
        tabela = resultado[i].split(";")
        tabela1[0][i] = int(tabela[1])  # Resíduo
        tabela1[1][i] = float(tabela[3])  # Probabilidade

    # Configuração do gráfico
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(
        tabela1[0], tabela1[1], "-o", color=cor_linha, label=f"{codigo_pdb}/{chain}", markersize=6
    )

    # Marcar os pontos de agregação forte e fraca
    for i in range(tamanho):
        if tabela1[1][i] >= limite_superior:
            ax.plot(
                tabela1[0][i], tabela1[1][i], "o", color=cor_forte, label="Agregação Forte" if i == 0 else "", markersize=10
            )
        elif tabela1[1][i] <= limite_inferior:
            ax.plot(
                tabela1[0][i], tabela1[1][i], "o", color=cor_fraca, label="Agregação Fraca" if i == 0 else "", markersize=10
            )

    # Adicionar linhas de limite e limiar
    ax.axhline(limite_superior, color="blue", linestyle="--", label="Limite Superior")
    ax.axhline(limite_inferior, color="orange", linestyle="--", label="Limite Inferior")
    ax.axhline(limiar, color="green", linestyle="-", label="Limiar")

    # Configurações do eixo e título
    ax.set_xlabel("Resíduo", fontsize=12, fontweight="bold")
    ax.set_ylabel("Probabilidade", fontsize=12, fontweight="bold")
    ax.set_title(f"Propensão de Agregação: {codigo_pdb} - {chain}", fontsize=14, fontweight="bold")
    ax.set_xlim(min(tabela1[0]) - 1, max(tabela1[0]) + 1)
    ax.set_ylim(0, 1.05)

    # Remover as bordas superior e direita
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Ajustar a legenda
    handles, labels = ax.get_legend_handles_labels()

    # Adicionar as legendas personalizadas na nova ordem
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cor_linha, markersize=6),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cor_forte, markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cor_fraca, markersize=10),
        plt.Line2D([0], [0], color='blue', linestyle='--', markersize=10),
        plt.Line2D([0], [0], color='orange', linestyle='--', markersize=10),
        plt.Line2D([0], [0], color='green', linestyle='-', markersize=10),
    ]
    labels = [
        f"{codigo_pdb}/{chain}",
        'Agregação Forte',
        'Agregação Fraca',
        'Limite Superior',
        'Limite Inferior',
        'Limiar'
    ]

    ax.legend(
        handles=handles,
        labels=labels,
        fontsize=10,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        borderaxespad=0.0,
    )

    st.pyplot(fig)

# Função principal
def main():
    st.title("Gerador de Gráfico de Agregação Estático")

    st.subheader("Como Funciona")
    st.markdown(
        """
        <p style='text-align: justify;'>O Gerador de Gráfico de Agregação Estático é uma ferramenta para visualizar 
        e analisar dados relacionados à agregação de proteínas. A agregação de proteínas é um fenômeno importante que pode 
        levar a doenças neurodegenerativas, como Alzheimer e Parkinson.</p>
        """,
        unsafe_allow_html=True,
    )

    codigo_pdb = st.text_input("Informe o código PDB da proteína:")

    if st.button("Gerar Gráfico"):
        if not codigo_pdb:
            st.warning("Por favor, insira um código PDB válido.")
            return

        global chain  # Declarar chain como global para atualizá-lo

        # Diretório base dinâmico, usa o diretório atual do script ou app
        base_dir = os.getcwd()

        # Monta os caminhos completos com os diretórios relativos "arquivos/pdb" e "arquivos/testes"
        pdb_dir = os.path.join(base_dir, "arquivos", "pdb")
        testes_dir = os.path.join(base_dir, "arquivos", "testes")

        # Garante que os diretórios existam
        os.makedirs(pdb_dir, exist_ok=True)
        os.makedirs(testes_dir, exist_ok=True)

        pdbr = f"https://files.rcsb.org/view/{codigo_pdb}.pdb"

        st.write("Código Solicitado:", codigo_pdb)

        # Download do arquivo PDB
        response = requests.get(pdbr)
        entrada = saida = os.path.join(pdb_dir, f"{codigo_pdb}.pdb")

        if response.status_code == 200:
            with open(saida, "w") as arq_saida:
                arq_saida.write(response.text)
        else:
            st.error("Não foi possível baixar o arquivo PDB.")
            return

        with open(entrada) as arq_entrada:
            pdbLines = arq_entrada.readlines()

        resultado = []
        for line in pdbLines:
            if line.startswith("ATOM") and line[21:22] != chain:
                chain = line[21:22]

                resultado = Prog_Funcoes1.movimenta(
                    f"{codigo_pdb}.pdb", pdbLines, "", chain, codigo_pdb, 1
                )
                resultado = Prog_Funcoes1.avalia_ruido(resultado)

                contatohpFilename = os.path.join(testes_dir, f"{codigo_pdb}_{chain}.csv")
                with open(contatohpFilename, "w") as contatohpFile:
                    contatohpFile.writelines(f"{line}\n" for line in resultado)

                # Evitar gerar gráfico repetido para o mesmo código e chain
                if (codigo_pdb, chain) in graficos_gerados:
                    st.warning("O gráfico já foi gerado para esta cadeia do código PDB.")
                    return

                Plota_Resultado(resultado, codigo_pdb, chain)
                graficos_gerados.add((codigo_pdb, chain))


if __name__ == "__main__":
    main()
