import streamlit as st
import requests
import matplotlib.pyplot as plt
import os
import Prog_Funcoes1

# Função para plotar o gráfico
def Plota_Resultado(resultado, codigo_pdb, chain):
    if not resultado:
        st.warning(f"Sem dados para plotar para a cadeia {chain} do PDB {codigo_pdb}.")
        return

    cor_linha = "black"
    cor_forte = "red"
    cor_fraca = "orange"
    limite_superior = 0.8
    limite_inferior = 0.2
    limiar = 0.5

    tamanho = len(resultado)
    tabela1 = [[0] * tamanho for _ in range(2)]

    for i in range(tamanho):
        tabela = resultado[i].split(";")
        tabela1[0][i] = int(tabela[1])
        tabela1[1][i] = float(tabela[3])

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(tabela1[0], tabela1[1], "-o", color=cor_linha, label=f"{codigo_pdb}/{chain}", markersize=6)

    for i in range(tamanho):
        if tabela1[1][i] >= limite_superior:
            ax.plot(tabela1[0][i], tabela1[1][i], "o", color=cor_forte, label="Agregação Forte" if i == 0 else "", markersize=10)
        elif tabela1[1][i] <= limite_inferior:
            ax.plot(tabela1[0][i], tabela1[1][i], "o", color=cor_fraca, label="Agregação Fraca" if i == 0 else "", markersize=10)

    ax.axhline(limite_superior, color="blue", linestyle="--", label="Limite Superior")
    ax.axhline(limite_inferior, color="orange", linestyle="--", label="Limite Inferior")
    ax.axhline(limiar, color="green", linestyle="-", label="Limiar")

    ax.set_xlabel("Resíduo", fontsize=12, fontweight="bold")
    ax.set_ylabel("Probabilidade", fontsize=12, fontweight="bold")
    ax.set_title(f"Propensão de Agregação: {codigo_pdb} - {chain}", fontsize=14, fontweight="bold")

    # Ajuste seguro do limite do eixo x
    if tabela1[0]:
        ax.set_xlim(min(tabela1[0]) - 1, max(tabela1[0]) + 1)
    else:
        ax.set_xlim(0, 1)

    ax.set_ylim(0, 1.05)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

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

    ax.legend(handles=handles, labels=labels, fontsize=10, loc="upper left", bbox_to_anchor=(1.05, 1))

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

    if "graficos_gerados" not in st.session_state:
        st.session_state["graficos_gerados"] = set()

    codigo_pdb = st.text_input("Informe o código PDB da proteína:").upper()

    if st.button("Gerar Gráfico"):
        if not codigo_pdb:
            st.warning("Por favor, insira um código PDB válido.")
            return

        base_dir = os.getcwd()
        pdb_dir = os.path.join(base_dir, "arquivos", "pdb")
        testes_dir = os.path.join(base_dir, "arquivos", "testes")

        os.makedirs(pdb_dir, exist_ok=True)
        os.makedirs(testes_dir, exist_ok=True)

        url_pdb = f"https://files.rcsb.org/view/{codigo_pdb}.pdb"
        entrada = saida = os.path.join(pdb_dir, f"{codigo_pdb}.pdb")

        st.write("Código solicitado:", codigo_pdb)

        response = requests.get(url_pdb)
        if response.status_code != 200:
            st.error("Não foi possível baixar o arquivo PDB.")
            return

        with open(saida, "w") as arq_saida:
            arq_saida.write(response.text)

        with open(entrada) as arq_entrada:
            pdbLines = arq_entrada.readlines()

        # Detectar todas as cadeias únicas
        cadeias = sorted(set(line[21] for line in pdbLines if line.startswith("ATOM")))

        if not cadeias:
            st.warning("Nenhuma cadeia encontrada.")
            return

        st.success(f"Cadeias detectadas: {', '.join(cadeias)}")

        resultados_por_chain = {}
        for chain in cadeias:
            resultado = Prog_Funcoes1.movimenta(f"{codigo_pdb}.pdb", pdbLines, "", chain, codigo_pdb, 1)
            resultado = Prog_Funcoes1.avalia_ruido(resultado)
            resultados_por_chain[chain] = resultado

        for chain, resultado in resultados_por_chain.items():
            contatohpFilename = os.path.join(testes_dir, f"{codigo_pdb}_{chain}.csv")
            with open(contatohpFilename, "w") as contatohpFile:
                contatohpFile.writelines(f"{line}\n" for line in resultado)

            if (codigo_pdb, chain) in st.session_state["graficos_gerados"]:
                st.info(f"Gráfico anterior para {codigo_pdb}/{chain} será substituído.")

            Plota_Resultado(resultado, codigo_pdb, chain)
            st.session_state["graficos_gerados"].add((codigo_pdb, chain))

            # Botão para download
            with open(contatohpFilename, "r") as f:
                csv_data = f.read()

            st.download_button(
                label=f"📥 Baixar CSV de resultados {codigo_pdb}_{chain}",
                data=csv_data,
                file_name=f"{codigo_pdb}_{chain}.csv",
                mime="text/csv"
            )

if __name__ == "__main__":
    main()
