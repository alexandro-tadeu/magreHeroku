import streamlit as st
import os
import plotly.graph_objects as go
import requests

import Prog_Funcoes1

# --- Configurações de arquivos ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.path.join(BASE_DIR, 'arquivos', 'pdb')
TESTES_DIR = os.path.join(BASE_DIR, 'arquivos', 'testes')
os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(TESTES_DIR, exist_ok=True)

def Plota_Resultado(resultado1, CodigoPDB1, chain1):
    cor_agregacao = 'red'
    cor_normal = 'black'
    cor_limite_inferior = 'orange'
    cor_limite_superior = 'blue'
    tamanho = len(resultado1)

    tabela1 = [[0] * tamanho for _ in range(2)]

    for i in range(tamanho):
        tabela = resultado1[i].split(";")
        tabela1[0][i] = int(tabela[1])
        tabela1[1][i] = float(tabela[3])

    limite_inferior = 0.2
    limite_superior = 0.8

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=tabela1[0], y=tabela1[1], mode='lines+markers',
        line=dict(color=cor_normal), name=f"{CodigoPDB1}/{chain1}",
        hovertemplate="Resíduo: %{x}<br>Probabilidade: %{y}<extra></extra>"
    ))

    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        mode='markers',
        marker=dict(color=cor_agregacao, size=10),
        name='Agregação Forte',
        hovertemplate="Resíduo: %{x}<br>Probabilidade: %{y}<extra></extra>"
    ))

    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        mode='markers',
        marker=dict(color=cor_limite_inferior, size=10),
        name='Agregação Fraca',
        hovertemplate="Resíduo: %{x}<br>Probabilidade: %{y}<extra></extra>"
    ))

    fig.add_trace(go.Scatter(
        x=tabela1[0], y=[limite_inferior] * tamanho,
        mode='lines', line=dict(color=cor_limite_inferior, dash='dash'), name='Limite Inferior'
    ))

    fig.add_trace(go.Scatter(
        x=tabela1[0], y=[limite_superior] * tamanho,
        mode='lines', line=dict(color=cor_limite_superior, dash='dash'), name='Limite Superior'
    ))

    fig.add_trace(go.Scatter(
        x=tabela1[0], y=[0.5] * tamanho,
        mode='lines', line=dict(color='green'), name='Limiar'
    ))

    fig.update_layout(
        title=f"Propensão à agregação/{CodigoPDB1}-{chain1}",
        xaxis=dict(title='Resíduo', showgrid=True),
        yaxis=dict(title='Probabilidade', showgrid=True),
        showlegend=True,
        hovermode='closest'
    )

    return fig

def main():
    st.title("Gráfico de Propensão à Agregação para Múltiplas Cadeias")

    if "graficos_gerados" not in st.session_state:
        st.session_state["graficos_gerados"] = set()

    codigo_pdb = st.text_input("Informe o código PDB da proteína:").upper()

    if st.button("Gerar Gráfico"):
        if not codigo_pdb:
            st.warning("Por favor, insira um código PDB válido.")
            return

        st.write("Código solicitado:", codigo_pdb)

        url_pdb = f'https://files.rcsb.org/view/{codigo_pdb}.pdb'
        response = requests.get(url_pdb)

        if response.status_code != 200:
            st.error(f"Erro ao baixar o arquivo PDB de {codigo_pdb}")
            return

        entrada = os.path.join(PDB_DIR, f"{codigo_pdb}.pdb")
        with open(entrada, 'w') as f:
            f.write(response.text)

        with open(entrada) as arq_entrada:
            pdbLines = arq_entrada.readlines()

        # Identificar todas as cadeias únicas
        cadeias = sorted(set(line[21] for line in pdbLines if line.startswith("ATOM")))

        if not cadeias:
            st.warning("Nenhuma cadeia encontrada no arquivo PDB.")
            return

        st.success(f"Cadeias encontradas: {', '.join(cadeias)}")

        resultados_por_chain = {}

        for chain in cadeias:
            resultado = Prog_Funcoes1.movimenta(f"{codigo_pdb}.pdb", pdbLines, "", chain, codigo_pdb, 1)
            resultado = Prog_Funcoes1.avalia_ruido(resultado)
            resultados_por_chain[chain] = resultado

        for chain, resultado in resultados_por_chain.items():
            fig = Plota_Resultado(resultado, codigo_pdb, chain)
            st.plotly_chart(fig)

            csv_path = os.path.join(TESTES_DIR, f"{codigo_pdb}_{chain}.csv")
            with open(csv_path, "w") as f:
                for line in resultado:
                    f.write(line + "\n")

            with open(csv_path, "r") as f:
                csv_data = f.read()

            st.download_button(
                label=f"📥 Baixar CSV: {codigo_pdb}_{chain}",
                data=csv_data,
                file_name=f"{codigo_pdb}_{chain}.csv",
                mime="text/csv"
            )

if __name__ == '__main__':
    main()
