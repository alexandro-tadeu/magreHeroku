import streamlit as st
import requests
import os
from pathlib import Path
import plotly.graph_objects as go
import Prog_Funcoes1

# Variáveis de controle e diretórios
graficos_gerados = set()

# Base dir dinâmico: tenta ler da variável ambiente MAGRE_BASE_DIR, senão usa a pasta atual
BASE_DIR = Path(os.getenv("MAGRE_BASE_DIR", ".")) / "arquivos"

# Diretórios específicos
arquivo_pdb = BASE_DIR / "pdb"
arquivo_testes = BASE_DIR / "testes"

# Certifique-se de que os diretórios existem
arquivo_pdb.mkdir(parents=True, exist_ok=True)
arquivo_testes.mkdir(parents=True, exist_ok=True)

# Função para baixar o arquivo PDB
def baixa_pdb(codigo_pdb):
    url = f"https://files.rcsb.org/view/{codigo_pdb}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        st.write("Arquivo PDB baixado com sucesso.")
        return response.text
    else:
        st.warning("Não foi possível baixar o arquivo PDB. Verifique o código.")
        return None

# Função para salvar o arquivo PDB
def salva_pdb(codigo_pdb, conteudo):
    caminho = arquivo_pdb / f"{codigo_pdb}.pdb"
    with open(caminho, "w") as f:
        f.write(conteudo)
    return caminho

# Função para processar uma cadeia específica no arquivo PDB
def processa_cadeia(codigo_pdb, linhas, cadeia):
    resultado = Prog_Funcoes1.movimenta(
        f"{codigo_pdb}.pdb", linhas, "", cadeia, codigo_pdb, 1
    )
    resultado = Prog_Funcoes1.avalia_ruido(resultado)
    return resultado

# Função para exibir sequência colorida
def exibe_sequencia(resultado):
    st.subheader("Detecção de Agregação")
    agregam_count = 0
    nao_agregam_count = 0

    sequence_elements = []
    for linha in resultado:
        dados = linha.split(";")
        if len(dados) < 5:
            continue
        posicao = int(dados[1])
        valor = float(dados[3])
        residuo = dados[2].strip()

        color = "red" if valor > 0.5 else "black"
        tooltip = "Este aminoácido agrega." if valor > 0.5 else "Este aminoácido não agrega."

        span_style = (
            f"color: {color}; "
            "font-size: 20px; "
            "background-color: #F0F0F0; "
            "border: 1px solid #CCCCCC; "
            "border-radius: 4px; "
            "padding: 2px 6px; "
            "margin: 2px;"
            "display: inline-block;"
            "text-align: center;"
            "width: 24px;"
        )
        sequence_elements.append(
            f'<span style="{span_style}" title="{tooltip}">{residuo}</span>'
        )

        if valor > 0.5:
            agregam_count += 1
        else:
            nao_agregam_count += 1

    sequence_html = "".join(sequence_elements)
    st.markdown(sequence_html, unsafe_allow_html=True)

    st.subheader("Resumo")
    st.write(f"Número de aminoácidos que agregam: {agregam_count}")
    st.write(f"Número de aminoácidos que não agregam: {nao_agregam_count}")

# Função para plotar resultado com Plotly
def plota_resultado_plotly(resultado, codigo_pdb, chain):
    tamanho = len(resultado)
    tabela1 = [[0] * tamanho for _ in range(2)]

    for i, linha in enumerate(resultado):
        tabela = linha.split(";")
        tabela1[0][i] = int(tabela[1])
        tabela1[1][i] = float(tabela[3])

    limite_inferior, limite_superior = 0.2, 0.8
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=tabela1[0], y=tabela1[1], mode="lines+markers",
                             line=dict(color="black"), name=f"{codigo_pdb}/{chain}"))
    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        mode="markers", marker=dict(color="red", size=10), name="Agregação Forte"))
    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        mode="markers", marker=dict(color="orange", size=10), name="Agregação Fraca"))
    fig.update_layout(
        title=f"Propensão de agregação/{codigo_pdb}-{chain}",
        xaxis_title="Resíduo",
        yaxis_title="Probabilidade"
    )
    st.plotly_chart(fig)

# Função principal
def main():
    st.title("Gerador de Gráfico e Detecção de Agregação")

    codigo_pdb = st.text_input("Informe o código PDB da proteína:").upper().strip()

    if st.button("Gerar Análise"):
        if not codigo_pdb:
            st.warning("Por favor, insira um código PDB válido.")
            return

        conteudo_pdb = baixa_pdb(codigo_pdb)
        if conteudo_pdb is None:
            return

        salva_pdb(codigo_pdb, conteudo_pdb)
        linhas = conteudo_pdb.splitlines()

        # Detectar todas as cadeias únicas no PDB
        cadeias = sorted(set(line[21] for line in linhas if line.startswith("ATOM")))

        if not cadeias:
            st.warning("Nenhuma cadeia encontrada no arquivo PDB.")
            return

        st.success(f"Cadeias detectadas: {', '.join(cadeias)}")

        for chain in cadeias:
            if (codigo_pdb, chain) in graficos_gerados:
                st.info(f"Análise para {codigo_pdb}/{chain} já gerada anteriormente.")
                continue

            resultado = processa_cadeia(codigo_pdb, linhas, chain)

            if not resultado:
                st.warning(f"Não foi possível processar a cadeia {chain}.")
                continue

            # Salvar CSV dos resultados
            arquivo_csv = arquivo_testes / f"{codigo_pdb}_{chain}.csv"
            with open(arquivo_csv, "w") as f:
                for linha in resultado:
                    f.write(linha + "\n")

            # Exibir gráfico e sequência
            plota_resultado_plotly(resultado, codigo_pdb, chain)
            exibe_sequencia(resultado)

            # Botão para baixar CSV (ajustado para não atualizar a página)
            with open(arquivo_csv, "r") as f:
                csv_data = f.read()

            with st.container():
                st.download_button(
                    label=f"📥 Baixar CSV da cadeia {chain}",
                    data=csv_data,
                    file_name=f"{codigo_pdb}_{chain}.csv",
                    mime="text/csv",
                    key=f"download_{codigo_pdb}_{chain}"  # chave única evita refresh
                )

            graficos_gerados.add((codigo_pdb, chain))

if __name__ == "__main__":
    main()
