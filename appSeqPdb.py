import requests
import streamlit as st
import os
import io  # para manipular arquivos em memória

def obter_codigo_fasta(codigo_pdb):
    url = f'https://www.rcsb.org/fasta/entry/{codigo_pdb}'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.splitlines()
    else:
        st.error(f"Erro ao obter o código FASTA {codigo_pdb}. Verifique o código e tente novamente.")
        return None

def extrair_sequencias_fasta_multicadeia(fasta_lines):
    sequencias = {}
    cadeia_atual = None
    for line in fasta_lines:
        if line.startswith(">"):
            parts = line.split("|")
            header = parts[0]
            if "_" in header:
                cadeia = header.split("_")[1].strip()
            else:
                cadeia = "N/A"
            cadeia_atual = cadeia
            sequencias[cadeia_atual] = ""
        else:
            if cadeia_atual is not None:
                sequencias[cadeia_atual] += line.strip()
    return sequencias

def baixar_arquivo_pdb(codigo_pdb):
    url_pdb = f"https://files.rcsb.org/view/{codigo_pdb}.pdb"
    response = requests.get(url_pdb)
    if response.status_code == 200:
        pasta_pdb = os.path.join(os.getcwd(), "arquivos", "pdb")
        os.makedirs(pasta_pdb, exist_ok=True)
        caminho_arquivo = os.path.join(pasta_pdb, f"{codigo_pdb}.pdb")
        with open(caminho_arquivo, "w") as f:
            f.write(response.text)
        return caminho_arquivo
    else:
        st.error("Não foi possível baixar o arquivo PDB.")
        return None

def identificar_cadeias_pdb(caminho_arquivo_pdb):
    cadeias = set()
    try:
        with open(caminho_arquivo_pdb, "r") as f:
            for linha in f:
                if linha.startswith("ATOM") or linha.startswith("HETATM"):
                    if len(linha) > 21:
                        cadeia = linha[21]
                        cadeias.add(cadeia)
        return sorted(list(cadeias))
    except Exception as e:
        st.error(f"Erro ao ler arquivo PDB: {e}")
        return []

import plotly.graph_objects as go

def plota_resultado_plotly(posicoes, valores, codigo_pdb, cadeia):
    limite_inferior, limite_superior = 0.2, 0.8
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=posicoes, y=valores, mode="lines+markers",
                             line=dict(color="black"), name=f"{codigo_pdb}/{cadeia}"))
    fig.add_trace(go.Scatter(
        x=[posicoes[i] for i, v in enumerate(valores) if v >= limite_superior],
        y=[v for v in valores if v >= limite_superior],
        mode="markers", marker=dict(color="red", size=10), name="Agregação Forte"))
    fig.add_trace(go.Scatter(
        x=[posicoes[i] for i, v in enumerate(valores) if v <= limite_inferior],
        y=[v for v in valores if v <= limite_inferior],
        mode="markers", marker=dict(color="orange", size=10), name="Agregação Fraca"))
    fig.update_layout(
        title=f"Propensão de agregação {codigo_pdb} - Cadeia {cadeia}",
        xaxis_title="Posição do Aminoácido",
        yaxis_title="Probabilidade"
    )
    st.plotly_chart(fig)

def processa_cadeia_para_plotagem(caminho_pdb, cadeia):
    import random
    tamanho = 50
    posicoes = list(range(1, tamanho+1))
    valores = [random.uniform(0, 1) for _ in range(tamanho)]
    seq = "".join(["A" for _ in range(tamanho)])
    return posicoes, valores, seq

def main():
    st.title("Extrator de Sequência de Aminoácidos do PDB sem aprendizado")

    codigo_pdb = st.text_input("Insira o código da entrada (por exemplo, 1XQ8):", "").upper()

    if st.button("Buscar Sequência"):
        if codigo_pdb:
            fasta_lines = obter_codigo_fasta(codigo_pdb)
            if fasta_lines:
                sequencias = extrair_sequencias_fasta_multicadeia(fasta_lines)
                if sequencias:
                    st.success("Sequências de Aminoácidos (FASTA) por cadeia:")
                    # Aqui exibimos como "Cadeia 1:", "Cadeia 2:", ...
                    for i, (cadeia, seq) in enumerate(sequencias.items(), start=1):
                        st.markdown(f"**Cadeia {i}:**")
                        st.code(seq)
                else:
                    st.warning("Nenhuma sequência FASTA encontrada.")

                caminho_pdb = baixar_arquivo_pdb(codigo_pdb)
                if caminho_pdb:
                    cadeias = identificar_cadeias_pdb(caminho_pdb)
                    if cadeias:
                        st.write(f"**Cadeias detectadas no arquivo PDB:** {', '.join(cadeias)}")

                        with open(caminho_pdb, "rb") as f:
                            pdb_bytes = f.read()
                        st.download_button(
                            label=f"📥 Baixar arquivo PDB {codigo_pdb}",
                            data=pdb_bytes,
                            file_name=f"{codigo_pdb}.pdb",
                            mime="chemical/x-pdb"
                        )

                        for c in cadeias:
                            st.subheader(f"Cadeia: {c}")

                            posicoes, valores, seq = processa_cadeia_para_plotagem(caminho_pdb, c)

                            sequence_elements = []
                            agregam_count = 0
                            nao_agregam_count = 0
                            for i, val in enumerate(valores):
                                cor = "red" if val > 0.5 else "black"
                                tooltip = "Este aminoácido agrega." if val > 0.5 else "Este aminoácido não agrega."
                                sequence_elements.append(
                                    f'<span style="color:{cor}; font-weight:bold;" title="{tooltip}">{seq[i]}</span>'
                                )
                                if val > 0.5:
                                    agregam_count += 1
                                else:
                                    nao_agregam_count += 1
                            st.markdown("".join(sequence_elements), unsafe_allow_html=True)

                            st.write(f"Nº aminoácidos que agregam: {agregam_count}")
                            st.write(f"Nº aminoácidos que não agregam: {nao_agregam_count}")

                            plota_resultado_plotly(posicoes, valores, codigo_pdb, c)

                            import csv
                            output = io.StringIO()
                            writer = csv.writer(output, delimiter=';')
                            writer.writerow(["Posição", "Valor", "Residuo"])
                            for p, v, r in zip(posicoes, valores, seq):
                                writer.writerow([p, f"{v:.3f}", r])

                            csv_data = output.getvalue()
                            output.close()

                            st.download_button(
                                label=f"📥 Baixar CSV da cadeia {c}",
                                data=csv_data,
                                file_name=f"{codigo_pdb}_{c}.csv",
                                mime="text/csv"
                            )
                    else:
                        st.warning("Nenhuma cadeia identificada no arquivo PDB.")
            else:
                st.warning("Não foi possível obter os dados FASTA.")
        else:
            st.warning("Por favor, insira um código PDB.")

if __name__ == "__main__":
    main()
