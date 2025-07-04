import streamlit as st
import pandas as pd
import numpy as np
import pickle
import requests
import os
import io
from sklearn.impute import SimpleImputer
import plotly.graph_objects as go

# --- Configurações de arquivos ---
MODEL_DIR = 'arquivos/modelos'
MODEL_FILENAME = 'modelo_random_forest.pkl'
MODEL_PATH = os.path.join(MODEL_DIR, MODEL_FILENAME)
DESCRIPTORS_PATH = 'arquivos/saida/tabelaDescritoresNew_28-07-2024.csv'
RESULTS_DIR = os.path.join('arquivos', 'testes')
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- Carregamento do modelo ---
clasificador_carregado = None
try:
    with open(MODEL_PATH, 'rb') as arquivo:
        clasificador_carregado = pickle.load(arquivo)
except FileNotFoundError:
    st.error(f"Erro: Arquivo do modelo não encontrado em '{MODEL_PATH}'.")
except Exception as e:
    st.error(f"Erro ao carregar modelo: {e}")

# --- Funções auxiliares ---

def preprocess_data(amino_seq):
    seq_df = pd.DataFrame(list(amino_seq), columns=["aminoacido"])
    dados_csv = pd.read_csv(DESCRIPTORS_PATH, sep=';')
    dados_csv = dados_csv.rename(columns={"Aminoacido": "aminoacido"})
    merge = pd.merge(seq_df, dados_csv, on='aminoacido', how='left')

    merge['Tamanho(3)'] = merge['Tamanho(3)'].map({"Grande": 3, "Medio": 2, "Pequeno": 1})
    merge['Carga(3)'] = merge['Carga(3)'].map({"Positivo": 3, "Neutro": 2, "Negativo": 1})
    df = merge.drop(columns=["aminoacido"])

    df['Prop-Alfa(2)'] = df['Prop-Alfa(2)'].astype(float)
    df['Prop-Beta(2)'] = df['Prop-Beta(2)'].astype(float)

    all_amino_acids = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 
                       'His', 'Ile', 'Leu', 'Lys', 'Met', 'Pre', 'Pro', 
                       'Ser', 'Thr', 'Trp', 'Typ', 'Val']
    df_encoded = pd.get_dummies(df['Simbolo'], prefix='Simbolo', dtype=int).reindex(
        columns=[f'Simbolo_{aa}' for aa in all_amino_acids], fill_value=0)
    df_final = pd.concat([df, df_encoded], axis=1).drop('Simbolo', axis=1)

    tamanhos_janela = [-6,-5,-4,-3,-2,-1,1,2,3,4,5,6]
    colunas_para_deslizamento = df_final.columns[df_final.columns != 'rotulo']
    resultados = pd.DataFrame()

    for coluna in colunas_para_deslizamento:
        for tamanho in tamanhos_janela:
            resultados[f"{coluna}_janela_{tamanho}"] = df_final[coluna].shift(tamanho)

    final_result = pd.concat([df_final, resultados], axis=1).fillna(99)
    imputer = SimpleImputer(strategy='mean')
    final_result = pd.DataFrame(imputer.fit_transform(final_result), columns=final_result.columns)
    return final_result

def exibe_sequencia(probs, seq):
    agregam_count = 0
    nao_agregam_count = 0
    sequence_elements = []

    for i, val in enumerate(probs):
        color = "red" if val > 0.5 else "black"
        tooltip = "Este aminoácido agrega." if val > 0.5 else "Este aminoácido não agrega."
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
        sequence_elements.append(f'<span style="{span_style}" title="{tooltip}">{seq[i]}</span>')

        if val > 0.5:
            agregam_count += 1
        else:
            nao_agregam_count += 1

    st.header("Sequência de Aminoácidos")
    st.markdown(" ".join(sequence_elements), unsafe_allow_html=True)
    st.header("Resumo")
    st.write(f"Número de aminoácido(s) que agregam: {agregam_count}")
    st.write(f"Número de aminoácido(s) que não agregam: {nao_agregam_count}")

def plota_resultado_plotly(probs, seq, titulo="Sequência"):
    tamanho = len(probs)
    posicoes = list(range(1, tamanho + 1))

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=posicoes, y=probs,
        mode="lines+markers",
        line=dict(color="black"),
        name=f"{titulo}"
    ))

    fig.add_trace(go.Scatter(
        x=[posicoes[i] for i in range(tamanho) if probs[i] > 0.5],
        y=[probs[i] for i in range(tamanho) if probs[i] > 0.5],
        mode="markers",
        marker=dict(color="red", size=10),
        name="Agregam (> 0.5)"
    ))

    fig.add_trace(go.Scatter(
        x=[posicoes[i] for i in range(tamanho) if probs[i] <= 0.5],
        y=[probs[i] for i in range(tamanho) if probs[i] <= 0.5],
        mode="markers",
        marker=dict(color="gray", size=8),
        name="Não agregam (≤ 0.5)"
    ))

    fig.add_trace(go.Scatter(
        x=posicoes,
        y=[0.3]*tamanho,
        mode="lines",
        line=dict(color="blue", dash="dash"),
        name="Limite Inferior"
    ))

    fig.add_trace(go.Scatter(
        x=posicoes,
        y=[0.7]*tamanho,
        mode="lines",
        line=dict(color="orange", dash="dash"),
        name="Limite Superior"
    ))

    fig.add_trace(go.Scatter(
        x=posicoes,
        y=[0.5]*tamanho,
        mode="lines",
        line=dict(color="green"),
        name="Limiar 0.5"
    ))

    fig.update_layout(
        title=f"Propensão à agregação / {titulo}",
        xaxis_title="Resíduo",
        yaxis_title="Probabilidade",
        yaxis=dict(range=[0,1]),
        hovermode="closest"
    )

    st.plotly_chart(fig)

# Dicionário para converter aminoácidos 3-letras para 1-letra
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# --- Aplicativo principal ---
def main():
    st.title("🧬 Classificador de Sequência de Aminoácidos")

    aba = st.radio("Escolha a forma de entrada da sequência:", ["Via código PDB", "Sequência manual", "Upload CSV"])

    if aba == "Via código PDB":
        codigo_pdb = st.text_input("Código PDB (ex: 1XQ8)").upper()
        if st.button("Buscar e Classificar"):
            if codigo_pdb:
                # Baixar arquivo PDB completo
                url_pdb = f'https://files.rcsb.org/view/{codigo_pdb}.pdb'
                response = requests.get(url_pdb)
                if response.status_code != 200:
                    st.error(f"Erro ao baixar o arquivo PDB de {codigo_pdb}")
                    return

                pdb_text = response.text
                pdb_lines = pdb_text.splitlines()

                # Extrair todas as cadeias únicas
                cadeias = sorted(set(line[21] for line in pdb_lines if line.startswith("ATOM")))

                if not cadeias:
                    st.warning("Nenhuma cadeia encontrada no arquivo PDB.")
                    return

                st.success(f"Cadeias encontradas: {', '.join(cadeias)}")

                for chain in cadeias:
                    # Extrair sequência da cadeia do arquivo PDB
                    seq_residues = []
                    last_res_seq = None
                    for line in pdb_lines:
                        if line.startswith("ATOM") and line[21] == chain:
                            res_seq = line[22:26].strip()
                            aa = line[17:20].strip()
                            # Evitar repetir resíduos (mesmo número sequencial)
                            if res_seq != last_res_seq:
                                seq_residues.append(aa)
                                last_res_seq = res_seq

                    seq = ''.join(three_to_one.get(aa, 'X') for aa in seq_residues)  # 'X' para desconhecidos

                    st.subheader(f"Cadeia {chain}")
                    st.code(seq, language='text')

                    try:
                        df_input = preprocess_data(seq)
                        probs = clasificador_carregado.predict_proba(df_input)[:, 1]

                        exibe_sequencia(probs, seq)
                        plota_resultado_plotly(probs, seq, f"{codigo_pdb}/{chain}")

                        # Preparar CSV para salvar e download
                        csv_data = pd.DataFrame({
                            "resíduo": list(range(1, len(probs) + 1)),
                            "probabilidade": probs
                        })

                        csv_path = os.path.join(RESULTS_DIR, f"{codigo_pdb}_{chain}.csv")
                        csv_data.to_csv(csv_path, index=False)

                        csv_bytes = csv_data.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label=f"📥 Baixar CSV: {codigo_pdb}_{chain}",
                            data=csv_bytes,
                            file_name=f"{codigo_pdb}_{chain}.csv",
                            mime="text/csv"
                        )

                    except Exception as e:
                        st.error(f"Erro ao processar cadeia {chain}: {e}")
            else:
                st.warning("Informe um código PDB.")

    elif aba == "Sequência manual":
        amino_seq = st.text_area("Digite a sequência de aminoácidos:")
        if st.button("Classificar"):
            if amino_seq:
                try:
                    df_input = preprocess_data(amino_seq)
                    probs = clasificador_carregado.predict_proba(df_input)[:, 1]
                    exibe_sequencia(probs, amino_seq)
                    plota_resultado_plotly(probs, amino_seq, "Sequência Manual")
                except Exception as e:
                    st.error(f"Erro ao classificar sequência: {e}")
            else:
                st.warning("Digite uma sequência.")

    elif aba == "Upload CSV":
        uploaded_file = st.file_uploader("Selecione um arquivo CSV", type=["csv"])
        if uploaded_file is not None:
            try:
                lines = []
                for line in uploaded_file.getvalue().decode("utf-8").split('\n'):
                    if len(line.strip().split(";")) == 5:
                        lines.append(line)

                df = pd.read_csv(io.StringIO("\n".join(lines)), sep=";", encoding="utf-8")

                # Extração da sequência da terceira coluna
                seq = "".join(df.iloc[:, 2])

                # Pré-processamento, predição e visualização
                df_preprocessed = preprocess_data(seq)
                probs = clasificador_carregado.predict_proba(df_preprocessed)[:, 1]

                exibe_sequencia(probs, seq)
                plota_resultado_plotly(probs, seq, titulo="Upload CSV")

            except Exception as e:
                st.error(f"Erro ao processar o arquivo CSV: {e}")

if __name__ == "__main__":
    main()
