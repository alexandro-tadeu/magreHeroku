import streamlit as st
import pandas as pd
import numpy as np
import pickle
from sklearn.impute import SimpleImputer

# Carregar o modelo salvo com pickle
with open('arquivos/modelos/modelo_random_forest.pkl', 'rb') as arquivo:
    clasificador_carregado = pickle.load(arquivo)
#st.write("Modelo carregado com sucesso usando pickle!")

# Função para processar e preparar os dados de entrada
def preprocess_data(amino_seq):
    # Convertendo a sequência em DataFrame
    seq_df = pd.DataFrame(list(amino_seq), columns=["aminoacido"])
    
    # Carregar os dados da tabela de descritores
    dados_csv = pd.read_csv('arquivos/saida/tabelaDescritoresNew_28-07-2024.csv', sep=';')
    dados_csv = dados_csv.rename(columns={"Aminoacido": "aminoacido"})
    
    # Juntar os dados da sequência com o arquivo de descritores
    merge = pd.merge(seq_df, dados_csv, on='aminoacido', how='left')

    # Mapear valores categóricos
    merge['Tamanho(3)'] = merge['Tamanho(3)'].map({"Grande": 3, "Medio": 2, "Pequeno": 1})
    merge['Carga(3)'] = merge['Carga(3)'].map({"Positivo": 3, "Neutro": 2, "Negativo": 1})
    
    # Removendo a coluna aminoacido
    df = merge.drop(columns=["aminoacido"])

    # Convertendo para tipo de dado float
    df['Prop-Alfa(2)'] = df['Prop-Alfa(2)'].astype(float)
    df['Prop-Beta(2)'] = df['Prop-Beta(2)'].astype(float)

    # Codificação one-hot dos aminoácidos
    all_amino_acids = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 
                       'His', 'Ile', 'Leu', 'Lys', 'Met', 'Pre', 'Pro', 
                       'Ser', 'Thr', 'Trp', 'Typ', 'Val']
    
    df_encoded = pd.get_dummies(df['Simbolo'], prefix='Simbolo', dtype=int).reindex(
        columns=[f'Simbolo_{aa}' for aa in all_amino_acids], fill_value=0)
    
    # Concatenar o DataFrame original com o codificado
    df_final = pd.concat([df, df_encoded], axis=1).drop('Simbolo', axis=1)

    # Aplicando a janela deslizante
    tamanhos_janela = [-6,-5,-4,-3,-2,-1,1,2,3,4,5,6]  # Tamanhos de janela
    colunas_para_deslizamento = df_final.columns[df_final.columns != 'rotulo']
    resultados = pd.DataFrame()

    for coluna in colunas_para_deslizamento:
        for tamanho in tamanhos_janela:
            janela_deslizante = df_final[coluna].shift(tamanho)
            novo_nome_coluna = f"{coluna}_janela_{tamanho}"
            resultados[novo_nome_coluna] = janela_deslizante

    # Concatenar colunas originais e janelas deslizantes
    final_result = pd.concat([df_final, resultados], axis=1).fillna(99)

    # Preencher valores nulos
    imputer = SimpleImputer(strategy='mean')
    final_result = pd.DataFrame(imputer.fit_transform(final_result), columns=final_result.columns)

    return final_result

# Função para mostrar a sequência colorida
def display_sequence_with_colors(amino_seq, agregam_count):
    colored_sequence = []
    agregam_count_value = 0
    nao_agregam_count_value = 0

    for i, aa in enumerate(amino_seq):
        # Condição para determinar se o aminoácido agrega ou não
        if agregam_count[i] == 1:  # Se o valor de agregação for 1 (Agrega)
            color = "red"  # Cor vermelha para quem agrega
            tooltip = "Este aminoácido agrega."
            agregam_count_value += 1
        else:  # Se o valor de agregação for 0 (Não Agrega)
            color = "black"  # Cor preta para quem não agrega
            tooltip = "Este aminoácido não agrega."
            nao_agregam_count_value += 1

        # Estilo extra para cada aminoácido
        span_style = f"color: {color}; font-size: 20px; background-color: #F0F0F0; border: 1px solid #CCCCCC; border-radius: 4px; padding: 2px 6px; margin: 2px"

        # Adiciona o aminoácido com o estilo e o tooltip
        colored_sequence.append(
            f'<span style="{span_style}" title="{tooltip}">{aa}</span>'
        )

    # Exibir a sequência com as cores e o estilo desejado
    st.markdown(''.join(colored_sequence), unsafe_allow_html=True)

    return agregam_count_value, nao_agregam_count_value

# Função principal
def main():
    # Título da aplicação
    st.title("Classificação de Sequência de Aminoácidos")

    # Entrada de texto para a sequência de aminoácidos
    amino_seq = st.text_area("Digite a sequência de aminoácidos:")

    # Botão para iniciar a classificação
    if st.button("Classificar"):
        if amino_seq:
            try:
                # Processar a sequência
                final_result = preprocess_data(amino_seq)
                
                # Fazer a previsão
                new_predictions = clasificador_carregado.predict(final_result)
                
                # Exibir os resultados
                st.text_area("Resultado da Classificação:", value=' '.join(map(str, new_predictions.flatten())))

                # Exibir a sequência com cores e contar as letras que agregam e não agregam
                agregam_count_value, nao_agregam_count_value = display_sequence_with_colors(amino_seq, new_predictions)

                # Exibir o resumo
                st.header("Resumo")
                st.write(f"Número de aminoácido(s) que agregam: {agregam_count_value}")
                st.write(f"Número de aminoácido(s) que não agregam: {nao_agregam_count_value}")

            except Exception as e:
                st.error(f"Ocorreu um erro durante a classificação: {e}")
        else:
            st.warning("Por favor, insira uma sequência de aminoácidos.")
