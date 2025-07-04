import streamlit as st
import pandas as pd
import io
import plotly.graph_objects as go
from Bio import SeqIO

def plota_resultado_plotly(df):
    tamanho = len(df)
    tabela1 = [[0] * tamanho for _ in range(2)]

    for i in range(tamanho):
        tabela1[0][i] = i + 1  # posição do aminoácido (índice + 1)
        tabela1[1][i] = float(df.iloc[i, 3])  # valor da propensão na 4ª coluna

    limite_inferior, limite_superior = 0.2, 0.8
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=tabela1[0],
        y=tabela1[1],
        mode="lines+markers",
        line=dict(color="black"),
        name="Propensão de agregação"
    ))
    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        mode="markers",
        marker=dict(color="red", size=10),
        name="Agregação Forte"
    ))
    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        mode="markers",
        marker=dict(color="orange", size=10),
        name="Agregação Fraca"
    ))
    fig.update_layout(
        title="Propensão de Agregação por Resíduo",
        xaxis_title="Posição do Aminoácido",
        yaxis_title="Probabilidade"
    )
    st.plotly_chart(fig)

def main():
    # Título da página
    st.title("Gerador de Sequência de Aminoácidos")

    st.markdown(
        """
        <p style='text-align: justify;'>
        O Gerador de Sequência de Aminoácidos é uma ferramenta essencial na pesquisa em biologia molecular e bioinformática. Aminoácidos são os blocos de construção fundamentais das proteínas e compreender suas sequências é crucial para o entendimento de processos biológicos complexos.
        </p>
        """,
        unsafe_allow_html=True,
    )

    # Subtítulo
    st.subheader("Como Funciona")

    st.markdown(
    """
    <p style='text-align: justify;'>
    O funcionamento do Gerador de Sequência de Aminoácidos é direto e eficiente. Os usuários seleciona um aquivo .csv proteína, contendo os dados preditos do PDB (Protein Data Bank) contendo a  sequência de aminoácidos. Em seguida, a ferramenta gera uma representação visual da sequência, destacando características importantes, como regiões de agregação propença a agregação.
    </p>
    """,
    unsafe_allow_html=True,
    )

    st.subheader("Visualização de Dados")

    st.markdown(
    """
    <p style='text-align: justify;'>
    A capacidade do Gerador de Sequência de Aminoácidos de criar representações visuais claras e informativas das sequências torna-o inestimável para pesquisadores e cientistas. Essas visualizações podem ser fundamentais para a identificação de padrões, análise de estruturas e tomada de decisões informadas em estudos envolvendo proteínas e suas funções.
    </p>
    """,
    unsafe_allow_html=True,
    )


    # Seletor de arquivo CSV
    uploaded_file = st.file_uploader("Selecione um arquivo CSV", type=["csv"])

    if uploaded_file is not None:
        try:
            # Ler o arquivo CSV e filtrar linhas com número incorreto de campos
            lines = []
            for line in uploaded_file.getvalue().decode().split('\n'):
                fields = line.strip().split(";")
                if len(fields) == 5:
                    lines.append(line)

            # Criar um DataFrame com as linhas válidas
            df = pd.read_csv(io.StringIO("\n".join(lines)), sep=";", encoding="utf-8")

            # Extraindo a terceira coluna do DataFrame como uma sequência de caracteres
            seq = "".join(df.iloc[:, 2])

            # Criação do componente SequenceViewer
            agregam_count = 0
            nao_agregam_count = 0

            # Criando uma lista de elementos HTML para representar cada letra com a cor apropriada e a dica de ferramenta (tooltip)
            sequence_elements = []
            for index, val in enumerate(df.iloc[:, 3]):
                # Define a cor preta para valor maior que 0.5 (Agrega) e vermelha para valor menor ou igual a 0.5 (Não Agrega)
                color = "red" if val > 0.5 else "black"
                tooltip = "Este aminoácido agrega." if val > 0.5 else "Este aminoácido não agrega."

                span_style = f"color: {color}; font-size: 20px; background-color: #F0F0F0; border: 1px solid #CCCCCC; border-radius: 4px; padding: 2px 6px; margin: 2px"

                sequence_elements.append(
                    f'<span style="{span_style}" title="{tooltip}">{seq[index]}</span>'
                )

                # Contando as letras que agregam e as que não agregam
                if val > 0.5:
                    agregam_count += 1
                else:
                    nao_agregam_count += 1

            # Exibindo a sequência gerada e o resumo
            st.header("Sequência de Aminoácidos")
            sequence_html = " ".join(sequence_elements)
            st.markdown(sequence_html, unsafe_allow_html=True)

            st.header("Resumo")
            st.write(f"Número de aminoácido(s) que agregam: {agregam_count}")
            st.write(f"Número de aminoácido(s) que não agregam: {nao_agregam_count}")

            # Chamada da função para gerar o gráfico
            plota_resultado_plotly(df)

        except Exception as e:
            st.error(f"Erro ao processar o arquivo CSV: {str(e)}")

if __name__ == "__main__":
    main()
