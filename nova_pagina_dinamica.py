import streamlit as st
import os
import plotly.graph_objects as go
import webbrowser
import requests

# Importe o módulo Prog_Funcoes1
import Prog_Funcoes1

# Variável de controle para verificar se o gráfico já foi gerado
graficos_gerados = set()

# Defina chain como uma variável global
chain = " "

# Função para plotar o gráfico com mais informações
def Plota_Resultado(resultado1, CodigoPDB1, chain1):
    cor_agregacao = 'red'
    cor_normal = 'black'
    cor_limite_inferior = 'orange'
    cor_limite_superior = 'blue'
    tamanho = len(resultado1)

    tabela1 = [[0] * tamanho for i in range(2)]

    for i in range(tamanho):
        tabela = resultado1[i].split(";")
        tabela1[0][i] = int(tabela[1])
        tabela1[1][i] = float(tabela[3])

    # Definindo os limites de agregação
    limite_inferior = 0.2
    limite_superior = 0.8
    agregacao_forte = [y for y in tabela1[1] if y >= limite_superior]
    agregacao_fraca = [y for y in tabela1[1] if y <= limite_inferior]

    fig = go.Figure()

    # Plotando o gráfico com a linha de agregação
    fig.add_trace(go.Scatter(x=tabela1[0], y=tabela1[1], mode='lines+markers',
                             line=dict(color=cor_normal), name=CodigoPDB1 + "/" + chain1,
                             hovertemplate="Resíduo: %{x}<br>Probabilidade: %{y}<extra></extra>"))

    # Marcando os pontos de agregação forte
    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] >= limite_superior],
        mode='markers',
        marker=dict(color=cor_agregacao, size=10),
        name='Agregação Forte',
        hovertemplate="Resíduo: %{x}<br>Probabilidade: %{y}<extra></extra>"  # Adiciona o hover aqui também
    ))

    # Marcando os pontos de agregação fraca
    fig.add_trace(go.Scatter(
        x=[tabela1[0][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        y=[tabela1[1][i] for i in range(tamanho) if tabela1[1][i] <= limite_inferior],
        mode='markers',
        marker=dict(color=cor_limite_inferior, size=10),
        name='Agregação Fraca',
        hovertemplate="Resíduo: %{x}<br>Probabilidade: %{y}<extra></extra>"  # Adiciona o hover aqui também
    ))

    # Linha de limiar inferior
    fig.add_trace(go.Scatter(x=tabela1[0], y=[limite_inferior] * tamanho,
                             mode='lines', line=dict(color=cor_limite_inferior, dash='dash'), name='Limite Inferior'))

    # Linha de limiar superior
    fig.add_trace(go.Scatter(x=tabela1[0], y=[limite_superior] * tamanho,
                             mode='lines', line=dict(color=cor_limite_superior, dash='dash'), name='Limite Superior'))

    # Linha do limiar padrão
    fig.add_trace(go.Scatter(x=tabela1[0], y=[0.5] * tamanho,
                             mode='lines', line=dict(color='green'), name='Limiar'))

    fig.update_layout(
        title="Propensão à agregação/" + CodigoPDB1 + "-" + chain1,
        xaxis=dict(title='Resíduo', showgrid=True,
                   zeroline=True, showline=True),
        yaxis=dict(title='Probabilidade', showgrid=True,
                   zeroline=True, showline=True),
        showlegend=True,
        hovermode='closest'
    )

    return fig

def main():
    st.title("Gerador de Gráfico de Agregação Dinâmico")
    st.markdown(
        """
        <p style='text-align: justify;'>O Gerador de Gráfico de Agregação Dinâmico é uma ferramenta poderosa usada em ciências biológicas, farmacêuticas e de pesquisa para visualizar e analisar dados relacionados à agregação de proteínas. 
        A agregação de proteínas é um fenômeno importante que pode levar a doenças neurodegenerativas, como Alzheimer, Parkinson e outras.
        </p>
        """,
        unsafe_allow_html=True,
    )

    st.subheader("Como Funciona")
    st.markdown(
    """
    <p style='text-align: justify;'>Os usuários começam inserindo seus dados no gerador, que pode incluir informações sobre código PDB. Em seguida, o gerador gera gráfico(s) que mostram como essas variáveis da agregação de proteínas acontecem.</p>
    """,
    unsafe_allow_html=True,
    )

    st.subheader("Visualização de Dados")
    st.markdown(
    """
    <p style='text-align: justify;'>Uma das características mais valiosas do Gerador de Gráfico de Agregação Dinâmico é sua capacidade de criar visualizações claras e informativas. 
    Os gráficos gerados podem incluir linhas de agregação de proteínas. Isso permite que os pesquisadores identifiquem tendências e padrões nos dados, o que é essencial para a tomada de decisões informadas.</p>
    """,
    unsafe_allow_html=True,
    )

    codigo_pdb = st.text_input("Informe o código PDB da proteína:")

    if st.button("Gerar Gráfico"):
        if not codigo_pdb:
            st.warning("Por favor, insira um código PDB válido.")
            return

        global chain  # Declarar chain como global para atualizá-lo

        arquivo = "C:/magre_dinamico/"
        arquivo1 = arquivo[0:17] + "/arquivos/pdb/"
        arquivo2 = arquivo[0:17] + "/arquivos/testes/"

        pdbr = 'https://files.rcsb.org/view/' + codigo_pdb + '.pdb'

        st.write("Código Solicitado:", codigo_pdb)

        response = requests.get(pdbr)

        entrada = saida = arquivo1 + codigo_pdb + ".pdb"
        A = response.text
        proteinArq = ""
        pdbLines = ""
        A3DLines = ""

        lin = 81
        fim = int(len(A) / lin)

        contpdb = 0
        desliga = 0
        proteinArq = codigo_pdb + ".pdb"
        chamador = 1
        resultado = []

        if response.status_code == 200:
            arq_saida = open(saida, 'w')
            for I in range(0, fim - 1, 1):
                words = A[lin * I:lin * (I + 1)]
                arq_saida.write(words)
            arq_saida.close()

        with open(entrada) as arq_entrada:
            pdbLines = arq_entrada.readlines()

        for line in pdbLines:
            words = line
            if len(words) < 1:
                pass
            elif (words[0:4] == "ATOM" and words[21:22] != chain):
                chain = words[21:22]

                resultado = Prog_Funcoes1.movimenta(
                    proteinArq, pdbLines, A3DLines, chain, codigo_pdb, chamador)

                resultado = Prog_Funcoes1.avalia_ruido(resultado)

                PastaArq = arquivo2
                contatohpFilename = PastaArq + codigo_pdb[0:4] + "_" + chain + ".csv"
                with open(contatohpFilename, "w") as contatohpFile:
                    for line in resultado:
                        print(line, file=contatohpFile)

                if (codigo_pdb, chain) in graficos_gerados:
                    st.warning("O(s) gráfico(s) já foi(ram) gerado(s) para esta cadeia do código PDB. Não é possível gerar novamente.")
                    return
                else:
                    # Exibe o gráfico mais informativo
                    fig = Plota_Resultado(resultado, codigo_pdb, chain)
                    st.plotly_chart(fig)

                    graficos_gerados.add((codigo_pdb, chain))

if __name__ == '__main__':
    main()
