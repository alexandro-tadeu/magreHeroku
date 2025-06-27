import streamlit as st
from stmol import showmol
import py3Dmol

def main():
    st.sidebar.title('Apresentação 3D da proteina')

    # Adicionando uma caixa de texto para o código PDB
    pdb_code = st.sidebar.text_input('Insira o código PDB e pressione Enter:')

    # Seletor de estilo
    style = st.sidebar.selectbox('Selecione o estilo de visualização', ['cartoon', 'stick', 'sphere'])

    # Verificando se um código PDB foi fornecido
    if pdb_code:
        # Função para renderizar a molécula
        def render_mol(pdb, style):
            xyzview = py3Dmol.view(query=f'pdb:{pdb}')
            xyzview.setStyle({style: {'color': 'spectrum'}})
            xyzview.setBackgroundColor('white')
            showmol(xyzview, height=500, width=800)

        # Chamando a função para renderizar a molécula com o código PDB inserido e o estilo selecionado
        render_mol(pdb_code, style)
    else:
        st.warning("Por favor, insira um código PDB na caixa de texto acima e pressione Enter.")
