import streamlit as st
from stmol import showmol
import py3Dmol
import requests

def obter_arquivo_pdb(codigo_pdb):
    url = f"https://files.rcsb.org/view/{codigo_pdb}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        st.error("Erro ao baixar o arquivo PDB. Verifique o código.")
        return None

def identificar_cadeias_pdb(pdb_texto):
    cadeias = set()
    for linha in pdb_texto.splitlines():
        if linha.startswith("ATOM") or linha.startswith("HETATM"):
            if len(linha) > 21:
                cadeia = linha[21]
                cadeias.add(cadeia)
    return sorted(list(cadeias))

def render_cadeia(pdb_texto, cadeia, style):
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_texto, 'pdb')
    view.setStyle({})  # Limpa todos os estilos anteriores (remove os fios vermelhos/cinza)
    view.setStyle({ "chain": cadeia }, { style: { "color": "spectrum" } })
    view.setBackgroundColor('white')
    view.zoomTo({ "chain": cadeia })
    showmol(view, height=500, width=800)

def main():
    st.sidebar.title('Apresentação 3D da proteína')

    # Input de código PDB
    codigo_pdb = st.sidebar.text_input('Insira o código PDB e pressione Enter:').upper().strip()

    # Estilo de visualização
    estilo = st.sidebar.selectbox('Selecione o estilo de visualização', ['cartoon', 'stick', 'sphere'])

    if codigo_pdb:
        pdb_texto = obter_arquivo_pdb(codigo_pdb)

        if pdb_texto:
            cadeias = identificar_cadeias_pdb(pdb_texto)

            if len(cadeias) == 1:
                st.write(f"Visualizando a única cadeia: {cadeias[0]}")
                render_cadeia(pdb_texto, cadeias[0], estilo)

            elif len(cadeias) > 1:
                cadeias_selecionadas = st.sidebar.multiselect("Selecione as cadeias a visualizar:", cadeias)
                if cadeias_selecionadas:
                    for cadeia in cadeias_selecionadas:
                        st.subheader(f"Visualização da cadeia {cadeia}")
                        render_cadeia(pdb_texto, cadeia, estilo)
                else:
                    st.info("Selecione ao menos uma cadeia para visualização.")
            else:
                st.warning("Nenhuma cadeia detectada no arquivo PDB.")
    else:
        st.warning("Por favor, insira um código PDB na caixa de texto.")

if __name__ == "__main__":
    main()
