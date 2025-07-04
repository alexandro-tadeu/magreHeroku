import streamlit as st
from nova_pagina import (
    Plota_Resultado,
    #main as NovaPaginaMain,
)  
from nova_pagina_dinamica import (
    main as NovaPaginaDinamicaMain,
)  
from nova_pagina_sequencia import main as NovaPaginaSequenciaMain
from nova_pagina_estrutura import main as NovaPaginaEstruturaMain  
from nova_pagina_contato import main as NovaPaginaContatoMain
from appSeqPdb import main as AppPdbMain
#from appSeqAmypro import main as AppAmyproMain
#from appTransforme import main as AppTransformeMain
#from appEstrutura3D import main as AppEstrutura3DMain
#from appMagreDois import main as AppMagreDoisMain
from appMagreTres import main as AppMagreTresMain  # Importação da nova página

# Define um menu de navegação
st.sidebar.title("Menu de Páginas")
pages = {
    "Inicial": "inicial",
    "Magre II Final": "grafico_dinamico",
    "Magre III Final": "appMagreTres",  # Nova opção adicionada ao menu
    #"Gráfico Estático": "grafico_estatico",
    "Sequência Aminoácidos": "nova_pagina_sequencia",
    "Obter Sequência Aminoácidos PDB": "appSeqPdb",
    #"Sequência Aminoácidos Amypro": "appSeqAmypro",
    #"Sequência Linear Aminoácido": "appTransforme",
    #"Estrutura 3D com Rotulagem": "appEstrutura3D",
    "Estrutura 3D sem Rotulagem": "nova_pagina_estrutura",
    #"Magre III WebServer": "appMagreDois",
    "Contato": "contato",
}
page = st.sidebar.selectbox("Ir para", list(pages.keys()))

# Função para criar um card
def create_card(title, text, image_url):
    st.subheader(title)
    st.image(image_url, use_container_width=True)
    st.write(text)

# Lista de imagens com links individuais (corrigido "url" para "image_url")
imagens = [
    {"title": "Card 1", "text": "Este é um subtexto explicativo para o Card 1.", "image_url": "img/amino.png"},
    {"title": "Card 2", "text": "Este é um subtexto explicativo para o Card 2.", "image_url": "img/beta.png"},
    {"title": "Card 3", "text": "Este é um subtexto explicativo para o Card 3.", "image_url": "img/helice.png"},
    {"title": "Card 4", "text": "Este é um subtexto explicativo para o Card 4.", "image_url": "img/helice.png"},
    {"title": "Card 5", "text": "Este é um subtexto explicativo para o Card 5.", "image_url": "img/helice.png"},
    {"title": "Card 6", "text": "Este é um subtexto explicativo para o Card 6.", "image_url": "img/helice.png"},
]

# Crie o aplicativo Streamlit
def main():
    if page == "Inicial":
        st.title("Proteínas e Agregação de Proteínas: Uso de Preditores")

        st.markdown(
            """
            <p style='text-align: justify;'>As proteínas desempenham um papel fundamental em inúmeras funções biológicas...</p>
            """,
            unsafe_allow_html=True,
        )

        col1, col2, col3 = st.columns(3)

        with col1:
            create_card(**imagens[0])
        with col2:
            create_card(**imagens[1])
        with col3:
            create_card(**imagens[2])

        with col1:
            create_card(**imagens[3])
        with col2:
            create_card(**imagens[4])
        with col3:
            create_card(**imagens[5])

    # elif page == "Gráfico Estático":
    #     NovaPaginaMain()

    elif page == "Magre II Final":
        NovaPaginaDinamicaMain()

    elif page == "Sequência Aminoácidos":
        NovaPaginaSequenciaMain()

    elif page == "Obter Sequência Aminoácidos PDB":
        AppPdbMain()

    # elif page == "Sequência Aminoácidos Amypro":
    #     AppAmyproMain()

    # elif page == "Sequência Linear Aminoácido":
    #     AppTransformeMain()

    # elif page == "Estrutura 3D com Rotulagem":
    #     AppEstrutura3DMain()

    elif page == "Estrutura 3D sem Rotulagem":
        NovaPaginaEstruturaMain()

    # elif page == "Magre III WebServer":
    #     AppMagreDoisMain()

    elif page == "Magre III Final":
        AppMagreTresMain()  # Chamada da nova página

    elif page == "Contato":
        NovaPaginaContatoMain()

if __name__ == "__main__":
    main()
