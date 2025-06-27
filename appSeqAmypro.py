import requests
import streamlit as st

def obter_sequencia_amypro(codigo_entry):
    # URL do arquivo .fasta
    url = f"https://www.bio2byte.be/amypro/data/entries/{codigo_entry}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        st.error(f"Erro ao obter dados do AmyPro. Status: {response.status_code}")
        return None

def extrair_informacoes(sequencia):
    # Separar as informações e a sequência de aminoácidos
    partes = sequencia.split("\n")
    # Separar a primeira linha com informações
    informacoes = partes[0].strip()
    # Inicializando um dicionário para armazenar dados
    dados = {}
    # A primeira parte contém a sequência de aminoácidos
    sequencia_aminoacidos = partes[1].strip() if len(partes) > 1 else ""
    # Processando a linha de informações
    for info in informacoes.split(" "):
        if "=" in info:
            chave, valor = info.split("=")
            dados[chave] = valor
        elif info.startswith(">"):
            dados["codigo"] = info[1:]  # Extrai o código da entrada
            # Extrai o nome da proteína e o organismo
            dados["nome_proteina"] = informacoes.split(" ")[1]  # O segundo elemento é o nome da proteína
        else:
            dados[info] = True  # Para informações que não têm um valor específico
    # A primeira parte da linha de informações é o nome da proteína
    # Extraindo regiões
    if "regions" in informacoes:
        start = informacoes.index("regions") + len("regions: ")
        end = informacoes.index("}", start) if "}" in informacoes[start:] else None
        dados["regioes"] = informacoes[start:end].strip().replace(" ", "").split(",") if end else []
    # Adicionando a sequência de aminoácidos ao dicionário
    dados["sequencia_aminoacidos"] = sequencia_aminoacidos
    # Extraindo o código PDB
    dados["codigo_pdb"] = dados.get("pdb", "N/A")
    return dados

def main():
    st.title("Extrator de Sequência de Aminoácidos do AmyPro")
    codigo_entry = st.text_input("Insira o código da entrada (por exemplo, AP00015):", "").upper()
    if st.button("Buscar Sequência"):
        if codigo_entry:
            sequencia = obter_sequencia_amypro(codigo_entry)
            if sequencia:
                informacoes_extraidas = extrair_informacoes(sequencia)
                # Exibir apenas as informações desejadas
                st.write(f"**Código:** {informacoes_extraidas['codigo']}")
                st.write(f"**Nome da Proteína:** {informacoes_extraidas.get('nome_proteina', 'N/A')}")
                st.write(f"**Código PDB:** {informacoes_extraidas.get('codigo_pdb', 'N/A')}")
                # # Exibir regiões de agregação formatadas
                # regions = informacoes_extraidas.get("regioes", [])
                # if regions:
                #     regions_formatadas = ", ".join(regions)  # Formata as regiões separadas por vírgulas
                #     st.write(f"**Regiões de Agregação:** {regions_formatadas}")
                # else:
                #     st.write("**Regiões de Agregação:** N/A")
                # Mostrando a sequência de aminoácidos separadamente
                st.success("Sequência de Aminoácidos:")
                st.code(informacoes_extraidas["sequencia_aminoacidos"])
            else:
                st.warning("Nenhuma sequência encontrada.")
        else:
            st.warning("Por favor, insira um código de entrada.")

if __name__ == "__main__":
    main()
