import requests
import streamlit as st

def obter_codigo_fasta(codigo_pdb):
    # URL para o formato FASTA do PDB
    url = f'https://www.rcsb.org/fasta/entry/{codigo_pdb}'
    response = requests.get(url)

    if response.status_code == 200:
        return response.text.splitlines()
    else:
        st.error(f"Erro ao obter o código FASTA {codigo_pdb}. Verifique o código e tente novamente.")
        return None

def extrair_sequencia_fasta(fasta_lines):
    sequencia = ""
    cadeia_info = ""
    for line in fasta_lines:
        # Ignora a linha de sequência e captura a linha de descrição
        if line.startswith(">"):
            cadeia_info = line.strip()  # Captura a linha de descrição
        else:
            sequencia += line.strip()  # Adiciona a linha de sequência sem espaços em branco
    return sequencia, cadeia_info

def obter_cadeia(cadeia_info):
    # Extrai as cadeias da linha de descrição
    if "|" in cadeia_info:
        partes = cadeia_info.split("|")
        if len(partes) > 1:
            return partes[1].replace("Chains ", "").strip()  # Remove "Chains " e espaços
    return "N/A"

# Função principal da aplicação Streamlit
def main():
    st.title("Extrator de Sequência de Aminoácidos do PDB")

    # Input do usuário para o código PDB
    codigo_pdb = st.text_input("Insira o código da entrada (por exemplo, 1XQ8):", "").upper()

    if st.button("Buscar Sequência"):
        if codigo_pdb:
            # Obter linhas do arquivo FASTA
            fasta_lines = obter_codigo_fasta(codigo_pdb)

            if fasta_lines:
                sequencia_aminoacidos, cadeia_info = extrair_sequencia_fasta(fasta_lines)

                if sequencia_aminoacidos:  # Verifica se a sequência não está vazia
                    # Obter a cadeia
                    cadeia = obter_cadeia(cadeia_info)
                    st.write(f"**Cadeia:** {cadeia}")
                    st.success("Sequência de Aminoácidos:")
                    st.code(sequencia_aminoacidos)
                    # st.write(sequencia_aminoacidos)
                else:
                    st.warning("Nenhuma sequência de aminoácidos encontrada.")
        else:
            st.warning("Por favor, insira um código PDB.")

if __name__ == "__main__":
    main()
