import streamlit as st

def main():
    # Título da página
    st.title("Gerador de Sequência de Aminoácidos")

    st.markdown(
        """
        <p style='text-align: justify;'>
        O Gerador de Sequência de Aminoácidos é uma ferramenta para visualizar e identificar regiões de agregação em proteínas.
        </p>
        """,
        unsafe_allow_html=True,
    )

    st.subheader("Como Funciona")
    st.markdown(
    """
    <p style='text-align: justify;'>
    Digite a sequência de aminoácidos e, em seguida, digite os valores de propensão de agregação (um valor para cada aminoácido, separados por espaço, vírgula ou ponto e vírgula). A ferramenta destacará onde cada aminoácido agrega ou não.
    </p>
    """,
    unsafe_allow_html=True,
    )

    # Entrada de texto para a sequência de aminoácidos
    amino_seq = st.text_area("Digite a sequência de aminoácidos (apenas letras):").strip()
    propensao_text = st.text_area("Digite os valores de propensão de agregação para cada aminoácido, separados por espaço, vírgula ou ponto e vírgula:").strip()

    # Verificação se as entradas estão preenchidas
    if amino_seq and propensao_text:
        try:
            # Convertendo a sequência de propensão em uma lista de floats
            delimiters = [',', ';', ' ']
            for delimiter in delimiters:
                if delimiter in propensao_text:
                    propensao_vals = list(map(float, propensao_text.split(delimiter)))
                    break
            else:
                propensao_vals = list(map(float, propensao_text.split()))  # Caso padrão: espaço
            
            # Validando o comprimento da sequência e dos valores de agregação
            if len(amino_seq) != len(propensao_vals):
                st.error("O número de aminoácidos e de valores de agregação não corresponde. Verifique e tente novamente.")
                return
            
            # Variáveis para contagem
            agregam_count = 0
            nao_agregam_count = 0

            # Criando elementos HTML para a sequência
            sequence_elements = []
            for index, (amino, val) in enumerate(zip(amino_seq, propensao_vals)):
                color = "red" if val > 0.5 else "black"
                tooltip = "Agrega" if val > 0.5 else "Não Agrega"

                span_style = f"color: {color}; font-size: 20px; background-color: #F0F0F0; border: 1px solid #CCCCCC; border-radius: 4px; padding: 2px 6px; margin: 2px"
                sequence_elements.append(
                    f'<span style="{span_style}" title="{tooltip}">{amino}</span>'
                )

                # Contagem
                if val > 0.5:
                    agregam_count += 1
                else:
                    nao_agregam_count += 1

            # Exibindo sequência e resumo
            st.header("Sequência de Aminoácidos")
            sequence_html = " ".join(sequence_elements)
            st.markdown(sequence_html, unsafe_allow_html=True)

            st.header("Resumo")
            st.write(f"Número de aminoácidos que agregam: {agregam_count}")
            st.write(f"Número de aminoácidos que não agregam: {nao_agregam_count}")

        except ValueError:
            st.error("Erro: Certifique-se de que os valores de propensão de agregação estão no formato numérico.")
    
if __name__ == "__main__":
    main()
