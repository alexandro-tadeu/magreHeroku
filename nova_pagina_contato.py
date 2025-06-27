# nova_pagina_contato.py

import streamlit as st

def main():
    st.title("Página de Contato")

    # Adicione campos de entrada
    nome = st.text_input("Nome")
    email = st.text_input("Email")
    mensagem = st.text_area("Mensagem")

    # Exemplo de como usar os valores inseridos nos campos
    if st.button("Enviar Mensagem"):
        st.write(f"Nome: {nome}")
        st.write(f"Email: {email}")
        st.write(f"Mensagem: {mensagem}")

if __name__ == "__main__":
    main()
