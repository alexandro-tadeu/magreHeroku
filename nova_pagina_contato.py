import streamlit as st
import smtplib
import ssl
from email.message import EmailMessage

# Substitua pelos seus dados
email_origem = "magreaufabc@gmail.com"
email_senha = "qjsm kqka zbap oixi"

def enviar_email(nome, email_usuario, mensagem_usuario):
    assunto = "Aplicativo magre web"

    # Corrigido: preparar mensagem com quebras de linha visíveis em HTML
    mensagem_formatada = mensagem_usuario.replace('\n', '<br>')

    corpo_html = f"""
    <html>
        <body>
            <h2>Mensagem recebida do formulário</h2>
            <p><strong>Nome:</strong> {nome}</p>
            <p><strong>Email:</strong> {email_usuario}</p>
            <p><strong>Mensagem:</strong><br>{mensagem_formatada}</p>
        </body>
    </html>
    """

    mensagem = EmailMessage()
    mensagem["From"] = email_origem
    mensagem["To"] = email_origem
    mensagem["Subject"] = assunto

    mensagem.set_content(f"Mensagem de {nome} <{email_usuario}>:\n\n{mensagem_usuario}")
    mensagem.add_alternative(corpo_html, subtype='html')

    contexto_seguro = ssl.create_default_context()
    with smtplib.SMTP_SSL('smtp.gmail.com', 465, context=contexto_seguro) as smtp:
        smtp.login(email_origem, email_senha)
        smtp.send_message(mensagem)

def main():
    st.title("Página de Contato")

    with st.form(key='myform', clear_on_submit=True):
        nome = st.text_input("Nome")
        email = st.text_input("Email")
        mensagem = st.text_area("Mensagem", height=300)
        enviar = st.form_submit_button("Enviar Mensagem")

        if enviar:
            if nome and email and mensagem:
                try:
                    enviar_email(nome, email, mensagem)
                    st.success("Mensagem enviada com sucesso!")
                except Exception as e:
                    st.error(f"Ocorreu um erro ao enviar o email: {e}")
            else:
                st.warning("Por favor, preencha todos os campos.")

if __name__ == "__main__":
    main()
