import Prog_Mod  # Certifique-se de importar corretamente o módulo Prog_Mod

def verifica_agregacao_sequencia(sequencia):
    janela = 1
    tabela3 = []

    for i in range(len(sequencia) - janela):
        amino1 = sequencia[i]
        amino2 = sequencia[i + 1]

        # Exemplo de lógica: verificar se amino1 e amino2 são propensos à agregação
        propensao_agregacao_amino1 = Prog_Mod.verifica_propensao_agregacao(amino1)
        propensao_agregacao_amino2 = Prog_Mod.verifica_propensao_agregacao(amino2)

        # Lógica de verificação de agregação (ajuste conforme necessário)
        if propensao_agregacao_amino1 and propensao_agregacao_amino2:
            tabela3.append(1)  # Indica possível agregação entre amino1 e amino2
        else:
            tabela3.append(0)  # Não há indicação de agregação entre amino1 e amino2

    return tabela3

# Exemplo de uso:
sequencia_aminoacidos = "ARNDC"
tabela_agregacao = verifica_agregacao_sequencia(sequencia_aminoacidos)

# Aqui, a tabela_agregacao contém 1 ou 0 indicando a agregação entre aminoácidos consecutivos
print(tabela_agregacao)
