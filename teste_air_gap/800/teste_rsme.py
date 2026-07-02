import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def calculate_rmse(file_path):
    """Lê o arquivo CSV e calcula o RMSE entre os fluxos de destilado."""
    df = pd.read_csv(file_path)
    
    # Nomes das colunas conforme o seu arquivo de exemplo
    col_pred = 'Distillate flow rate (L/h)'
    col_exp = 'Exp. distillate flow rate (L/h)'
    
    # Verifica se as colunas existem no arquivo
    if col_pred in df.columns and col_exp in df.columns:
        # Cálculo do RMSE: Raiz da média dos erros ao quadrado
        rmse = np.sqrt(((df[col_pred] - df[col_exp]) ** 2).mean())
        return rmse
    else:
        print(f"Erro: Colunas não encontradas no arquivo {file_path}")
        return None

# Valores de x de 0.1 até 2.0
# Usamos round para evitar problemas de precisão de ponto flutuante no Python (ex: 0.30000000000000004)
x_values = [round(x, 1) for x in np.arange(0.1, 2.1, 0.1)]

rmse_list = []
valid_x = []

print("Calculando o RMSE para os arquivos...")

for x in x_values:
    # Formata o nome do arquivo para garantir 1 casa decimal (ex: results_0.1.csv)
    file_name = f'results_{x:.1f}.csv'
    
    # Verifica se o arquivo existe na pasta antes de tentar ler
    if os.path.exists(file_name):
        rmse = calculate_rmse(file_name)
        if rmse is not None:
            rmse_list.append(rmse)
            valid_x.append(x)
            print(f"Arquivo: {file_name} | RMSE: {rmse:.4f}")
    else:
        print(f"Aviso: Arquivo não encontrado -> {file_name}")

# Plotagem do gráfico geral se houverem dados processados
if valid_x and rmse_list:
    plt.figure(figsize=(10, 6))
    plt.plot(valid_x, rmse_list, marker='o', linestyle='-', color='b', label='RMSE do arquivo')
    plt.xlabel('Valor de x')
    plt.ylabel('RMSE - Distillate flow rate (L/h)')
    plt.title('Avaliação do RMSE para os arquivos results_x.x')
    plt.xticks(x_values) # Garante que todos os passos de 0.1 em 0.1 apareçam no eixo x
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Identifica e destaca o menor RMSE no gráfico
    min_rmse = min(rmse_list)
    min_index = rmse_list.index(min_rmse)
    min_x = valid_x[min_index]
    
    plt.plot(min_x, min_rmse, marker='*', markersize=15, color='r', label=f'Menor RMSE: {min_rmse:.4f} (x={min_x:.1f})')
    
    plt.legend()
    plt.tight_layout()
    
    # Salva o gráfico em uma imagem e também o exibe na tela
    plt.savefig('grafico_rmse_geral.png', dpi=300)
    plt.show()
    
    print(f"\n--- CONCLUSÃO ---")
    print(f"O arquivo que promove o menor RMSE é o 'results_{min_x:.1f}.csv' com um erro de {min_rmse:.4f}.")
else:
    print("\nNenhum dado válido para plotar. Verifique se os arquivos estão na mesma pasta que o script.")