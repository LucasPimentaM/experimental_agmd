import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, fsolve

# ==========================================
# 1. DEPENDÊNCIAS E PROPRIEDADES
# ==========================================
def rho(temperature, salinity):
    a = [9.999e2, 2.034e-2, -6.162e-3, 2.261e-5, -4.657e-8]
    b = [8.020e2, -2.001, 1.677e-2, -3.060e-5, -1.613e-5]

    temperature_part = a[0] + a[1] * temperature
    temperature_part += a[2] * temperature * temperature
    temperature_part += a[3] * temperature * temperature * temperature
    temperature_part += a[4] * temperature * temperature * temperature * temperature

    salinity_part = b[0] * salinity
    salinity_part += b[1] * salinity * temperature
    salinity_part += b[2] * salinity * temperature * temperature
    salinity_part += b[3] * salinity * temperature * temperature * temperature
    salinity_part += b[4] * salinity * salinity * temperature * temperature

    return temperature_part + salinity_part

def salinity_equation(x, salinity, temperature):
    return rho(temperature, x) * x - salinity

# ==========================================
# 2. FUNÇÃO OBJETIVO (O PROBLEMA INVERSO)
# ==========================================
def funcao_residuos_vagmd(params, cases):
    """
    Roda o simulador vagmd0Dmodel para todos os casos experimentais
    utilizando o air_gap atual testado pelo Levenberg-Marquardt.
    """
    air_gap_estimado = params[0] # O air gap que o otimizador quer testar
    residuos = []
    
    print(f"\n---> [Otimizador] Testando Air Gap: {air_gap_estimado:.6f} m")

    for index, elem in cases.iterrows():
        # Extração de variáveis do caso atual
        hot_feedwater_inlet_temperature = elem['Feed temperature at the inlet (°C)']
        cold_feedwater_inlet_temperature = elem['Coolant temperature at the inlet (°C)']
        vacuum_pressure = elem['Vacuum pressure (Pa)']
        
        initial_guess = elem['Feed salinity at the inlet (g/L)'] / 1000.0
        feed_salinity = fsolve(salinity_equation, initial_guess, args=(elem['Feed salinity at the inlet (g/L)'], elem['Feed temperature at the inlet (°C)']))[0]
        
        initial_guess = elem['Coolant salinity at the inlet (g/L)'] / 1000.0
        cool_salinity = fsolve(salinity_equation, initial_guess, args=(elem['Coolant salinity at the inlet (g/L)'], elem['Coolant temperature at the inlet (°C)']))[0]
        
        BaCl2_concentration = elem['Barium concentration at the inlet (g/L)']
        feed_mass_flow_rate = rho(hot_feedwater_inlet_temperature, feed_salinity) * elem['Feed flow rate (L/h)'] / 3600000.0
        cool_mass_flow_rate = rho(cold_feedwater_inlet_temperature, cool_salinity) * elem['Coolant flow rate (L/h)'] / 3600000.0
        
        membrane_area = elem['Membrane area (m²)']
        membrane_thickness = elem['Membrane thickness (m)']
        membrane_porosity = elem['Membrane porosity']
        pore_diameter = elem['Pore diameter (m)']
        polymer_conductivity = elem['Polymer conductivity (W/mK)']
        feed_channel_height = elem['Feed channel height (m)']
        cold_channel_height = elem['Cold channel height (m)']
        channel_width = elem['Channel width (m)']
        spacer_porosity = elem['Spacer porosity']
        gap_spacer_porosity = elem['Gap spacer porosity']
        wall_thickness = elem['Wall thickness (m)']
        spacer_conductivity = elem['Spacer thermal conductivity (W/mK)']
        wall_conductivity = elem['Condensing wall thermal conductivity (W/mK)']
        number_channels = elem['Number of channels']
        
        # Fluxo experimental alvo[cite: 1]
        fluxo_experimental = elem['Distilled water flux (L/h)']

        # RODA O SIMULADOR EXECUTÁVEL (Note que o air_gap_thickness foi substituído pela variável do otimizador)
        comando = ('./bin/vagmd0Dmodel -entry_temperature_feed ' + str(hot_feedwater_inlet_temperature) +
                  ' -entry_temperature_cool ' + str(cold_feedwater_inlet_temperature) +
                  ' -feed_mass_flow_rate ' + str(feed_mass_flow_rate) +
                  ' -cool_mass_flow_rate ' + str(feed_mass_flow_rate) + # Note: mantive feed_mass_flow_rate aqui pois estava assim no seu código original
                  ' -vacuum_pressure ' + str(vacuum_pressure) +
                  ' -BaCl2_concentration ' + str(BaCl2_concentration) +
                  ' -entry_salinity_feed ' + str(feed_salinity) +
                  ' -entry_salinity_cool ' + str(cool_salinity) +
                  ' -membrane_area ' + str(membrane_area) +
                  ' -number_channels ' + str(number_channels) +
                  ' -wall_conductivity ' + str(wall_conductivity) +
                  ' -spacer_conductivity ' + str(spacer_conductivity) +
                  ' -wall_thickness ' + str(wall_thickness) +
                  ' -air_gap_thickness ' + str(air_gap_estimado) + # <--- PARÂMETRO A SER OTIMIZADO
                  ' -gap_spacer_porosity ' + str(gap_spacer_porosity) +
                  ' -spacer_porosity ' + str(spacer_porosity) +
                  ' -channel_width ' + str(channel_width) +
                  ' -feed_channel_height ' + str(feed_channel_height) +
                  ' -cold_channel_height ' + str(cold_channel_height) +
                  ' -polymer_conductivity ' + str(polymer_conductivity) +
                  ' -pore_diameter ' + str(pore_diameter) +
                  ' -membrane_porosity ' + str(membrane_porosity) +
                  ' -membrane_thickness ' + str(membrane_thickness))
        
        # Oculta a saída padrão do C++ no terminal para não poluir, ou deixe os.system(comando) se quiser ver
        os.system(comando + " > /dev/null 2>&1")

        # LÊ O RESULTADO DA SIMULAÇÃO[cite: 1]
        try:
            report = pd.read_csv('./results/report.csv', header=None)
            distillate_flow_rate_simulado = float(report[report[0] == 'Distillate flow rate ='][1].values[0])
        except Exception as e:
            print(f"Erro ao ler report.csv no caso {elem['Case']}: {e}")
            distillate_flow_rate_simulado = 0.0

        # CÁLCULO DO RESÍDUO
        residuo = distillate_flow_rate_simulado - fluxo_experimental
        residuos.append(residuo)
        
    # Calcula e exibe o RMSE atual para você acompanhar o avanço
    rmse_atual = np.sqrt(np.mean(np.array(residuos)**2))
    print(f"RMSE Atual: {rmse_atual:.6f} | Erros: {residuos[:3]}...") # Mostra amostra dos erros

    return np.array(residuos)


# ==========================================
# 3. ROTINA PRINCIPAL (OTIMIZAÇÃO)
# ==========================================
def estimar_air_gap():
    print("Carregando dados experimentais...")
    cases = pd.read_excel('Experimento_bancada.xlsx', 'experimentos_adaptado')
    
    # ATENÇÃO: Como rodar o C++ para os 42 casos a cada iteração pode ser demorado, 
    # recomendo descomentar a linha abaixo na primeira vez para testar com apenas 3 casos:
    #cases = cases.head(3)

    # Chute inicial em metros (ex: 2 mm = 0.002 m)
    chute_inicial = [0.002]
    
    print("Iniciando Levenberg-Marquardt...")
    
    resultado = least_squares(
        fun=funcao_residuos_vagmd, 
        x0=chute_inicial, 
        method='trf', 
        bounds=(0.0001, 0.05), # Limites físicos: 0.1 mm até 50 mm
        args=(cases,) 
    )

    if resultado.success:
        air_gap_otimo = resultado.x[0]
        print(f"\n" + "="*40)
        print(f"✅ OTIMIZAÇÃO CONCLUÍDA!")
        print(f"O valor estimado para o Air Gap é: {air_gap_otimo:.6f} m")
        print(f"Custo Final (Soma dos quadrados): {resultado.cost:.6f}")
        print("="*40 + "\n")
    else:
        print("\n❌ A otimização falhou.")
        print(resultado.message)

# Executa o script
if __name__ == '__main__':
    estimar_air_gap()