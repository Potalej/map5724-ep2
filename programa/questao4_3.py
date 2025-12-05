"""
MAP5724 - EP 2: Metodo ADI para EDPs Parabolicas

Nome: Octavio Augusto Potalej
NUSP: 12558676

> QUESTAO 4.3
  O ponto central desta investigacao eh verificar se a convergencia numerica
  observada condiz com a consistencia teorica dos metodos.
  - Definam um tempo final T (e.g., T=0.5)
  - Para um refinamento N, definam dx = dy = h = pi/N
  - Sigam a diretriz dt = dx = dy, ou seja, k = h.
  - Calcule o erro global no tempo T usando a norma do maximo:
      E(h) = max|u_e(T,x,y) - v(T,x,y)|
  - Construam uma tabela para cada metodo (PR e DR) contendo h, o erro E(h),
    e a taxa de convergencia p (calculada por p = log_2(E(2h)/E(h))). Utilizem
    uma serie de refinamentos (e.g., N=10, 20, 40, 80).
"""
import helpers
import numpy as np
from time import time
from tabulate import tabulate
from questao4_2 import peaceman_rachford, douglas_rachford
import os

def questao_4_3_convergencia (Ms:list, saida_pr:str, saida_dr:str, dts_pr:list=[], dts_dr:list=[])->None:
  """
  Teste de convergencia dos metodos Peaceman-Rachford e Douglas-Rachford
  para o problema enunciado na questao 4. Fiz uma funcao para poder usa-la
  na questao seguinte.
  
  Parametros
  ----------
  Ms : list
    Lista com os valores de M (tamanhos da grade quadrada)
  dts_pr : list = []
    Tamanhos de passo temporal para o PR. Se passado vazio, dts = pi/Ms.
  dts_dr : list = []
    Tamanhos de passo temporal para o DR. Se passado vazio, dts = pi/Ms.
  """
  # Diretorio para salvar as tabelas em LaTeX
  if not os.path.isdir('latex'): os.mkdir('latex')
  
  # Parametros do problema
  b1, b2 = 1.0, 1.0   # Parametros da equacao
  eixos = [0, np.pi]  # Intervalo da grade
  T = 0.5             # Tempo final
  fronteira = 0.0     # Condicao de contorno
  checkpoints = 8
  
  # Valor inicial
  u0 = lambda z: np.sin(z[0])*np.sin(z[1])
  # Solucao
  solucao = lambda t, sol0: np.exp(-2*t)*sol0
  
  # Para salvar os erros
  erros_pr = []
  erros_dr = []
  for i, M in enumerate(Ms):
    print(f"Rodando para M={M}...", end=" ")
    timer = time()
    
    # Tamanhos de passo espacial e temporal
    h = np.pi/M
    dt_pr = h if len(dts_pr) == 0 else dts_pr[i]
    dt_dr = h if len(dts_dr) == 0 else dts_dr[i]

    # Valores iniciais
    x = np.linspace(0, np.pi, M+1)
    grid = np.meshgrid(x,x)
    v0 = u0(grid)

    """
    Primeiro teste: Peaceman-Rachford
    """
    metodo = peaceman_rachford
    sol_pr = helpers.rodar(v0, dt_pr, h, h, b1, b2, metodo, fronteira, T, checkpoints)

    tf = max(sol_pr.keys())
    erro_ultimo = np.max(np.abs(sol_pr[tf] - solucao(tf, v0)))
    erros_pr.append(erro_ultimo)

    """
    Segundo teste: Douglas-Rachford
    """
    metodo = douglas_rachford
    sol_dr = helpers.rodar(v0, dt_dr, h, h, b1, b2, metodo, fronteira, T, checkpoints)

    tf = max(sol_dr.keys())
    erro_ultimo = np.max(np.abs(sol_dr[tf] - solucao(tf, v0)))
    erros_dr.append(erro_ultimo)
    
    print(f"({(time() - timer):.2f}s)")
    
  # P-estimativa de convergencia
  ps_erros_pr = [np.log2(erros_pr[i]/erros_pr[i+1]) for i in range(len(erros_pr)-1)]
  ps_erros_dr = [np.log2(erros_dr[i]/erros_dr[i+1]) for i in range(len(erros_dr)-1)]

  # Montando a tabela para exibicao
  hs = [np.pi/M for M in Ms]
  cabecalho = ["$h$", "$E(h)$", "$p(h)$"]

  # Peaceman-Rachford
  tabela_pr = [
    [f"{h:.4f}" for h in hs], 
    [f"{e:.2e}" for e in erros_pr],
    ['-'] + [f"{p:.4f}" for p in ps_erros_pr]
  ]
  tabela_pr = list(zip(*tabela_pr))
  print("\n", tabulate(tabela_pr, cabecalho))

  # Douglas-Rachford
  tabela_dr = [
    [f"{h:.4f}" for h in hs],
    [f"{e:.2e}" for e in erros_dr],
    ['-'] + [f"{p:.4f}" for p in ps_erros_dr]
  ]
  tabela_dr = list(zip(*tabela_dr))
  print("\n", tabulate(tabela_dr, cabecalho))

  # Salvando em LaTeX para colocar no relatorio
  with open("latex/" + saida_pr, 'w') as arquivo:
    arquivo.write(tabulate(tabela_pr, cabecalho, colalign=("left","left","left"), tablefmt="latex_raw"))
  with open("latex/" + saida_dr, 'w') as arquivo:
    arquivo.write(tabulate(tabela_dr, cabecalho, colalign=("left","left","left"), tablefmt="latex_raw"))

"""
============================================================
TESTES PEDIDOS COM k=h

As tabelas sao exibidas no terminal e salvas em LaTeX na
pasta "latex".
============================================================
"""
if __name__ == '__main__':
  """
  Vamos usar dt = h para ambos os metodos, e os valores de M = 10, 20, 40, 80.
  """
  expoentes = range(0,4)
  Ms = [10 * (2**i) for i in expoentes]
  dts_pr = []
  dts_dr = []
  questao_4_3_convergencia(Ms, 'tabela_pr.txt', 'tabela_dr.txt', dts_pr, dts_dr)