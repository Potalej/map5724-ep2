"""
MAP5724 - EP 2
Metodo ADI para EDPs Parabolicas

Nome: Octavio Augusto Potalej
NºUSP: 12558676

> QUESTAO 4.2
  Implementem computacionalmente ambos os esquemas (PR) e (DR) para
  resolver o problema modelo definido acima. Utilizem um algoritmo
  eficiente (e.g, Algoritmo de Thomas) para a inversao dos sistemas
  tridiagonais em cada estagio.
  Gerem uma visualizacao (e.g., um mapa de calor ou grafico de
  superficie) da solucao numerica v_{i,j}^n e do erro absluto
  |u_e - v| em um tempo t = T (ex., T=0.5).
"""
import numpy as np
from numba import jit
import helpers

"""
============================================================
IMPLEMENTACAO DO ALGORITMO DE THOMAS DOS METODOS ADI

Para acelerar o processo, usei o Numba.
============================================================
"""
@jit
def thomas (a:float, b:float, c:float, d:np.array)->np.array:
  """
  Algoritmo de Thomas para resolucao de sistemas tridiagonais.
  """
  x = np.zeros_like(d)
  n = len(d)
  x[-1] = d[-1]

  c_ = np.zeros(n-1)
  c_[0] = c / b
  d[0] = d[0] / b

  for i in range(1, n):
    if i < n - 1: c_[i] = c / (b - c_[i-1] * a)
    d[i] = (d[i] - d[i-1] * a) / (b - c_[i-1] * a)

  x[-1] = d[-1]
  for i in range(n-2, -1, -1):
    x[i] = d[i] - c_[i] * x[i+1]

  return x

@jit
def peaceman_rachford (v:np.array, delta_t:float, delta_x:float, delta_y:float, b_1:float, b_2:float)->np.array:
  """
  Integracao numerica da equacao do calor usando o metodo Peaceman-Rachford.
  """
  M, N = v.shape

  # Constantes
  mu_x = delta_t / (delta_x**2)
  mu_y = delta_t / (delta_y**2)

  # Vetor de estados intermediario
  v_tilde = np.zeros((M,N))

  # Primeiro fazemos a integracao em x a partir de v
  # Termos das diagonais principal e secundarias
  diag_pri = 1. + b_1 * mu_x
  diag_sec = -.5 * b_1 * mu_x

  # Integracao
  for b in range(1, N-1):
    lado_direito = v[1:-1,b] + 0.5 * b_2 * mu_y * (v[1:-1,b-1] - 2 * v[1:-1,b] + v[1:-1,b+1])
    v_tilde[1:-1,b] = thomas(diag_sec, diag_pri, diag_sec, lado_direito)
  
  # Agora fazemos a integracao em y a partir de v_tilde
  # Termos das diagonais principal e secundarias
  diag_pri = 1. + b_2 * mu_y
  diag_sec = -.5 * b_2 * mu_y

  # Integracao
  for a in range(1, M-1):
    lado_direito = v_tilde[a,1:-1] + 0.5 * b_1 * mu_x * (v_tilde[a-1,1:-1] - 2 * v_tilde[a,1:-1] + v_tilde[a+1,1:-1])
    v[a,1:-1] = thomas(diag_sec, diag_pri, diag_sec, lado_direito)
  
  return v

@jit
def douglas_rachford (v:np.array, delta_t:float, delta_x:float, delta_y:float, b_1:float, b_2:float)->np.array:
  """
  Integracao numerica da equacao do calor usando o metodo Douglas-Rachford.
  """
  M, N = v.shape

  # Constantes
  mu_x = delta_t / (delta_x**2)
  mu_y = delta_t / (delta_y**2)

  # Vetor de estados intermediario
  v_tilde = np.zeros((M,N))

  # Primeiro fazemos a integracao em x a partir de v
  # Termos das diagonais principal e secundarias
  diag_pri = 1. + 2. * b_1 * mu_x
  diag_sec = - b_1 * mu_x
  # Integracao
  for b in range(1, N-1):
    lado_direito = v[1:-1,b] + b_2 * mu_y * (v[1:-1,b-1] - 2 * v[1:-1,b] + v[1:-1,b+1])
    v_tilde[1:-1,b] = thomas(diag_sec, diag_pri, diag_sec, lado_direito)
  
  # Agora fazemos a integracao em y a partir de v_tilde
  # Termos das diagonais principal e secundarias
  diag_pri = 1. + 2. * b_2 * mu_y
  diag_sec = -b_2 * mu_y
  v_prox = np.zeros((M,N))

  # Integracao
  for a in range(1, M-1):
    lado_direito = v_tilde[a,1:-1] - b_2 * mu_y * (v[a,0:-2] - 2.*v[a,1:-1] + v[a,2:])
    v_prox[a, 1:-1] = thomas(diag_sec, diag_pri, diag_sec, lado_direito)
  
  return v_prox

"""
============================================================
VISUALIZACAO DAS RESOLUCOES E DOS ERROS NUMERICOS

Escolhi usar uma grade 50 x 50 e um tamanho de passo
temporal dt = 0.005, apenas por ser parecido com o utilizado
na questao 1.2.

Os mapas de calor sao salvos na pasta "img".
============================================================
"""
if __name__ == '__main__':
  
  # Parametros do problema
  b1, b2 = 1.0, 1.0   # Parametros da equacao
  M, N = 50, 50       # Tamanho da grade
  eixos = [0, np.pi]  # Intervalo da grade
  dt = 0.005          # Tamanho de passo
  T = 0.5             # Tempo final
  fronteira = 0.0     # Condicao de contorno
  checkpoints = 8
  dx, dy = (eixos[1]-eixos[0])/M, (eixos[1] - eixos[0])/N

  # Valores iniciais
  u0 = lambda z: np.sin(z[0])*np.sin(z[1])
  x = np.linspace(0, np.pi, M+1)
  y = np.linspace(0, np.pi, N+1)
  grid = np.meshgrid(x,y)
  v0 = u0(grid)

  # Solucao do problema
  solucao = lambda t: np.exp(-2*t)*v0

  """
  Primeiro teste: Peaceman-Rachford
  """
  metodo = peaceman_rachford
  print("Rodando Peaceman-Rachford...", end=" ")
  sols_pr = helpers.rodar(v0, dt, dx, dy, b1, b2, metodo, fronteira, T, checkpoints)
  # Plota o mapa de calor da solucao
  helpers.plotar(sols_pr, dt, dx, dy, b1, b2, eixos, eixos, 0, 1, "PR", "q4.2_PR")
  print("Imagens salvas: img/q4.2_PR.png, ", end="")

  # Agora plotando o erro absoluto
  erro_absoluto = {t: np.abs(solucao(t) - sols_pr[t]) for t in sols_pr}
  helpers.plotar(erro_absoluto, dt, dx, dy, b1, b2, eixos, eixos, 0, 1e-4, "PR - Erro absoluto", "q4.2_PR_erro")
  print("img/q4.2_PR_erro.png")


  """
  Segundo teste: Douglas-Rachford
  """
  metodo = douglas_rachford
  print("Rodando Douglas-Rachford...", end=" ")
  sols_dr = helpers.rodar(v0, dt, dx, dy, b1, b2, metodo, fronteira, T, checkpoints)
  helpers.plotar(sols_dr, dt, dx, dy, b1, b2, eixos, eixos, 0, 1, "DR", "q4.2_DR")
  print("Imagens salvas: img/q4.2_DR.png, ", end="")

  # Agora plotando o erro absoluto
  erro_absoluto = {t: np.abs(solucao(t) - sols_dr[t]) for t in sols_dr}
  helpers.plotar(erro_absoluto, dt, dx, dy, b1, b2, eixos, eixos, 0, 0.003, "DR - Erro absoluto", "q4.2_DR_erro")
  print("img/q4.2_DR_erro.png")