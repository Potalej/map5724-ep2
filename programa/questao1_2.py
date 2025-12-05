"""
MAP5724 - EP 2
Metodo ADI para EDPs Parabolicas

Nome: Octavio Augusto Potalej
NºUSP: 12558676

> QUESTAO 1.2
  Implemente computacionalmente o esquema explicito para resolver a
  Equacao 1. Defina um dominio computacional, escolha valores
  adequados para b1 e b2 e estabeleca condicoes iniciais u(0,x,y)=f(x,y)
  e condicoes de fronteira (e.g., Dirichlet homogeneas) apropriadas.
  Verifique empiricamente a condicao de estabilidade obtida na questao 1.1.
"""
import numpy as np
import helpers

"""
============================================================
IMPLEMENTACAO DO METODO FTCS
============================================================
"""
# Operadores de diferencas centradas de segunda ordem
# A condicao de fronteira eh = 0
def A_x (v:np.array, delta_x:float)->np.array:
	vet = np.zeros_like(v)
	vet[1:-1] = v[2:] - 2 * v[1:-1] + v[:-2]
	return vet / (delta_x*delta_x)

def A_y (v:np.array, delta_y:float)->np.array:
	vet = np.zeros_like(v)
	vet[1:-1,1:-1] = v[1:-1, 2:] - 2 * v[1:-1, 1:-1] + v[1:-1, :-2]
	return vet / (delta_y*delta_y)

# Metodo FTCS
def ftcs (v0:np.array, delta_t:float, delta_x:float, delta_y:float, b_1:float, b_2:float, fronteira=0.0)->np.array:
  # Aplicando os operadores de diferencas centradas
  Ax_v0 = A_x(v0, delta_x)[1:-1,1:-1]
  Ay_v0 = A_y(v0, delta_y)[1:-1,1:-1]

  # Condicao de fronteira
  v = np.zeros_like(v0) + fronteira
  # Metodo ftcs
  v[1:-1,1:-1] = v0[1:-1,1:-1] + delta_t * (b_1 * Ax_v0 + b_2 * Ay_v0)

  return v

"""
============================================================
VERIFICACAO DA CONDICAO DE ESTABILIDADE

Os mapas de calor sao salvos na pasta "img".
============================================================
"""
# Parametros do problema
b1, b2 = .25, .25     # Parametros da equacao
dx, dy = 1/50, 1/50   # Tamanho dos intervalos
T = 0.4               # Tempo final
metodo = ftcs         # Metodo
eixos = [0, 1]        # Intervalo da grade
fronteira = 0.0       # Condicao de contorno
checkpoints = 8

# Condicao inicial: o grid inteiro = 1
M,N = int(1/dx), int(1/dy)
v0 = np.zeros((M,N)) + 1.0

# Calculo o delta_t conforme a condicao de estabilidade
dx2, dy2 = dx*dx, dy*dy
dt_estavel = dx2 * dy2 / (2*b1*dy2 + 2*b2*dx2)
print(f"Rodando FTCS Estavel (dt={dt_estavel:.6f})...", end=" ")
sols = helpers.rodar(v0, dt_estavel, dx, dy, b1, b2, metodo, fronteira, T, checkpoints)
helpers.plotar(sols, dt_estavel, dx, dy, b1, b2, eixos, eixos, 0, 1, "FTCS Estável", "q1.2_estavel")
print("Imagem salva: img/q1.2_estavel.png")

# Usando um delta_t ligeiramente maior
dt_instavel = dt_estavel * 1.02
print(f"Rodando FTCS Instavel (dt={dt_instavel:.6f})...", end=" ")
sols = helpers.rodar(v0, dt_instavel, dx, dy, b1, b2, metodo, fronteira, T, checkpoints)
helpers.plotar(sols, dt_instavel, dx, dy, b1, b2, eixos, eixos, 0, 1, "FTCS Instável", "q1.2_instavel")
print("Imagem salva: img/q1.2_instavel.png")