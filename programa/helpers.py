"""
MAP5724 - EP 2
Metodo ADI para EDPs Parabolicas

Nome: Octavio Augusto Potalej
NºUSP: 12558676

Funcoes de ajuda para rodar os metodos ate que t=T e para
fazer os plots de calor.
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def plotar (solucoes:dict, dt:float, dx:float, dy:float, b1:float, b2:float, 
                xlim:list=[0,1], ylim:list=[0,1], vmin=0.0, vmax=0.5, 
                titulo:str="", saida:str=""):
  """
  Plota uma mapa de calor para um dado dicionario {t: u(t)}.
  """
  # Diretorio para salvar imagens
  if not os.path.isdir('img'): os.mkdir('img')
  
  # Configuracoes da figura
  fig, ax = plt.subplots(1, len(solucoes), figsize=(20,3))
  sup = r" Equação do calor com $\Delta t={:.5f}$, $\Delta x={}$, $\Delta y={}$, $b_1={}$, $b_2={}$".format(dt, dx, dy, b1, b2)
  if len(titulo): titulo_sup = f"({titulo}) {sup}"
  else: titulo_sup = sup
  fig.suptitle(titulo_sup, fontsize=18)
  
  # Grade
  x = np.arange(xlim[0], xlim[1]+dx, dx)
  y = np.arange(ylim[0], ylim[1]+dy, dy)
  
  # Plots
  for i, t in enumerate(solucoes):
    axpcm = ax[i].pcolormesh(x, y, solucoes[t], cmap=plt.cm.jet, vmin=vmin, vmax=vmax)
    ax[i].set_title(f"t={t:.2f}")
    fig.colorbar(axpcm, ax=ax[i])
    
    if i > 0: ax[i].set_yticks([]) # Tirando os ticks verticais dos ax[i>0]
  
  # Ajuste visual para aparecer o titulo
  plt.tight_layout(rect=[0,0,1,1])
  
  # Salvando a imagem
  if saida == "":
    titulo = titulo.replace(" ", "_")
    plt.savefig(f"img/{titulo}.png")
  else:
    plt.savefig(f"img/{saida}.png")

def rodar (v0:np.array, dt:float, dx:float, dy:float, b1:float, b2:float, metodo, 
          fronteira:float=0.0, tf:float=1.0, qntd_checkpoints:int=10):
  """
  Essa funcao vai chamando o metodo de integracao ate que t=T.
  """
  # Tamanho do grid
  M,N = v0.shape

  # Integracao
  qntd_iteracoes = int(np.ceil(tf / dt))
  instantes = np.arange(0, tf+dt, dt)
  
  # Pontos que serao salvos
  pontos_checkpoints = np.linspace(0, len(instantes)-1, qntd_checkpoints).round().astype(int)
  inst_checks = [instantes[i] for i in pontos_checkpoints]

  # Dicionario de solucoes
  sols = dict()
  v = v0.copy()

  for t in range(qntd_iteracoes+1):
    # Salvando no dicionario
    if t*dt in inst_checks:
      sols[t*dt] = v.copy()
      
    # Integra
    v = metodo(v, dt, dx, dy, b1, b2)
    
  return sols