"""
MAP5724 - EP 2: Metodo ADI para EDPs Parabolicas

Nome: Octavio Augusto Potalej
NUSP: 12558676

> QUESTAO 4.4 - Discussao
  Analisem os resultados. A taxa de convergencia p observada para k = h
  eh consistente com as ordens de truncamento teoricas? Discutam qualquer
  discrepancia ou confirmacao.
"""
from questao4_3 import questao_4_3_convergencia

"""
Vamos fixar o tamanho de passo de ambos...
"""
expoentes = range(0,4)
Ms = [10 * (2**i) for i in expoentes]
dt = 1e-4
dts_pr = [dt for _ in Ms]
dts_dr = dts_pr
questao_4_3_convergencia(Ms, 'tabela_pr_fixo.txt', 'tabela_dr_fixo.txt', dts_pr, dts_dr)

print()

"""
Por fim, vou usar h^2
"""
expoentes = range(0,4)
Ms = [10 * (2**i) for i in expoentes]
dts_pr = [(1./M)**2 for M in Ms]
dts_dr = dts_pr
questao_4_3_convergencia(Ms, 'tabela_pr_h2.txt', 'tabela_dr_h2.txt', dts_pr, dts_dr)