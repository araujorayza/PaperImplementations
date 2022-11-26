# [Ultimate boundedness sufficient conditions for nonlinear systems using TS fuzzy modelling](https://doi.org/10.1016/j.fss.2018.03.010)

Recebido de Michele Valentino:

```
Estou te enviando os códigos no anexo que rodam o exemplo do meu artigo e vou tentar te explicar um pouquinho sobre eles aqui. Você irá usar somente os três códigos citados abaixo. Para rodá-los você precisará instalar no Matlab os pacotes Yalmip e Sedumi, os quais também estão no anexo.

No arquivo exemplo_2_ODE.m:
- Nas linhas  6 a 13, construímos os dois subsistemas do exemplo
- Na linha 19 temos a nossa função V
- Na linha 20 e 21 temos a derivada de V com relação ao subsistema 1 e subsistema 2
- Nas linhas 22 a 62 temos a lei de chaveamento, na qual acionamos fora de \Omega_l o subsistema que tem a menor derivada e dentro de omega_l ficamos chaveando de uma para outra

No arquivo lmi_chav_fuzzychaveado.m:
- Descrevemos ai todas as nossas LMIs

No arquivo Exemplo_2.m:
- Até a linha 24 descrevemos nossos parâmetros e nossas matrizes dos modelos locais
- Na linha 33 verificamos a factibilidade das LMIs
- Nas linhas 45 a 66 descrevemos a derivada da função V com relação a combinação convexa dos subsistemas
- Nas linhas 67 a 78 verificamos se a derivada comentada no item anterior é positiva, em caso positivo pedimos para plotar o ponto onde isso ocorre (assim formamos os pontinhos cinzas da Fig 1 (pontos onde a derivada é maior ou igual a zero))
- Na linha 79 pegamos o máximo valor que a V assume nesse conjunto onde a derivada assume valores positivos
- Na linha 80, escolhemos LM=0.13 que é maior que o valor comentado acima
- Nas linhas 82 e 83 construímos as curvas de nível que aparecem na Fig 1

Acredito que a parte que tenha mais interesse, é esta da linha 79.
```