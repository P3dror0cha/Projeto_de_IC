Imagem da pipeline do MG-RAST:

![image](https://github.com/user-attachments/assets/e2b11016-9d1e-47e6-9c47-96ae206863f1)

Resumo das etapas realizadas:
**Pré-professamento:** Regiões de baixa qualidade são retiradas usando o programa SolexaQA (Cox, Peterson, and Biggs 2010). Reads maiores/menores que dois desvios padrões da média também são retirados.
**De-replication:** Retirada de reads duplicadas do conjunto de dados.
**DRISEE:** Montagem do DRISEE Score que quantifica as amostras duplicadas.
**Screening:** Retira contaminação de outras espécies utilizando o programa Bowtie.
**Gene Calling:** 
