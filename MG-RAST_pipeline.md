Imagem da pipeline do MG-RAST:

![image](https://github.com/user-attachments/assets/e2b11016-9d1e-47e6-9c47-96ae206863f1)

Resumo das etapas realizadas:

**Pré-professamento:** Regiões de baixa qualidade são retiradas usando o programa SolexaQA (Cox, Peterson, and Biggs 2010). Reads maiores/menores que dois desvios padrões da média também são retirados.

**De-replication:** Retirada de reads duplicadas do conjunto de dados.

**DRISEE:** Montagem do DRISEE Score que quantifica as amostras duplicadas.

**Screening:** Retira contaminação de outras espécies utilizando o programa Bowtie.

**Gene Calling:** Utilizam a ferramenta FragGeneScan (Rho, Tang, and Ye 2010) baseada em Machine Learning.

**AA Clustering:** Junção de proteínas semelhantes em grupos com até 90% de semelhança.

**Protein Identification:** Cada cluster gera uma sequência representativa que é analisada. Utilizam o BLAST.

**Abundance Profile:** Perfil de abundância taxonômica ou funcional das sequências analisadas em um metagenoma. Ele mostra quais organismos (ou genes/funções) estão presentes na amostra e em que quantidade. Usa o NCBI Taxonomy para predizer a taxonomia. Alguns parâmetros da anotação das proteínas podem ser modificados para gerar o perfil de abundância: **E-value:** Mede a confiança do alinhamento. Quanto MENOR, MAIS PRECISO. **Percent Identity:** Porcentagem de similaridade entre a sua sequência e a referência. **Minimum Alignment Length:** Tamanho mínimo da região que precisa alinhar. **Minimum Abundance:** Número mínimo de vezes que uma anotação precisa aparecer pra ser incluída no profile.
