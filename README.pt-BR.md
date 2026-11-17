# MEGAN_SC: emissões biogênicas com o MEGANv3.21 e o Google Earth Engine, para qualquer região e período

O MEGAN_SC é uma implementação em Python do código de emissões MEGANv3.21 (Guenther, 2024b), com uma camada de dados
no Google Earth Engine. Para uma região dada por um arquivo de limites e um período dado em meses, ele calcula numa
grade regular de latitude e longitude, a cada 3 h, a emissão das 19 classes do MEGAN: isopreno, monoterpenos,
sesquiterpenos, COV oxigenados e leves, NO do solo e CO. A meteorologia vem do ERA5-Land, o índice de área foliar (IAF)
do MODIS, a cobertura da terra do MapBiomas (Brasil) ou do MODIS (qualquer país) e os fatores de emissão dos arquivos de
entrada do MEGAN 3.2. Os procedimentos de preparação das entradas e de execução (FIX-1 a FIX-8) estão em
`src/megan_sc/fixes.py`, cujo SHA-256 fica gravado em cada arquivo de saída, e o motor, tal como usado na produção,
reproduz etapa por etapa os executáveis Fortran oficiais na precisão de float32.

O repositório traz um exemplo completo: Santa Catarina, 2003-2025, em grade de 0,025° (144 × 224 células, das quais
13.970 dentro do Estado). As saídas do modelo, todas as tabelas, números e figuras calculados a partir delas e um
comando que recalcula tudo acompanham o exemplo. Qualquer pessoa pode conferi-lo no próprio computador e aplicar o mesmo
procedimento a outra região.

> [README in English](README.md) · **Tutorial ilustrado em PDF: [`TUTORIAL_MEGAN_SC_pt.pdf`](TUTORIAL_MEGAN_SC_pt.pdf)** ([English](TUTORIAL_MEGAN_SC_en.pdf)) · [Capítulos técnicos (em inglês)](docs/tutorial/README.md)

## O que o repositório permite

| Objetivo | Requisitos | Tempo |
|---|---|---|
| **Reproduzir o exemplo**: recalcular todos os números, tabelas e figuras do exemplo de Santa Catarina a partir das saídas do modelo distribuídas e compará-los com as cópias publicadas | Python 3.10 ou mais recente, cerca de 1,2 GB em disco, sem contas | alguns minutos de processamento |
| **Executar o modelo para a sua região e período**: arquivo de limites, assistente, execução, verificações e resultados | conta Google com um projeto do Earth Engine, Google Colab e Google Drive | minutos para um teste; horas de requisições ao Earth Engine por ano simulado |
| **Executar de novo o exemplo a partir dos dados brutos**, com a comparação com o Fortran, a avaliação com o TROPOMI e os testes de sensibilidade | como no item anterior | como no item anterior; 41 EECU-horas no caderno do TROPOMI |

## Reprodução do exemplo em cinco comandos

```bash
git clone https://github.com/jessicajcss/MEGAN_SC.git MEGAN_SC && cd MEGAN_SC
pip install -r requirements-reproduce.txt          # ou: conda env create -f environment.yml && conda activate megan-sc
python tools/data_release.py fetch https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0   # os arquivos de dados (DATA_AVAILABILITY.md)
python reproduce.py check-data                     # SHA-256 de cada arquivo de dados
python reproduce.py all                            # todas as etapas; a última compara cada arquivo com results/v10
```

O endereço acima é o da versão (*release*) v1.0.0 do GitHub, que guarda os arquivos de dados. Uma pasta local com os
arquivos zip funciona do mesmo modo. Tudo é gravado em `reproduced/`, sem alterar `data/` nem `results/`, e a última
etapa compara cada arquivo recalculado com a cópia publicada. Em 16 de novembro de 2026, numa máquina Linux com as versões de `requirements-lock.txt`, as 16 etapas levaram 3 min 37 s
e todos os arquivos CSV coincidiram com a cópia publicada (56 de 61 byte a byte, 5 dentro de uma tolerância relativa de
10⁻⁶). Dos 46 arquivos PNG, 34 ficaram idênticos pixel a pixel. As figuras dependem das fontes e das versões das
bibliotecas, e em outra máquina algumas podem mudar em poucos pixels ou na posição dos painéis, sem mudança nos números.
O registro completo, arquivo por arquivo, está em `docs/tutorial/reproduction_log.txt`.

## A sua região e o seu período

Tudo no Google Colab, com o repositório copiado para o Google Drive como `MyDrive/MEGAN_SC` (o capítulo 3 do tutorial
mostra como copiar e como registrar um projeto do Earth Engine):

1. `run_your_region/region_setup.ipynb`: informe o arquivo de limites, o primeiro e o último mês, o projeto do Earth
   Engine e a resolução. O assistente monta a grade, lê a cobertura da terra, roda um teste curto com o motor completo
   e grava `config/<região>_run.yaml` ([START_HERE.md](run_your_region/START_HERE.md));
2. `notebooks/10_production_run.ipynb` com `RUN_CONFIG = 'config/<região>_run.yaml'`: a execução do modelo, mês a mês,
   retomada automaticamente depois de uma desconexão;
3. `notebooks/11_verify_run.ipynb`: seis verificações automáticas de cada ano;
4. `notebooks/25_results_any_region.ipynb`: totais por classe e por grupo, ciclo sazonal, mapas de intensidade,
   variáveis de controle, correlações e tendências.

O capítulo 8 do tutorial lista as escolhas a revisar em outra região, e o capítulo 9 explica como ler os resultados.

## Cadernos

Os cadernos Jupyter (*notebooks*) de `notebooks/` rodam no Google Colab com o Google Drive, ou no Jupyter num computador
com credenciais do Earth Engine. Todos começam pelo mesmo bloco, que encontra o repositório e, no Colab, conecta o
Drive. Os cadernos 10, 11 e 25 servem a qualquer região. Os demais acompanham o exemplo: 12 (comparação com os
executáveis Fortran), 13 (tabelas, números e figuras de Santa Catarina, com mesorregiões, municípios e o inventário
estadual anterior), 14 (avaliação com o formaldeído do TROPOMI, Earth Engine), 15 (testes de sensibilidade), 16
(comparação da queda do IAF do MODIS em 2024 e 2025 com o IAF do VIIRS, Earth Engine), 17 (2024 e 2025 com o IAF
corrigido, Earth Engine), 18 (comparação com as emissões do MEGAN no BRAIN, 2019), 19 (instalação de 2024 e 2025 com o
IAF corrigido e reprodução com os dois IAF, executado uma vez na preparação do exemplo), 20 (novas execuções de 2019 com
uma configuração alterada por vez, Earth Engine), 21 (ERA5-Land contra estações do INMET, Earth Engine), 22 (IAF do
VIIRS do NOAA-20 e do NOAA-21, com login do Earthdata) e 23 (testes das explicações das diferenças em relação ao BRAIN).
Cada um abre no Colab pelo endereço https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/<nome do caderno>.ipynb,
desde que o repositório esteja no Drive. [notebooks/README.md](notebooks/README.md) descreve o que cada caderno precisa e
produz, a célula de configurações e o tempo de execução.

## Onde está cada coisa

* **Tutorial ilustrado** (`TUTORIAL_MEGAN_SC_pt.pdf`, em português, e `TUTORIAL_MEGAN_SC_en.pdf`, em inglês; fontes
  HTML em `docs/tutorial/`): o MEGANv3.21 em palavras simples, os dados, o repositório, o exemplo, o Google Drive, o
  Colab e o Earth Engine, a sua região, a leitura dos resultados, glossário e lista de conferência.
* **Capítulos técnicos** (`docs/tutorial/`, em inglês): 12 capítulos, das equações do MEGAN à execução para outra região,
  à leitura dos resultados e à forma de identificar o código, os dados e os resultados.
* **Configurações** (`config/`): `santa_catarina_2003_2025.yaml`, a do exemplo; o assistente grava ali as novas.
* **Resultados do exemplo** (`results/v10/`): tabelas, números e figuras da execução de Santa Catarina (caderno 10), com
  o IAF de 2024 e 2025 corrigido (cadernos 17 e 19). Os resultados com o IAF do MODIS sem correção estão em
  `results/v10_modis_lai_2024_2025/`.
* **Dados** (`data/`): limites do IBGE, tabelas do MEGAN e ecótipo ficam no git. Cinco arquivos zip chegam com
  `tools/data_release.py` (ver `DATA_AVAILABILITY.md`): as saídas do modelo (2003-2025), os anos de 2024 e 2025
  calculados com o IAF do MODIS, os campos mensais do TROPOMI, as duas execuções de teste do passo de tempo e o IAF
  mensal do VIIRS e do MODIS (caderno 16).
* **Validação** (`VALIDATION.md`): comparação com os executáveis Fortran, verificações automáticas da execução,
  avaliação com o TROPOMI e testes de sensibilidade.

## Citação e licenças

Os metadados de citação estão em `CITATION.cff`, que vai receber o DOI depois do depósito. O código está sob a licença
MIT; documentos, tabelas e figuras, sob CC-BY-4.0; os dados de terceiros mantêm as próprias licenças
(`LICENSE-docs-data.md`). O código MEGANv3.21 (Guenther, 2024b; Zenodo, DOI 10.5281/zenodo.10526206) fica fora deste
repositório. As tabelas de fatores de emissão e de especiação usadas pelo motor (`data/megan_inputs/`, sem alteração) e o
mapa de ecótipos recortado para Santa Catarina (`data/sc_ecotype.nc`) vêm dos arquivos de entrada do MEGAN 3.2
(Guenther, 2024a; Zenodo, DOI 10.5281/zenodo.10939297; licença CC BY 4.0). Cite também esses dois registros.
