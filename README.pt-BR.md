# MEGAN_SC: emissões biogênicas de Santa Catarina, 2003-2025

Código, saídas do modelo e análises do artigo *Emissões biogênicas em Santa Catarina, Brasil, 2003-2025:
inventário com o MEGANv3.21 forçado por MapBiomas, MODIS e ERA5-Land* e do relatório *Inventário das emissões
biogênicas na qualidade do ar em Santa Catarina* (ano-base 2025, série 2003-2025). Os dois documentos citam este
repositório pelo endereço https://github.com/jessicajcss/MEGAN_SC.

> [README in English](README.md) · **Tutorial ilustrado em PDF: [`TUTORIAL_MEGAN_SC_pt.pdf`](TUTORIAL_MEGAN_SC_pt.pdf)** ([English](TUTORIAL_MEGAN_SC_en.pdf)) · [Capítulos técnicos (em inglês)](docs/tutorial/README.md)

O MEGAN_SC é uma implementação em Python do código de emissões MEGANv3.21 (Guenther, 2024b), com uma camada de dados
no Google Earth Engine. Para Santa Catarina, calcula em grade de 0,025° (144 × 224 células, das quais 13.970 dentro do
Estado), a cada 3 h de 2003 a 2025, a emissão das 19 classes do MEGAN: isopreno, monoterpenos, sesquiterpenos, COV
oxigenados e leves, NO do solo e CO. Os procedimentos de preparação das entradas e de execução (FIX-1 a FIX-8) estão em
`src/megan_sc/fixes.py`, cujo SHA-256 fica gravado em cada arquivo de saída, e o motor, tal como usado na produção,
reproduz etapa por etapa os executáveis Fortran oficiais na precisão de float32. Ao ler a execução, as análises
multiplicam o isopreno pelo fator de perda e produção no dossel, ρ = 0,96 (Guenther et al., 2006), que o código
MEGANv3.21 não aplica.

## Três níveis de uso

| Nível | Objetivo | Requisitos | Tempo |
|---|---|---|---|
| **A. Reproduzir** | recalcular todos os números, tabelas e figuras do artigo e do relatório a partir das saídas do modelo distribuídas com o repositório e compará-los com as cópias publicadas | Python 3.10 ou mais recente, cerca de 1,2 GB em disco, sem contas | cerca de 5 min de processamento |
| **B. Replicar** | executar de novo o modelo a partir dos dados brutos (ERA5-Land, MODIS, MapBiomas), as verificações, a comparação com o Fortran e a avaliação com o TROPOMI | conta Google com Earth Engine e Google Colab | horas de requisições ao Earth Engine por ano simulado; 41 EECU-horas no caderno do TROPOMI |
| **C. Estender** | montar o mesmo inventário para outro estado, província ou país | como em B, mais um arquivo de limites | como em B |

## Reprodução em quatro comandos (nível A)

```bash
git clone https://github.com/jessicajcss/MEGAN_SC megan-sc && cd megan-sc
pip install -r requirements-reproduce.txt          # ou: conda env create -f environment.yml && conda activate megan-sc
python tools/data_release.py fetch https://github.com/jessicajcss/MEGAN_SC/releases/download/v1.0.0   # os arquivos de dados (DATA_AVAILABILITY.md)
python reproduce.py check-data                     # SHA-256 de cada arquivo de dados
python reproduce.py all                            # 16 etapas; a última compara cada arquivo com results/v10
```

O endereço acima é o da versão (*release*) do GitHub, marcada v1.0.0, que vai receber os arquivos de dados quando o
repositório for publicado. A cópia arquivada no Zenodo, depois do depósito, ou uma pasta local com os arquivos zip
funcionam do mesmo modo. Tudo é gravado em `reproduced/`, sem alterar `data/` nem `results/`. A última etapa compara cada arquivo recalculado
com a cópia publicada. Numa execução de 15/11/2026, numa máquina Linux com versões de bibliotecas mais novas que as de
`requirements-lock.txt`, em uma única passada numa pasta vazia (5 min 22 s), os 58 arquivos CSV ficaram idênticos (26) ou
iguais dentro da tolerância (32), sem nenhum diferente, e 17 das 59 figuras PNG ficaram idênticas pixel a pixel
(`docs/tutorial/reproduction_log.txt`). Quatro dos arquivos CSV idênticos e essas 17 figuras foram gravados em
`results/v10` nessa máquina, e a comparação não os testa de forma independente. As demais figuras foram desenhadas no
Google Colab e diferem na renderização, por causa de fontes e versões de bibliotecas, sem mudança nos números, e os
arquivos iguais dentro da tolerância diferem no máximo em 9,5 × 10⁻⁷ em termos relativos.

## Cadernos no Google Colab

Os catorze cadernos Jupyter (*notebooks*) de `notebooks/` são executados no Google Colab com o Google Drive (níveis B e C), nesta ordem: 12 (comparação com os executáveis Fortran), 10 (execução de produção, Earth Engine), 11 (verificações C0 a C5), 13 (números, tabelas e figuras), 14 (avaliação com o formaldeído do TROPOMI, Earth Engine), 15 (testes de sensibilidade), 16 (comparação da queda do IAF do MODIS em 2024 e 2025 com o IAF do VIIRS, Earth Engine), 17 (2024 e 2025 com o IAF corrigido pela anomalia do VIIRS, Earth Engine), 18 (comparação com as emissões do MEGAN no BRAIN, 2019), 19 (verificação e instalação de 2024 e 2025 com o IAF corrigido e reprodução com os dois IAF), 21 (ERA5-Land contra estações do INMET e gradiente térmico vertical lido pelo caderno 20, Earth Engine), 20 (novas execuções de 2019 com uma configuração alterada por vez, Earth Engine), 22 (IAF do VIIRS do NOAA-20 e do NOAA-21, com login do Earthdata) e 23 (testes das explicações das diferenças em relação ao BRAIN). Cada um abre diretamente no Colab pelo endereço https://colab.research.google.com/github/jessicajcss/MEGAN_SC/blob/main/notebooks/<nome do caderno>.ipynb (ramo `main`); [notebooks/README.md](notebooks/README.md) descreve o que cada caderno precisa e produz, a célula de configurações e o tempo de execução. Informe o seu projeto do Earth Engine (`EE_PROJECT`) na célula "3 · Earth Engine" dos cadernos 10, 17 e 20, na célula "2 · Earth Engine" do caderno 21 e na célula de configurações dos cadernos 14 e 16. O caderno 19 instalou na execução os anos de 2024 e 2025 com o IAF corrigido e foi executado uma vez pelos autores; o `reproduce.py` repete as suas comparações.

## Onde está cada coisa

* **Tutorial ilustrado** (`TUTORIAL_MEGAN_SC_pt.pdf`, em português, e `TUTORIAL_MEGAN_SC_en.pdf`, em inglês; fontes
  HTML em `docs/tutorial/`): do que é o inventário à leitura e à verificação dos resultados, com glossário e lista de
  conferência.
* **Capítulos técnicos** (`docs/tutorial/`, em inglês): 13 capítulos, das equações do MEGAN ao mapa que liga cada número dos
  documentos ao arquivo e à coluna de onde vem (capítulo 8), cada figura ou tabela ao script que a produz (capítulo 9)
  e a forma de identificar o código, os dados e os resultados (capítulo 13).
* **Resultados publicados** (`results/v10/`): tabelas, números e figuras da execução de produção (caderno 10), com o IAF
  de 2024 e 2025 corrigido (cadernos 17 e 19). Os resultados com o IAF do MODIS sem correção estão em
  `results/v10_modis_lai_2024_2025/`.
* **Dados** (`data/`): limites do IBGE, tabelas do MEGAN e ecótipo ficam no git. Cinco arquivos zip chegam com
  `tools/data_release.py` (ver `DATA_AVAILABILITY.md`): as saídas do modelo (2003-2025), os anos de 2024 e 2025
  calculados com o IAF do MODIS, os campos mensais do TROPOMI, as duas execuções de teste do passo de tempo e o IAF
  mensal do VIIRS e do MODIS (caderno 16).
* **Validação** (`VALIDATION.md`): comparação com os executáveis Fortran, verificações automáticas da execução,
  avaliação com o TROPOMI e testes de sensibilidade.

## Citação e licenças

Os metadados de citação estão em `CITATION.cff`, que vão receber o DOI depois do depósito. O código está sob a licença
MIT; documentos, tabelas e figuras, sob CC-BY-4.0; os dados de terceiros mantêm as próprias licenças
(`LICENSE-docs-data.md`). O código MEGANv3.21 (Guenther, 2024b; Zenodo, DOI 10.5281/zenodo.10526206) fica fora deste
repositório. As tabelas de fatores de emissão e de especiação usadas pelo motor (`data/megan_inputs/`, sem alteração) e o
mapa de ecotipos recortado para Santa Catarina (`data/sc_ecotype.nc`) vêm dos arquivos de entrada do MEGAN 3.2
(Guenther, 2024a; Zenodo, DOI 10.5281/zenodo.10939297; licença CC BY 4.0). Cite também esses dois registros.
