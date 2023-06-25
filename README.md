# DifferentialExpression

## About

Set of R scripts to run DESeq2 using a variety of different methods

## Installation

Dependencies can be installed using conda:

```bash
conda env create -f deseq.yml
```

## Input files

counts.txt

|  | sample_1 | sample_2 |  sample_3 | sample_4 |
| --------------- | --------------- | --------------- | --------------- | --------------- |
| GENE_A | ## | ## | ## | ## |
| GENE_B | ## | ## | ## | ## |
| GENE_C | ## | ## | ## | ## |
| GENE_D | ## | ## | ## | ## |

```tsv
    sample_1    sample_2    sample_3    ... sample_n
GENE_A 15   6   12  ... 15
GENE_B 32   21   17  ... 15
GENE_C 37   43   20  ... 5
GENE_D 24   25   21  ... 35
```

samples.txt (for 'Wald' and 'LRT')

|  | condition | sample |
| --------------- | --------------- | --------------- |
| sample_A_1 | condition1 | sample_A_1 |
| sample_A_2 | condition1 | sample_A_2 |
| sample_B_1 | condition2 | sample_B_1 |
| sample_B_2 | condition2 | sample_B_2 |


samples.txt (for 'RIP')

|  | condition | assay |sample |
| --------------- | --------------- | --------------- | --------------- |
| sample_A_1 | condition1 | Input | sample_A_1 |
| sample_A_2 | condition1 | Input | sample_A_2 |
| sample_A_1 | condition1 | IP | sample_A_1 |
| sample_A_2 | condition1 | IP | sample_A_2 |
| sample_B_1 | condition2 | Input | sample_B_1 |
| sample_B_2 | condition2 | Input | sample_B_2 |
| sample_B_1 | condition2 | IP | sample_B_1 |
| sample_B_2 | condition2 | IP | sample_B_2 |
