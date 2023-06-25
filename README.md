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

|  | sample_A_1 | sample_A_2 |  sample_B_1 | sample_B_2 |
| --------------- | --------------- | --------------- | --------------- | --------------- |
| GENE_A | ## | ## | ## | ## |
| GENE_B | ## | ## | ## | ## |
| GENE_C | ## | ## | ## | ## |
| GENE_D | ## | ## | ## | ## |

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
