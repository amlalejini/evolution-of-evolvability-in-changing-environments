# Experiments

## Relevant Avida parameters

Avida parameters most relevant to these experiments.

- General
  - `RANDOM_SEED`
- Population structure
  - `BIRTH_METHOD` (4, mass action)
- Variable mutation rates
  - `META_COPY_MUT` (1.0)
  - `META_STD_DEV` (0.01)
  - `MUT_RATE_SOURCE` (1)
- Fixed-length genomes
  - `DIVIDE_INS_PROB` (0.0)
  - `DIVIDE_DEL_PROB` (0.0)
  - `OFFSPRING_SIZE_RANGE` (1.0)
  - `STERILIZE_UNSTABLE` (1)
- Mutation rates
  - `COPY_MUT_PROB` (0.01)
- Costly h-copy
  - `COSTLY_HEAD_COPY`
  - `MAX_HEAD_COPY_COST` (0.01)

## Relevant Avida output files

Non-standard output important for these experiments.

- `PrintMutationRateData`

## HPC

Running on MSU's HPC

```
module load CMake/3.18.4
module load GCC/11.2.0
CXX="g++"
```