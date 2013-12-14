This is test code (scala) for parallel Gibbs sampling.

Input: factor graph, in Denny's format
  - file directory: eval
  - files: varialbes.tsv, weights.tsv, factors.tsv 
Output: histogram file
Usage: FactorTest arg0 arg1
arg0: int, number of samples
arg1: int, number of threads


```
sbt ~compile
```

Example call:

```
sbt "run-main GibbsSampler -v src/main/resources/eval/denny/variables.tsv -f src/main/resources/eval/denny/factors.tsv -w src/main/resources/eval/denny/weights.tsv -n 10 -t 2 -o target/out.out"
```
