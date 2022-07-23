# ElectrochemicalKinetics.jl benchmarks
The files in this folder define a benchmark suite with the tools provided by
[PkgJogger](https://github.com/awadell1/PkgJogger.jl) and
[BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl). PkgJogger is similar to [PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl), but it plays nicer with [Revise.jl](https://github.com/timholy/Revise.jl), so you can change source code and benchmark again without having to reload everything!

To run the benchmarks "from scratch", execute:

```julia
using PkgJogger
@jog ElectrochemicalKinetics
results = JogElectrochemicalKinetics.benchmark()
```
`results` will be of type `BenchmarkGroup`.

(Note that neither PkgJogger nor BenchmarkTools are a direct dependency of ElectrochemicalKinetics, so you will have to have them installed in another environment that is also loaded (I keep them in my base environment).)

To save results to file, either pass `save=true` to the `benchmark` function above, or do:

```
julia> JogElectrochemicalKinetics.save_benchmarks(filepath)
```

and read them back in later with
```
julia> results = PkgJogger.load_benchmarks(filename)
```

## Comparing benchmarks
Given two `BenchmarkGroup` objects, you can compare their medians using:
```julia
julia> PkgJogger.judge(results1, results2)
```
This also works with filenames of saved results, and you can pass different functions to `metric` (default is `median`) by which to compare.

## Pretuning parameters
Tuning the benchmark parameters (number of samples etc.) can often take longer than actually running the benchmarks, so it can be advantageous to pretune and save/load parameters rather than benchmarking with `retune=true` every time. To do so, do:

```julia
# pretune and save
suite = JogElectrochemicalKinetics.suite()
PkgJogger.tune!(suite)
BenchmarkTools.save("params.json", params(suite))

# reload later
loadparams!(JogElectrochemicalKinetics.suite(), BenchmarkTools.load("params.json")[1], :evals, :samples);
```

## Other options
### Running more exhaustive benchmarks
Want better stats? Change the value of `time_mult` to be something larger. Similarly, make it smaller to run more quickly, but note that you may get very few samples on some of the benchmarks, leading to noisy results.

### Running a subset of benchmarks
To just run certain tags,  do e.g.
```julia
using BenchmarkTools
suite = JogElectrochemicalKinetics.suite()[@tagged "tag1" && "tag2"] # for example
run(suite) # or pretune, etc. first
```
