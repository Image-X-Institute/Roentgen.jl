# Benchmarks

This subdirectory contains scripts to run various benchmarks on the code.
It uses [PkgBenchmarks.jl](https://juliaci.github.io/PkgBenchmark.jl/stable/), which is based on [BenchmarkTools.jl](https://juliaci.github.io/BenchmarkTools.jl/stable/), which has useful methods to re-run benchmarks, compare with previous versions and export the results to markdown.
The benchmarks can be found in `benchmarks.jl`. 

See `latest_comparison.md` for a comparison with `master`.
The benchmark results below and in `latest_comparison.md` were generated using `run_benchmarks.jl`.

**Currently implemented benchmarks:**

- Dose-Fluence Matrix Construction
- Fluence calculation from MLC aperture
- Dose Calculation

## Usage

Start the Julia REPL in the root directory, and activate the local environment `] activate .`

Run the benchmarks with
```julia
julia> results = benchmarkpkg(".");
```

To view the results of a given benchmark, you can index into the results
```julia
julia> results.benchmarkgroup["dose-calculation"]
BenchmarkTools.Trial: 652 samples with 1 evaluation.
 Range (min … max):  6.314 ms … 19.932 ms  ┊ GC (min … max): 0.00% … 54.46%
 Time  (median):     7.140 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   7.495 ms ±  1.262 ms  ┊ GC (mean ± σ):  0.78% ±  4.38%

   ▁██▁▁▂    ▁   
  ▄██████▇████▇▅▅▅▅▄▅▇▇▆▄▃▃▂▃▃▄▃▃▃▃▃▂▃▃▃▃▂▂▃▂▂▂▂▂▂▂▃▃▂▁▂▂▂▁▂ ▃
  6.31 ms        Histogram: frequency by time        11.2 ms <

 Memory estimate: 566.45 KiB, allocs estimate: 2.
```

Finally, [PkgBenchmarks.jl](https://juliaci.github.io/PkgBenchmark.jl/stable/) has a useful `export_markdown` function, which will format the results and save it to a markdown file.
The output of this is shown in the Results section

# Results

## Job Properties
* Time of benchmark: 1 Nov 2021 - 17:7
* Package commit: 6c3b99
* Julia commit: 1b93d5
* Julia command flags: None
* Environment variables: None

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                 | time            | GC time | memory          | allocations |
|------------------------------------|----------------:|--------:|----------------:|------------:|
| `["dose-calculation"]`             |   7.157 ms (5%) |         | 566.45 KiB (1%) |           2 |
| `["dose-fluence-matrix-creation"]` |  51.382 ms (5%) |         |  90.16 MiB (1%) |          10 |
| `["fluence-from-mlc"]`             | 606.700 μs (5%) |         | 266.89 KiB (1%) |          19 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `[]`

## Julia versioninfo
```
Julia Version 1.6.2
Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
      Microsoft Windows [Version 10.0.19043.1288]
  CPU: Intel(R) Core(TM) i7-8650U CPU @ 1.90GHz: 
              speed         user         nice          sys         idle          irq
       #1  2112 MHz   20385046            0     24968703    557006500      2417171  ticks
       #2  2112 MHz   12122234            0     12175234    578062468       287640  ticks
       #3  2112 MHz   27103484            0     21972031    553284406       253015  ticks
       #4  2112 MHz   18430828            0     11908937    572020156       166687  ticks
       #5  2112 MHz   25640265            0     19157046    557562609       219453  ticks
       #6  2112 MHz   17212250            0     12047171    573100500       143671  ticks
       #7  2112 MHz   26792906            0     18500453    557066562       189125  ticks
       #8  2112 MHz   25851640            0     14004109    562504171       134437  ticks
       
  Memory: 31.881214141845703 GB (17000.8984375 MB free)
  Uptime: 979837.0 sec
  Load Avg:  0.0  0.0  0.0
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, skylake)
```