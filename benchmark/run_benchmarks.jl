using PkgBenchmark

results = benchmarkpkg(".")
export_markdown("benchmark/latest_benchmark.md", results)

master_results = benchmarkpkg(".", "master")


comparison = judge(results, master_results)
export_markdown("benchmark/latest_comparison.md", comparison)
