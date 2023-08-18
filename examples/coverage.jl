
using Pkg, Printf
using Coverage
using PrettyTables

clean_folder(".")

Pkg.activate("."); Pkg.test(coverage=true)

function process_coverage()
    # process '*.cov' files
    coverage = process_folder(); # defaults to src/; alternatively, supply the folder name as argument
    # process '*.info' files, if you collected them
    coverage = merge_coverage_counts(coverage, filter!(
        let prefixes = (joinpath(pwd(), "src", ""),
                        joinpath(pwd(), "deps", ""))
            c -> any(p -> startswith(c.filename, p), prefixes)
        end,
        LCOV.readfolder("test")));
end

function coverate_stats(coverage)
    covered_lines, total_lines = get_summary(coverage)
    covered_lines, total_lines, 100*covered_lines/total_lines
end
getfilename(path) = splitdir(path)[end]

function summarise(coverage)
    println("Total: ", (@sprintf "%i/%i, %0.1f" coverate_stats(coverage)...), "%")

    filenames = getfilename.(getfield.(coverage, :filename))

    data = hcat(collect.(coverate_stats.(coverage))...)'
    data = hcat(filenames, data)

    s = reverse(sortperm(data[:, 3]))

    pretty_table(data[s, :], header=["Filename", "No. Lines Covered", "Total No. Lines", "Coverage (%)"],
                formatters = ft_printf("%.1f", 4), )
end

coverage = process_coverage();
summarise(coverage)
