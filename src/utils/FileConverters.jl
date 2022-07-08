using Glob
#export read_eclipse_kernel_data, read_eclipse_measured_data

"""
    File Converters

This file contains methods to convert unfriendly filetypes into JSON.
"""

"""
    eclipse_to_SPBK(input_datadir)

Read in Eclipse exported files, extracting Single Pencil Beam Kernel data

Requires `SinglePencilBeamKernel{depth}` and `IntensityProf{depth}` files in 
`input_datadir` for multiple depths. Returns a dict with kernel and intensity
profile data.
"""
function eclipse_to_SPBK(input_datadir)
    
    kernel_filenames = glob(input_datadir*"SinglePencilBeamKernel*")
    kernel_data = read_eclipse_datafile.(kernel_filenames)
    
    intensityprof_filenames = glob(input_datadir*"IntensityProf*")
    intensityprof_data = read_eclipse_datafile.(intensityprof_filenames)

    kernel = Dict()
    kernel["Radius"], kernel["Depth"], kernel["Value"] = combine_filedata(kernel_data, "Radial", "Kernel Value")
    
    intensityprof = Dict()
    intensityprof["Radius"], intensityprof["Depth"], intensityprof["Dose Ratio"] = combine_filedata(intensityprof_data, "Radial", "Dose Ratio")

    data = Dict()
    data["version"] = "1.0"
    data["Algorithm"] = "Single Pencil Beam Kernel"
    data["Kernel"] = kernel
    data["Intensity Profile"] = intensityprof

    data
end

"""
    eclipse_to_SIK(input_datadir)

    Read in Eclipse exported files, extracting Scaled Isoplane Kernel data

Requires `SinglePencilBeamKernel{depth}` for multiple depths and `MeasuredProfile`
files in `input_datadir`. Returns a dict with kernel and depth dose data.
"""
function eclipse_to_SIK(input_datadir, depth=100., fieldsize=100.)
    
    kernel_data = read_eclipse_datafile(input_datadir*"SinglePencilBeamKernel$(convert(Int, depth))")

    depth_dose_data = read_eclipse_depth_dose_data(input_datadir*"MeasuredDepthDose")

    kernel = Dict()
    kernel["Radius"] = kernel_data["Radial"]
    kernel["Kernel Value"] = kernel_data["Kernel Value"]

    i = findfirst(x->x==fieldsize, depth_dose_data["Field Size [mm]"])

    depth, dose = truncate_zeros(depth_dose_data["Depth [mm]"], depth_dose_data["Dose [%]"][:, i])

    depth_dose = Dict()
    depth_dose["Field Size"] = depth_dose_data["Field Size [mm]"][i]
    depth_dose["Depth"] = depth
    depth_dose["Dose"] = dose

    data = Dict()
    data["version"] = "1.0"
    data["Algorithm"] = "ScaledIsoplaneKernel"
    data["Kernel"] = kernel
    data["Depth Dose"] = depth_dose

    data
end

function combine_filedata(data, x_name, value_name)
    value = hcat(getindex.(data, Ref(value_name))...)

    x = data[1][x_name]
    depth = getindex.(getindex.(data, Ref("depth")), 3)

    sort_order = sortperm(depth)

    depth = depth[sort_order]
    value = value[:, sort_order]

    x, depth, value
end

"""
    read_eclipse_datafile(filename::String, truncate=true)

Read an Eclipse data from file.

Reads Eclipse data from files such as `SinglePencilBeamKernel`, `IntensityProf`
and `EnvelopeProfile` files. The files are generally of the form:

    ```
    %attr1_name: attr1_value
    %attr2_name: attr2_value
        ...
    <x1 y1>
    <x2 y2>
        ...
    <xn yn>
    ```

Args:
- `filename`: File path
- `truncate`: Whether to remove zeros, default `true` (see `truncate_zeros`)
"""
function read_eclipse_datafile(filename::String, truncate=true)

    data = Dict()
    nskip = 0

    open(filename, "r") do file
        for (i, line) in enumerate(eachline(file))
            if(startswith(line, "%") && ':' in line)

                name, string_value = split(line, ": ")
                data[name[2:end]] = parse_value(string_value)
            end

            if(startswith(line, "<"))
                nskip = i-1
                break
            end
        end
    end

    dataset = readdlm(filename, skipstart=nskip, comments=true, comment_char='$')
    x = [parse(Float64, d[2:end]) for d in dataset[:, 1]]
    y = [parse(Float64, d[1:end-1]) for d in dataset[:, 2]]

    if(truncate)
        x, y = truncate_zeros(x, y)
    end
    data[data["axis legend"]] = x
    data[data["data legend"]] = y

    data
end

"""
    read_eclipse_depth_dose_data(filename)

Read an Eclipse exported depth dose file.

This file is of the form:

    ```
    %attr1_name: attr1_value
    %attr2_name: attr2_value
        ...
        ,    col1,    col2, ...,    colm
    row1,     x11,     x12, ...,    col1m
    row2,     x21,     x22, ...,    col2m
        ...
    rown, ...
    ```

Returns a dictionary with attributes, the column header (as given by `column legend`),
the row index (`row legend`), and a matrix containing the data.
"""
function read_eclipse_depth_dose_data(filename)

    data, nskip = read_attributes(filename)

    dataset, headers = readdlm(filename, ',', Float64, skipstart=nskip, header=true)

    data[data["row legend"]] = dataset[:, 1]
    data[data["data legend"]] = dataset[:, 2:end]
    data[data["column legend"]]  = parse.(Float64, headers[2:end]) 

    data
end

"""
    read_attr(line)

Read attribute name and value
"""
read_attr(line) = split(line, ": ")

"""
    read_attributes(filename)

Read attributes from Eclipsed exported measured data file.
See read_eclipse_depth_dose and read_eclipse_profile.
"""
function read_attributes(filename)

    data = Dict()
    nskip = 0

    # Read the attributes
    open(filename, "r") do file
        for (i, line) in enumerate(eachline(file))
            !in(':', line) && break
            
            name, string_value = read_attr(line)
            data[name] = parse_value(string_value)
            nskip += 1

        end
    end
    data, nskip
end

"""
    parse_value(string_value)

Tries to parse a string, returning string if parsing failed
"""
function parse_value(string_value)

    string_array = split(string_value)

    value = nothing
    for T in [Int, Float64]
        if(tryparse(T, string_array[1]) != nothing)
            value = parse.(T, string_array)
            length(value) == 1 && return value[1]
            return value
        end
    end
    
    value == nothing && return string_value
end

"""
    truncate_zeros(x, y)

Removes every element after first zero in `y` and corresponding `x`
"""
function truncate_zeros(x, y)
    for i in eachindex(y)
        if(y[i]<=0.)
            return x[1:i-1], y[1:i-1]
        end
    end
    x, y
end

"""
    read_eclipse_profile_data(filename)

Read an Eclipse exported profile file.

Similar to read_eclipse_depth_dose_data, except profile files can have multiple
datasets at different depths. This file is of the form:

    ```
    %attr1_name: attr1_value
    %attr2_name: attr2_value
        ⁞
    Profile at d1
        ,    col1,    col2, ...,    colm
    row1,     x11,     x12, ...,    col1m
        ⁞
    rown, ...
    Profile at d2
        ,    col1,    col2, ...,    colm
    row1,     x11,     x12, ...,    col1m
        ⁞
    rown, ...
        ⁞
    ```

Returns a dictionary with attributes and profile data. Profile data is stored
in a vector under data["Profiles"].
"""
function read_eclipse_profile_data(filename)

    data, _ = read_attributes(filename)

    lines = Vector{String}(undef, 0)
    open(filename, "r") do file
        lines = readlines(file)
    end

    dataset_linenumbers = get_dataset_line_numbers(lines)
    depths = [parse(Float64, read_attr(lines[linenumber])[end]) for linenumber in dataset_linenumbers[1:end-1]]
    lines_in_each_dataset = split_lines_into_datasets(lines, dataset_linenumbers)

    datasets = read_dataset.(lines_in_each_dataset)

    profiles = []
    for (depth, dataset) in zip(depths, datasets)
        profile = Dict()
        rows, cols, dset = dataset
        profile[data["row legend"]] = rows
        profile[data["column legend"]] = cols
        profile[data["data legend"]] = dset
        profile["Depth [mm]"] = depth
        
        push!(profiles, profile)
    end
    data["Profiles"] = profiles

    data
end

"""
get_dataset_line_numbers(lines)

Get the line numbers of each data section.

Looks for the string "Profiles" in each line, then adds the line number to the
list of dataset line numbers, which is returned. See read_eclipse_profile_data
"""
function get_dataset_line_numbers(lines)
    sections = Int[]
    for (linenumber, line) in enumerate(lines)
        if(startswith(line, "Profiles"))
            push!(sections, linenumber)
        end
    end
    push!(sections, length(lines))
    sections
end

"""
    split_lines_into_datasets(lines, sections)

Split the lines into each dataset. See read_eclipse_profile_data
"""
function split_lines_into_datasets(lines, sections)
    [lines[start+1:stop-1] for (start, stop) in zip(sections[1:end-1], sections[2:end])]
end


"""
    read_data_line(line; delim=',')

Read a line of data of the form:
    ```x1, x2, ...```
where ',' is the delimiter
"""
read_data_line(line; delim=',') = tryparse.(Float64, split(line, delim))

"""
    read_dataset(lines)

Read each line in `lines` and extract the header row, row index and data.
"""
function read_dataset(lines)
    cols = read_data_line.(lines[1])[2:end]
    data = hcat(read_data_line.(lines[2:end])...)
    rows = data[1, :]
    data = data[2:end, :]

    rows, cols, permutedims(data, (2, 1))
end

function load_depth_data(filename_glob, grid)
    files = glob(filename_glob)

    # Read each data file
    data = [read_data_file(filename) for filename in files]

    # Extract Data
    depth = getindex.(data, "depth")

    for d in data
        d["radius"], d["value"] = truncate_zeros(d["radius"], d["value"])
        d["value"] = log.(d["value"])
    end

    Δx, Δy = spacing(grid)

    i = eachindex(grid)
    r = @. √( (getindex(grid, 1) + 0.5*Δx)^2 + (getindex(grid, 2) + 0.5*Δy)^2)

    # Interpolate each profile over a common radius
    values = @. interp(getindex(data, "radius"), getindex(data, "value"), r[:]')
    values = reshape(values, length(depth), size(grid)...)
    # values = collect(reshape(values', size(grid)..., length(depth)))

    # Sort Depths
    i = sortperm(depth)
    depth = depth[i]
    values = values[i, :, :]

    prepend!(depth, 0)
    # values = cat(zeros(1, size(values, 2), size(values, 3)), values, dims=1);
    values = cat(0.5 .*reshape(values[end, :, :], 1, size(values, 2), size(values, 3)),
                 values, dims=1);

    depth, values
end
