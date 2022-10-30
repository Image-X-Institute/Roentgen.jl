
#--- HDF5

for T in (MultiLeafCollimator, MultiLeafCollimatorSequence)
    eval(quote
        function save(file::HDF5.H5DataStore, mlc::$T)
            file["mlc/positions"] = mlc.positions
            file["mlc/edges"] = collect(mlc.edges)
            nothing
        end

        function load(::Type{$T}, file::HDF5.H5DataStore)
            positions = read(file["mlc/positions"])
            edges = read(file["mlc/edges"])
            $T(positions, edges)
        end
    end)
end
