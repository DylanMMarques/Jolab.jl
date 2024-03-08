volume(grid::CartesianGrid{3}, ind) = prod(grid.spacing)
volume(grid::Domain, ind) = Meshes.volume(element(grid, ind))

measure(grid::CartesianGrid, ind) = prod(grid.spacing)
measure(grid::Domain, ind) = Meshes.measure(element(grid, ind))
