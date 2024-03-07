volume(grid::CartesianGrid{3}, ind) = prod(grid.spacing)
volume(grid::Domain, ind) = Meshes.volume(element(grid, ind))