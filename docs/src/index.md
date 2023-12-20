![GraviPet logo](figures/GraviPet-light-background.jpg)

# GraviPet.jl

GraviPet is the **G**eneral **R**elativistic **A**strophysics
**V**isualization, **I**nitialization, and **P**ostprocessing
**E**fficient **T**oolkit.

```@autodocs
Modules = [GraviPet]
Private = false
```

## Domains
```@autodocs
Modules = [GraviPet.Domains]
Private = false
```

### Intervals
```@autodocs
Modules = [GraviPet.Intervals]
Private = false
```
### Boxes
```@autodocs
Modules = [GraviPet.Boxes]
Private = false
```

## Categories
```@autodocs
Modules = [GraviPet.Categories]
Private = false
```

### Julia Functions
```@autodocs
Modules = [GraviPet.JuliaFunctions]
Private = false
```
### 1D Grid Functions
```@autodocs
Modules = [GraviPet.GridFunctions1D]
Private = false
```
### Multi-dimensional Grid Functions
```@autodocs
Modules = [GraviPet.GridFunctions]
Private = false
```
### Blocked Functions (Domain Decomposition)
```@autodocs
Modules = [GraviPet.BlockFunctions]
Private = false
```
### Threaded Functions
```@autodocs
Modules = [GraviPet.ThreadedFunctions]
Private = false
```

## Miscellaneous
```@docs
Base.:(==)(::Domain, ::Domain)
Base.last(::Domain)
Base.isdisjoint(::Domain, ::Domain)
Base.eltype(::Domain)
Base.ndims(::Domain)
Base.issubset(::Domain, ::Domain)
Base.in(::Any, ::Domain)
Base.first(::Domain)
```

## Public Helpers
```@autodocs
Modules = [GraviPet.Funs]
```
