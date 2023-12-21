# Categories
```@autodocs
Modules = [GraviPet.Categories]
Private = false
```

## Julia Functions
```@autodocs
Modules = [GraviPet.JuliaFunctions]
Private = false
```
## 1D Grid Functions
```@autodocs
Modules = [GraviPet.GridFunctions1D]
Private = false
```
## Multi-dimensional Grid Functions
```@autodocs
Modules = [GraviPet.GridFunctions]
Private = false
```
## Blocked Functions (Domain Decomposition)
```@autodocs
Modules = [GraviPet.BlockFunctions]
Private = false
```
## Threaded Functions
```@autodocs
Modules = [GraviPet.ThreadedFunctions]
Private = false
```
## Distributed Functions
```@autodocs
Modules = [GraviPet.DistributedFunctions]
Private = false
```
## Kernel Functions
```@autodocs
Modules = [GraviPet.KernelFunctions]
Private = false
```

# Miscellaneous
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
