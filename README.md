# GPDPS

`GPDPS` is a [Julia](http://julialang.org) package that provides the reference implementation of the generalized primal-dual proximal splitting approach described in the paper "[Primal-dual proximal splitting and generalized conjugation in non-smooth non-convex optimization](https://arxiv.org/abs/1901.02746)" by Christian Clason, Stanislav Mazurenko, and Tuomo Valkonen.

## Usage

To use the provided test codes (compatible with Julia version 1.0 and above, tested on macOS and Linux x86_64):

* start Julia from this directory with `julia --project=.` (or with the relative path to the `GPDPS.jl` directory from somewhere else)
* do `]instantiate` to download all dependencies (only required the first time, make sure you go back to standard Julia prompt with backspace afterwards)
* (highly recommended: `using Revise` so you can make changes in the code without having to restart Julia; if not install, start `julia` (without `--project`) and `]add Revise`)
* load the package with `using GPDPS`

To run the example for the elliptic Nash equilibrium problem:

* do `test_enep(N)`, where `N` controls the discretization (number of nodes per coordinate, default `N=128`)

To run the example for the Huber-Potts segmentation model:

* do `test_potts(alpha,gamma,keyword=value)`, where `keyword` is one or more of the following (comma seperated, order insensitive, may be removed before submission) with default value if omitted:

- `image`: test image; default is `"blobs"` (size 256x254), other images in `.tif` format can be specified if placed in the `img` folder
- `isotropic`: use isotropic (value `true`, default) or anisotropic (value `false)` Potts model
- `maxit`: maximum number of iterations (default 500000)

If `alpha` or `gamma` are not specified, the default values `1` and `1e-3` are used.

Possible issues:
* on macOS, it may also be required to `]add QuartzImageIO` if it is not included in the default environment.

## Reference

If you find this code useful, you can cite the paper as

    @article{GPDPS,
        author = {Clason, Christian and Mazurenko, Stanislav and Valkonen, Tuomo},
        title = {Primal-dual proximal splitting and generalized conjugation in non-smooth non-convex optimization},
        year = {2019},
        eprinttype = {arxiv},
        eprint = {1901.02746},
    }


