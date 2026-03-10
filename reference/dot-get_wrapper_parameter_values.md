# Get matching parameters of a target function with current context

This is a helper function for wrapping functions such as plotting. It
gets the formal arguments of the target function func, and gets their
values from the context given in current, that will have the same name.
It returns a named list argument - value that can be passed to func
using do.call()

## Usage

``` r
.get_wrapper_parameter_values(func, current)
```

## Arguments

- func:

  Target function

- current:

  Current context as a named list

## Value

Named list argument - value
