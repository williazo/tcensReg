# tcensReg

Estimate parameters from a linear model when the data comes from a truncated normal distribution with censoring. Maximum likelihood values are returned derived from Newton-Raphson algorithm using analytic values of the gradient and hessian.

## Example

Some common examples where this type of problem may arise is when there is a natural truncation imposed by the structure of the data. For instance several applications have an implied zero truncation such as product lifetimes, age, or detection thresholds.

``` r
## basic example code
```
