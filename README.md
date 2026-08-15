# s*: minimum MinHash sketch size
Run: `chmod +x 00_orchestrator.sh && ./00_orchestrator.sh` (add `--force` to ignore checkpoints).
Outputs: `results/`, `figures/`, `rds/`, `logs/`.

## s*: with and without a saturation floor (R, Python, C, Rust)

R:
```r
s_star <- function(J, rho = 0.10, a = 0.05) {
    return qnorm(1 - a/2)^2 * (1 - J)/(rho^2 * J)
}

s_star_sat <- function(J, rho = 0.10, a = 0.05, delta = 0.05) {
    return max(ceiling(s_star(J, rho, a)), ceiling(log(delta)/log(J)))
}
```

Python:
```python
from math import log, ceil 
from scipy.stats import norm

def s_star(J, rho=0.10, a=0.05):
    return norm.ppf(1-a/2)**2*(1-J)/(rho**2*J)
    
def s_star_sat(J, rho=0.10, a=0.05, delta=0.05):
    return max(ceil(s_star(J,rho,a)), ceil(log(delta)/log(J)))
```

C:
```c
double s_star(double J, double rho, double z){ 
    return z * z * (1-J)/(rho * rho * J); 
}

double s_star_sat(double J, double rho, double z, double delta){ 
    double a = ceil(s_star(J,rho,z)), b = ceil(log(delta)/log(J)); 
    return a > b ? a : b; 
}
```

Rust:
```rust
fn s_star(j: f64, rho: f64, z: f64) -> f64 { 
    z * z * (1.0 - j)/(rho * rho * j) 
}

fn s_star_sat(j: f64, rho: f64, z: f64, delta: f64) -> f64 { 
    s_star(j, rho, z).ceil().max((delta.ln()/j.ln()).ceil()) 
}
```

Examples (R):
```r
s_star(J = 0.13) # 2570.8, no saturation floor
s_star_sat(J = 0.98, delta = 0.05) # 149, saturation floor applied
```
