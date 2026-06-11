# Method Notes for Public Reproducibility

## Velocity and LOS Convention

Velocity components are stored as east, north, and up in mm/yr. LOS geometry
rows contain the matching east, north, and up unit-vector components. LOS
observations are computed as:

```text
d_los = g_e * v_e + g_n * v_n + g_u * v_u + noise
```

## Baselines

WLS is the unconstrained three-component weighted least-squares inversion.

APFM uses a strict two-basis model: horizontal downslope motion plus an
independent vertical component.

SPFM uses a strict terrain-parallel two-basis model:

```text
v_u = dz_dx * v_e + dz_dy * v_n
```

## SATFM Weak Constraints

The public implementation adds weak normal equations for:

1. cross-slope horizontal motion,
2. terrain-parallel vertical consistency,
3. edge-aware spatial smoothing.

Soft-zone membership probabilities control the local directional and
terrain-parallel penalties:

```text
lambda_d(i) = gamma * r(i) * sum_k mu_k(i) * lambda_d_base(k)
lambda_t(i) = gamma * r(i) * sum_k mu_k(i) * lambda_t_base(k)
```

where `mu_k(i)` is the soft-zone membership and `r(i)` is the reliability
weight.

## Edge-Aware Smoothing

The manuscript Eq. (12) and Table 2 reproduction use:

```text
w_ij = exp(-dm2/tau_m^2) * exp(-ds/tau_s^2) * exp(-q_ij*da/tau_a^2)
```

Here `q_ij` is the edge gate. In the manuscript Table 2 reproduction,
`q_ij = c_dir(i) c_dir(j)`, where `c_dir` is the direction-confidence score.
This linear-`da` kernel with direction-confidence gating is labeled
`edge_kernel="legacy_linear"` in the run diagnostics.
