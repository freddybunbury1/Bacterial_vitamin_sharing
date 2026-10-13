# Run from the release root. No production outputs are changed.
source('R/utils/b12_model_structure_comparison.R')
stopifnot(requireNamespace('deSolve', quietly=TRUE))
times <- c(240, 0, 0.25, 5, 20, 60, 120, 5)
grid <- expand.grid(k=c(1e-8, 0.05, 5, 1000), capacity=c(10, 1000, 1e6),
                    od=c(0, 0.05, 0.5), b0=c(0, 500), release=c(0, 10))
max_error <- 0
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  pars <- c(k_per_od=g$k, Qmax_per_od=g$capacity)
  a <- quota_predict_state(times,g$b0,g$od,pars,g$release)
  b <- quota_predict_state(times,g$b0,g$od,pars,g$release,
                           rtol=1e-10,atol=1e-10,force_adaptive=TRUE)
  error <- max(abs(a$B-b$B)); max_error <- max(max_error,error)
  stopifnot(error < 1e-4+1e-7*max(1,abs(b$B)),
            max(abs(a$B+g$od*a$Q-(g$b0+g$od*g$release*times))) <
              1e-4+1e-7*max(1,g$b0+g$od*g$release*times),
            a$B[2]==g$b0, a$Q[2]==0, a$B[4]==a$B[8])
}
# Capacity equals the initial concentration, and almost equals it.
for (delta in c(0,-1e-10,1e-10)) {
  p <- c(k_per_od=5,Qmax_per_od=(500+delta)/0.5)
  a <- quota_predict_state(times,500,0.5,p)
  b <- quota_predict_state(times,500,0.5,p,rtol=1e-10,atol=1e-10,force_adaptive=TRUE)
  stopifnot(max(abs(a$B-b$B))<1e-4)
}
# Nonzero initial uptake, already full quota, and zero-time-only inputs.
p <- c(k_per_od=5,Qmax_per_od=1000)
a <- quota_predict_state(times,500,0.5,p,q0=200)
b <- quota_predict_state(times,500,0.5,p,q0=200,force_adaptive=TRUE)
stopifnot(max(abs(a$B-b$B))<1e-4,
          all(quota_predict_state(times,500,0.5,p,q0=1000)$B==500),
          all(quota_predict_state(c(0,0),500,0.5,p,release_per_od=1)$B==500))
# Regression: fixed-step RK4 previously went below this culture's capacity floor.
p <- c(k_per_od=5.4467557,Qmax_per_od=882.0815)
a <- quota_predict_state(times,500,0.1657396,p)
stopifnot(min(a$B)>=500-0.1657396*p[['Qmax_per_od']]-1e-8)
# Infinite-capacity limit approaches first-order uptake.
p <- c(k_per_od=0.05,Qmax_per_od=1e12)
a <- quota_predict_state(times,500,0.2,p)
stopifnot(max(abs(a$B-500*exp(-0.05*0.2*times)))<1e-5)
# Invalid inputs must be rejected instead of silently clipping them.
for (call in list(
  quote(quota_predict_state(-1,500,0.5,p)),
  quote(quota_predict_state(1,-1,0.5,p)),
  quote(quota_predict_state(1,500,-0.5,p)),
  quote(quota_predict_state(1,500,0.5,c(k_per_od=1,Qmax_per_od=0))),
  quote(quota_predict_state(1,500,0.5,c(k_per_od=1,Qmax_per_od=1),q0=2)))) {
  stopifnot(inherits(tryCatch(eval(call),error=identity),'error'))
}
cat('PASS:',nrow(grid),'quota parameter combinations; analytical/adaptive agreement,',
    'capacity, conservation, limiting cases and invalid-input guards. Max B difference:',max_error,'pM.\n')
