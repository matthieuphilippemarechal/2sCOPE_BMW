local p = 1000 
local sample_sizes 250 500 1000 2000 4000

local beta = 0.7
local gamma_ = 0.4
local delta = 1
local sigma_x = 3        // x ~ N(0, sigma_x^2)
local lambda_e = 1/3     // e ~ Exp(lambda_e)

set matsize 8000
tempfile acumulado

local primero = 1
foreach nb of local sample_sizes {
    set seed 20261109                       
    matrix results = J(7 * `p', 12, .)

    forvalues k = 1/7 {
        use corr_z_u_versus_rho, clear
        quietly mean rho if _n == `k'
        matrix b = e(b)
        local rho = b[1,1]
        quietly mean corr_z_u if _n == `k'
        matrix b = e(b)
        local corr_z_u_teorico = b[1,1]

        forvalues i = 1/`p' {
            local fila = `i' + (`k'-1)*`p'
            matrix results[`fila',1]  = `rho'
            matrix results[`fila',2]  = `corr_z_u_teorico'
            matrix results[`fila',12] = `nb'
            clear
            quietly set obs `nb'

            generate x = `sigma_x' * rnormal()
            generate e = -(1 / `lambda_e')*log(1-runiform(0,1))
            generate z = `delta' * x + e
            generate eta = invnormal(1 - exp(-`lambda_e' * e))
            generate eps = rnormal()
            generate u = `rho' * eta + eps
            generate y = `beta' * x + `gamma_' * z + u

            quietly correlate z u
            matrix C = r(C)
            matrix results[`fila',3] = C[1,2]

            quietly regress z x
            matrix b = e(b)
            generate e_est = z - b[1,1] * x - b[1,2]
            sort e_est
            generate eta_e = invnormal(_n / (_N + 1))
            generate id = _n
            sort z
            generate eta_z = invnormal(_n / (_N + 1))
            sort x
            generate eta_x = invnormal(_n / (_N + 1))
            quietly regress eta_z eta_x
            matrix b = e(b)
            generate e_est_star = eta_z - b[1,1] * eta_x - b[1,2]
            sort id

            quietly regress y x z
            matrix b = e(b)
            matrix results[`fila',4] = b[1,1]
            matrix results[`fila',5] = b[1,2]

            quietly regress y x z e_est_star
            matrix b = e(b)
            matrix results[`fila',6] = b[1,1]
            matrix results[`fila',7] = b[1,2]
            matrix results[`fila',8] = b[1,3]

            quietly regress y x z eta_e
            matrix b = e(b)
            matrix results[`fila',9]  = b[1,1]
            matrix results[`fila',10] = b[1,2]
            matrix results[`fila',11] = b[1,3]
        }
    }

    clear
    svmat results
    if `primero' == 0 append using `acumulado'
    save `acumulado', replace
    local primero = 0
    display "Listo n = `nb'"
}

use `acumulado', clear
rename results1 rho
rename results2 corr_z_u_teorico
rename results3 corr_z_u
rename results4 beta_sin_correccion
rename results5 gamma_sin_correccion
rename results6 beta_con_cupula1
rename results7 gamma_con_cupula1
rename results8 rho_con_cupula1
rename results9 beta_con_cupula2
rename results10 gamma_con_cupula2
rename results11 rho_con_cupula2
rename results12 N
sort N rho, stable
save simulations_results_under_BMW_hypothesis, replace
