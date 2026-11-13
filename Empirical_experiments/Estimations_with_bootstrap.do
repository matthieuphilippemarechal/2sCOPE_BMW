* Estimations_with_bootstrap.do
* - The four estimators (OLS, 2sCOPE, BMW and IV via control function) are estimated in the SAME
*   replications: each replication repeats the whole procedure (ranks, first and second stage).
*   This gives the standard error of the differences between methods (contrasts).
* - If the dataset has a variable "cluster" (BLP, airfare), the bootstrap resamples clusters.
* - Three specifications are reported for each dataset (column "especificacion" of the results file):
*   "lineal"    : linear model.
*   "con z2"    : z^2 added to the equation of y. The control term is built for z only (z^2 does not
*                 enter the first stage). The effect of z at the mean is reported,
*                 efz = b(z) + 2*b(z2)*mean(z).
*   "cambio z2" : in each bootstrap replication BOTH specifications are estimated on the same resampled
*                 data, and the change efz(z2) - b(z)(linear) is computed for each method, as well as the
*                 change in rho. The standard error of the change thus accounts for the covariance
*                 between both specifications. Under the assumptions of each method, adding z^2 should
*                 not change the correction: a large t(change) indicates that the correction is not
*                 robust to the functional form. Contrasts of the change between methods are also
*                 reported (e.g. change 2sCOPE - change IV).
*   If z is binary (MEPS), the specifications with z^2 are skipped (z^2 = z).
* - All bootstraps use the same seed, so the replications draw the same resamples in all specifications.
* Results are saved in Results/Resultados_Estimaciones_bootstrap.xlsx (one sheet per dataset).

* Paths (relative to the root folder of the repository: "cd" to that folder before running this do-file)
local dir_db  "Cleaned databases"   // input datasets (.dta)
local dir_res "Results"             // results (.xlsx)
capture mkdir "`dir_res'"

local B = 1000              // number of bootstrap replications
local semilla = 20261110    // seed for reproducibility

* Auxiliary program: builds the coefficient vector of a method with readable names
capture program drop armar_b
program define armar_b, rclass
    args eq cf z2 zbar
    tempname m
    matrix `m' = e(b)
    local nombres : colnames e(b)
    if "`z2'" != "" {
        matrix `m' = `m', J(1, 1, _b[z] + 2*_b[`z2']*`zbar')
        local nombres `nombres' efz
    }
    if "`cf'" != "" local nombres : subinstr local nombres "`cf'" "rho", word
    if "`z2'" != "" local nombres : subinstr local nombres "`z2'" "z2", word
    matrix colnames `m' = `nombres'
    matrix coleq `m' = `eq'
    return matrix b = `m'
end

* Program that computes OLS, 2sCOPE, BMW and IV (control function) on the data in memory
capture program drop ifcf
program define ifcf, eclass
    syntax , x(varlist) [iv(varlist) cuadratico]
    tempname b1 b2 b3 b4 b
    tempvar z2 rz etaz cf2 e rb cfb v touse

    local z2var ""
    if "`cuadratico'" != "" {
        quietly generate double `z2' = z^2
        local z2var `z2'
    }
    quietly summarize z, meanonly
    local zbar = r(mean)

    // OLS
    quietly regress y z `z2var' `x'
    armar_b MCO "" "`z2var'" `zbar'
    matrix `b1' = r(b)

    // 2sCOPE
    local eta_controles_exogenos ""
    foreach var of local x {
        tempvar r_aux eta_aux
        quietly egen `r_aux' = rank(`var')
        quietly generate `eta_aux' = invnormal(`r_aux' / (_N + 1))
        local eta_controles_exogenos `eta_controles_exogenos' `eta_aux'
    }
    quietly egen `rz' = rank(z)
    quietly generate `etaz' = invnormal(`rz' / (_N + 1))
    quietly regress `etaz' `eta_controles_exogenos'
    quietly predict `cf2', residual
    quietly regress y z `z2var' `x' `cf2'
    armar_b tsCOPE "`cf2'" "`z2var'" `zbar'
    matrix `b2' = r(b)

    // BMW
    quietly regress z `x'
    quietly predict `e', residual
    quietly egen `rb' = rank(`e')
    quietly generate `cfb' = invnormal(`rb' / (_N + 1))
    quietly regress y z `z2var' `x' `cfb'
    armar_b BMW "`cfb'" "`z2var'" `zbar'
    matrix `b3' = r(b)

    // IV via control function (in the linear model, same coefficient as 2SLS)
    if "`iv'" != "" {
        quietly regress z `x' `iv'
        quietly predict `v', residual
        quietly regress y z `z2var' `x' `v'
        armar_b IV "`v'" "`z2var'" `zbar'
        matrix `b4' = r(b)
        matrix `b' = `b1', `b2', `b3', `b4'
    }
    else {
        matrix `b' = `b1', `b2', `b3'
    }
    quietly generate byte `touse' = 1
    ereturn post `b', esample(`touse')
end

* Program that estimates, on the same data, the linear model and the model with z^2, and returns the change
capture program drop ifcf_cambio
program define ifcf_cambio, eclass
    syntax , x(varlist) [iv(varlist)]
    tempname bl bq out sl sq
    local opcion_iv ""
    if "`iv'" != "" local opcion_iv "iv(`iv')"

    ifcf, x(`x') `opcion_iv'
    matrix `bl' = e(b)
    ifcf, x(`x') `opcion_iv' cuadratico
    matrix `bq' = e(b)

    local metodos MCO tsCOPE BMW
    if "`iv'" != "" local metodos `metodos' IV
    local nombres ""
    foreach m of local metodos {
        // effect of z: linear = b(z); with z^2 = efz
        scalar `sl' = el(`bl', 1, colnumb(`bl', "`m':z"))
        scalar `sq' = el(`bq', 1, colnumb(`bq', "`m':efz"))
        matrix `out' = nullmat(`out'), (`sl', `sq', `sq' - `sl')
        local nombres `nombres' `m':z_lineal `m':efz_z2 `m':cambio
        // coefficient of the control function
        if "`m'" != "MCO" {
            scalar `sl' = el(`bl', 1, colnumb(`bl', "`m':rho"))
            scalar `sq' = el(`bq', 1, colnumb(`bq', "`m':rho"))
            matrix `out' = `out', (`sl', `sq', `sq' - `sl')
            local nombres `nombres' `m':rho_lineal `m':rho_z2 `m':cambio_rho
        }
    }
    matrix colnames `out' = `nombres'
    tempvar touse
    quietly generate byte `touse' = 1
    ereturn post `out', esample(`touse')
end

* Program that saves coefficients, bootstrap standard errors and t-statistics in the results file
capture program drop guardar
program define guardar
    args db especificacion
    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local nombres : colnames `b'
    local ecuaciones : coleq `b'
    local k = 0
    foreach v of local nombres {
        local ++k
        local metodo : word `k' of `ecuaciones'
        if "`metodo'" == "tsCOPE" local metodo "2sCOPE"
        local se = sqrt(`V'[`k',`k'])
        post resultados ("`db'") ("`especificacion'") ("`metodo'") ("`v'") (`b'[1,`k']) (`se') (`b'[1,`k']/`se')
    }
end

* Program that saves contrasts between methods (differences with their bootstrap standard error)
capture program drop contrastes
program define contrastes
    args db especificacion p hay_iv
    local lista "tsCOPE-MCO BMW-MCO tsCOPE-BMW"
    if "`hay_iv'" == "1" local lista "`lista' IV-MCO tsCOPE-IV BMW-IV"
    foreach c of local lista {
        gettoken m1 m2 : c, parse("-")
        local m2 : subinstr local m2 "-" ""
        quietly lincom [`m1']`p' - [`m2']`p'
        local etiqueta "`m1' - `m2'"
        local etiqueta : subinstr local etiqueta "tsCOPE" "2sCOPE", all
        post resultados ("`db'") ("`especificacion'") ("Contraste") ("`etiqueta' (`p')") (r(estimate)) (r(se)) (r(estimate)/r(se))
    }
end

capture postclose resultados
tempfile res
postfile resultados str40 db str10 especificacion str10 metodo str32 variable double(coef se t) using `res', replace

local archivos : dir "`dir_db'" files "*.dta"
local list_db ""
foreach archivo of local archivos {
    local nombre = subinstr("`archivo'", ".dta", "", .)
    local list_db `list_db' `nombre'
}

foreach db of local list_db {
    display "`db'"
    use "`dir_db'/`db'", clear

    capture unab controles_exogenos : X*

    local variables_instrumentales ""
    capture unab variables_instrumentales : iv*
    local opcion_iv ""
    local hay_iv = 0
    if "`variables_instrumentales'" != "" {
        display "Instruments in `db': `variables_instrumentales'"
        local opcion_iv "iv(`variables_instrumentales')"
        local hay_iv = 1
    }
    else {
        display "Instruments in `db': no instrumental variables found"
    }

    local opcion_cluster ""
    capture confirm variable cluster
    if _rc == 0 {
        local opcion_cluster "cluster(cluster)"
    }

    // Linear model
    bootstrap _b, reps(`B') seed(`semilla') `opcion_cluster' nodots: ifcf, x(`controles_exogenos') `opcion_iv'
    guardar "`db'" "lineal"
    contrastes "`db'" "lineal" z `hay_iv'

    // Check with z^2 (skipped if z is binary)
    quietly count if !inlist(z, 0, 1)
    if r(N) > 0 {
        bootstrap _b, reps(`B') seed(`semilla') `opcion_cluster' nodots: ifcf, x(`controles_exogenos') `opcion_iv' cuadratico
        guardar "`db'" "con z2"
        contrastes "`db'" "con z2" efz `hay_iv'

        // Change when z^2 is added: both specifications in the same replication
        bootstrap _b, reps(`B') seed(`semilla') `opcion_cluster' nodots: ifcf_cambio, x(`controles_exogenos') `opcion_iv'
        guardar "`db'" "cambio z2"
        contrastes "`db'" "cambio z2" cambio `hay_iv'
    }
    else {
        display "`db': z is binary, the check with z^2 is skipped"
    }
}

postclose resultados
use `res', clear
* One sheet per dataset
capture erase "`dir_res'/Resultados_Estimaciones_bootstrap.xlsx"
levelsof db, local(bases)
foreach b of local bases {
    preserve
    keep if db == "`b'"
    drop db
    export excel using "`dir_res'/Resultados_Estimaciones_bootstrap.xlsx", sheet("`b'") firstrow(variables) sheetreplace
    restore
}
