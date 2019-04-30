/*
Version 1.0.0 (April 6, 2019)
Author:  Moritz Marbach, moritz.marbach@gess.ethz.ch
*/

program ivdesc,  rclass
	version 12

	syntax varlist(default=none min=1 max=3) [if/] [in/] [, NOboot VARiance fmt(passthru) Reps(integer 1000)]

	local X: word 1 of `varlist'
	local D: word 2 of `varlist'
	local Z: word 3 of `varlist'

	preserve 

		if !missing("`if'") {
			quietly: keep if `if' & `X'!=. & `D'!=. & `Z'!=.
		}
		else {
			quietly: keep if `X'!=. & `D'!=. & `Z'!=. 
		}

		if missing("`noboot'") {
			
			ivdesc_calc `X' `D' `Z', `variance'

			local mu = r(mu)
			local mu_co = r(mu_co)
			local mu_nt = r(mu_nt)
			local mu_at = r(mu_at)

			local pi_co = r(pi_co)
			local pi_nt = r(pi_nt)
			local pi_at = r(pi_at)

			if !missing("`variance'") {

				local v_co = r(v_co)
				local v_nt = r(v_nt)
				local v_at = r(v_at)
				local v = r(v)

			}

			bootstrap mu_co=r(mu_co) mu_at=r(mu_at) mu_nt=r(mu_nt) mu=r(mu) /// 
								pi_co=r(pi_co) pi_at=r(pi_at) pi_nt=r(pi_nt), reps(`reps') notable nolegend nowarn noheader: ivdesc_calc `X' `D' `Z' 
			
			matrix bootse = e(se)

			local se_mu = bootse[1,colnumb(bootse,"mu")]
			local se_mu_co = bootse[1,colnumb(bootse,"mu_co")]
			local se_mu_at = bootse[1,colnumb(bootse,"mu_at")]
			local se_mu_nt = bootse[1,colnumb(bootse,"mu_nt")]

			local se_pi_co = bootse[1,colnumb(bootse,"pi_co")]
			local se_pi_at = bootse[1,colnumb(bootse,"pi_at")]
			local se_pi_nt = bootse[1,colnumb(bootse,"pi_nt")]

			if missing("`variance'") {

				matrix input res = ( `mu', `se_mu', 1, 0 \ `mu_co', `se_mu_co', `pi_co', `se_pi_co' \ `mu_at', `se_mu_at', `pi_at', `se_pi_at' \ `mu_nt', `se_mu_nt', `pi_nt', `se_pi_nt' )
				matrix colnames res = "Mean" "Boot.-SE" "Proportion" "Boot.-SE"

			} 

			else {

				matrix input res = ( `mu', `se_mu', `v', 1, 0 \ `mu_co', `se_mu_co', `v_co', `pi_co', `se_pi_co' \ `mu_at', `se_mu_at', `v_at',  `pi_at', `se_pi_at' \ `mu_nt', `se_mu_nt', `v_nt',  `pi_nt', `se_pi_nt' )
				matrix colnames res = "Mean" "Boot.-SE" "Variance" "Proportion" "Boot.-SE"

			}

		} 

		else {

			ivdesc_calc `X' `D' `Z', `variance'
	
			if missing("`variance'") {

				matrix input res = ( `r(mu)', `r(se_mu)', 1, 0 \ `r(mu_co)', ., `r(pi_co)', `r(se_pi_co)' \ `r(mu_nt)', `r(se_mu_nt)', `r(pi_nt)', `r(se_pi_nt)' \ `r(mu_at)', `r(se_mu_at)', `r(pi_at)', `r(se_pi_at)'  )
				matrix colnames res = "Mean" "SE" "Proportion" "SE"

			} 

			else {

				matrix input res = ( `r(mu)', `r(se_mu)', `r(v)', 1, 0 \ `r(mu_co)', ., `r(v_co)', `r(pi_co)', `r(se_pi_co)' \ `r(mu_nt)', `r(se_mu_nt)', `r(v_nt)', `r(pi_nt)', `r(se_pi_nt)' \ `r(mu_at)', `r(se_mu_at)', `r(v_at)',  `r(pi_at)', `r(se_pi_at)'  )
				matrix colnames res = "Mean" "SE" "Variance" "Proportion" "SE"

			}

		}

		matrix rownames res = "whole sample" "complier" "never-taker" "always-taker"
			
		estout matrix(res, `fmt'), mlabels(none) title("Variable: " `X')

		return mat ivdesc = res

	restore

end



