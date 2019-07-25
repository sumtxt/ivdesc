/*
Version 1.0.1 (April 30, 2019)
Author:  Moritz Marbach, moritz.marbach@gess.ethz.ch
URL: https://github.com/sumtxt/ivdesc
*/

program ivdesc_calc,  rclass
	version 12 
	
	syntax varlist(default=none min=1 max=3) [, VARiance]

	local X: word 1 of `varlist'
	local D: word 2 of `varlist'
	local Z: word 3 of `varlist'

	* Estimation 

	quietly: { 

		summarize `D' if `Z'==1
		local pi_co1 = r(mean)
		local v_pi_co1 = r(Var)

		summarize `D' if `Z'==0
		local pi_co2 = r(mean)
		local v_pi_co2 = r(Var)

		summarize `Z'
		local N_z1 = r(sum)
		local N_z0 = r(N)-r(sum)

		local pi_co = `pi_co1'-`pi_co2'
		local se_pi_co = sqrt( (`v_pi_co1'/`N_z1') + (`v_pi_co2'/`N_z0') )

		if `pi_co'<0 {
			display as error "First-stage is negative. Please reverse coding of Z."
    	exit 42
    }

		summarize `D' if `Z'==1
		local pi_nt = 1-r(mean)
		local se_pi_nt = sqrt( r(Var)/`N_z1' )

		summarize `D' if `Z'==0 
		local pi_at = r(mean)
		local se_pi_at = sqrt( r(Var)/`N_z0' )

		summarize `X' if `Z'==1 & `D'==0
		local v_nt = r(Var)
		local mu_nt = r(mean)

		summarize `X' if `Z'==0 & `D'==1
		local v_at = r(Var)
		local mu_at = r(mean)

		summarize `X' 
		local mu = r(mean)
		local v = r(Var)

		count if `Z'==1 & `D'==0
		local k_nt = r(N) 

		count if `Z'==0 & `D'==1
		local k_at = r(N) 

		if `k_at'<2 {

			local mu_co = (1/`pi_co') * `mu' - (`pi_nt'/`pi_co') * `mu_nt' 

			if !missing("`variance'") {

				local v_co1 = (`v_nt' * `pi_nt') 

				local v_co2 = ( `mu_co' * `mu_co' * `pi_co' * (1-`pi_co') ) + ( `mu_nt' *`mu_nt' * `pi_nt' * (1-`pi_nt') ) 

				local v_co3 = (`mu_nt' * `pi_nt' * `mu_co' * `pi_co') 

				local v_co = (1/`pi_co')*`v' - (1/`pi_co')*(`v_co1' + `v_co2' -2 * `v_co3')

			}

		}
		
		else if `k_nt'<2 {

			local mu_co = (1/`pi_co') * `mu' - (`pi_at'/`pi_co') * `mu_at'

			if !missing("`variance'") {

				local v_co1 = (`v_at' * `pi_at')

				local v_co2 = ( `mu_co' * `mu_co' * `pi_co' * (1-`pi_co') ) + ( `mu_at' * `mu_at' * `pi_at' * (1-`pi_at') ) 

				local v_co3 = (`mu_at' * `pi_at' * `mu_co' * `pi_co')

				local v_co = (1/`pi_co')*`v' - (1/`pi_co')*(`v_co1' + `v_co2' -2 * `v_co3')

			}
		
		} 

		else {

			local mu_co = (1/`pi_co') * `mu' - (`pi_nt'/`pi_co') * `mu_nt' - (`pi_at'/`pi_co') * `mu_at'

			if !missing("`variance'") {

				local v_co1 = (`v_nt' * `pi_nt') + (`v_at' * `pi_at')

				local v_co2 = ( `mu_co' * `mu_co' * `pi_co' * (1-`pi_co') ) + ( `mu_nt' *`mu_nt' * `pi_nt' * (1-`pi_nt') ) + ( `mu_at' * `mu_at' * `pi_at' * (1-`pi_at') ) 

				local v_co3 = (`mu_nt' * `pi_nt' * `mu_co' * `pi_co') + (`mu_at' * `pi_at' * `mu_co' * `pi_co') + (`mu_at' * `pi_at' * `mu_nt' * `pi_nt') 

				local v_co = (1/`pi_co')*`v' - (1/`pi_co')*(`v_co1' + `v_co2' -2 * `v_co3')

			}

		}

		local se_mu = sqrt(`v')/sqrt(_N)
		local se_mu_nt = sqrt(`v_nt')/sqrt(`k_nt')
		local se_mu_at = sqrt(`v_at')/sqrt(`k_at')
	
	}

	return scalar mu = `mu'

	return scalar mu_nt = `mu_nt'
	return scalar mu_at = `mu_at'
	return scalar mu_co = `mu_co'

	return scalar se_mu = `se_mu'
	return scalar se_mu_nt = `se_mu_nt'
	return scalar se_mu_at = `se_mu_at'

	return scalar pi_co = `pi_co'
	return scalar pi_nt = `pi_nt'
	return scalar pi_at = `pi_at'
	
	return scalar se_pi_co = `se_pi_co'
	return scalar se_pi_nt = `se_pi_nt'
	return scalar se_pi_at = `se_pi_at'
	
	if !missing("`variance'") {

		return scalar v_nt = `v_nt'
		return scalar v_at = `v_at'
		return scalar v_co = `v_co'
		return scalar v = `v'

	}	

end
