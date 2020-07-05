{smcl}
{* *! version 1.0.1  April 30, 2019 @ 18:24:26}{...}
{cmd:help ivdesc}
{hline}

{title:Title}

{p2colset 5 20 22 2}{...}  {p2col :{hi:ivdesc} {hline 2}}
 Profiling compliers and non-compliers for instrumental variable analysis {p_end} {p2colreset}{...}

{title:Syntax}

{p 8 17 2}{cmd:ivdesc} 
{it:{help varname:covariate}} 
{it:{help varname:treatment}} 
{it:{help varname:instrument}}
{ifin} 
[{cmd:,} {it:options}]


{marker options}{...}
{title:Options}

{phang}{opt no:boot} report analytical standard errors instead of bootstrap standard errors 

{phang}{opt r:eps(#)}  perform # bootstrap replications; default is reps(1000)

{phang}{opt var:iance} also report the variance of the covariate for each subgroup 

{phang}{opt nobal:ance} skip the balance test 

{phang}{cmd:fmt(}{it:subopts}{cmd:)} passed through to {cmd:estout matrix(,subopts)} that handles the display of the estimates


{marker description}{...}
{title:Description}

{phang} {cmd:ivdesc} estimates the mean and the associated standard error of a covariate for the complier, never-taker and always-taker subpopulation within a sample where some, but not all, units are encouraged by instrument to take the treatment. 

{title:Remarks}

{phang} Observations with missing values in either {it:covariate}, {it:treatment}, 
or {it:instrument} are deleted before estimation (listwise deletion). 

{phang} One-sided noncompliance is supported. The mean for the always-/never-taker subpopulation will only be computed if there are at least two observed units in these subpopulations. 

{phang} If {cmd:noboot}, analytical standard errors are calculated for the mean of the whole sample as well as the never-taker and always-taker subpopulation. For the complier subpopulation no analytical estimator for the standard error is available.  

{phang} The balance test is a t-test allowing for unequal variances. 


{title:Examples}

{pstd}Load the JTPA data 

{p 4 8 2}{stata "use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear":. use "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta",clear}{p_end}

{p 4 8 2}{stata "ivdesc age training assignmt":. ivdesc age training assignmt}{p_end}
{p 4 8 2}{stata "ivdesc hispanic training assignmt":. ivdesc hispanic training assignmt}{p_end}

{pstd}Plot the results

{p 4 8 2}{stata "matrix C = r(ivdesc)'":. matrix C = r(ivdesc)'}{p_end}
{p 4 8 2}{stata "coefplot matrix(C), se(C[2])":. coefplot matrix(C), se(C[2])}{p_end}


{title:Saved results}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(ivdesc)}} all estimates {p_end}

{p2col 5 15 19 2: Scalar}{p_end}
{synopt:{cmd:r(pval)}} p-value from balance test {p_end}


{title:Reference}
{p 4 8 2}

{pstd}Moritz Marbach and Dominik Hangartner. (2020). Profiling Compliers and Non-compliers for Instrumental Variable Analysis. {it:Political Analysis} 28(3), 435-444. {p_end}

{title:Authors}

	Moritz Marbach (Maintainer), moritz.marbach@gess.ethz.ch
	ETH Zurich, Switzerland
