{smcl}
{* *! version 0.1  30jul2022}{...}
{vieweralsosee "[R] cumul" "help cumul"}{...}
{viewerjumpto "Syntax" "latentcs##syntax"}{...}
{viewerjumpto "Description" "latentcs##description"}{...}
{viewerjumpto "Options" "latentcs##options"}{...}
{viewerjumpto "Stored results" "latentcs##results"}{...}
{viewerjumpto "Remarks" "latentcs##remarks"}{...}
{viewerjumpto "Examples" "latentcs##examples"}{...}
{viewerjumpto "Author" "latentcs##author"}{...}
{viewerjumpto "References" "latentcs##references"}{...}

{title:Title}

{phang}
{* phang is short for p 4 8 2}
{bf:latentcs} {hline 2}  Compare two latent distributions given ordinal data


{marker syntax}{...}
{title:Syntax}

{phang}
{cmd:latentcs}
{varname}
{ifin}
{cmd:,}
{opth "by(varname:groupvar)"}
[{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{* syntab:Main}{p2coldent :* {opth "by(varname:groupvar)"}}binary variable defining two groups{p_end}
{synopt :{opt sw:itch}}switch direction of comparison{p_end}
{* syntab:Reporting}{synopt :{opt l:evel(#)}}set confidence level; default is 90 or {cmd:c(level)} as set by {cmd:set level}{p_end}
{* syntab:Graphical}{synopt :{opt noplot}}suppress plot{p_end}
{* syntab:Advanced}{synopt :{opt seed(#)}}set random-number seed to #{p_end}
{synopt :{opt ndraw(#)}}set number of simulation draws to #; default is 10,000{p_end}
{synoptline}
{pstd}
* {opt by(groupvar)} is required.{p_end}
{pstd}
{bf:by} (the prefix) is allowed; see {help by} or {manlink D by}{p_end}
{p 4 6 2}{cmd:aweight}s and {cmd:fweight}s are allowed; see {help weight}.{p_end}

{pstd}
As with {help ksmirnov}, {it:{help varname:groupvar}} must have exactly two distinct values.
If it takes more, then {bf:{help if}} can be used.
For example, {input:latentcs y if g==0 | g==1 , by(g)}


{marker description}{...}
{title:Description}

{pstd}
{* pstd is short for p 4 4 2}
{cmd:latentcs} compares two latent distributions given two samples of ordinal data.
The {varname} is the ordinal variable for comparison, like self-reported health or level of agreement with a political statement.
The {it:{help varname:groupvar}} is a variable taking only two distinct values that defines two groups, like an indicator/dummy for "male" (equal to 1 for male, equal to 0 for non-male).
It is assumed the groups are sampled independently, and although weights are accounted for if provided, sampling is otherwise assumed iid from the two respective group population distributions.
(The general methodology from Kaplan and Zhao (2022) can allow non-iid sampling, but this particular implementation does not; see {help "latentcs##references":References}.)

{pstd}
Interest is in whether the latent distribution for the larger value of {it:groupvar} tends to have higher values than the distribution for the smaller value of {it:groupvar}.
Specifically, interest is in the set of tau for which the tau-quantile is higher.
This command stores and reports both the point estimate (our best guess) of that set of tau as well as an inner confidence set that is contained within the true set of tau with high probability, specifically the requested confidence level.
That is, there is a high probability the true set of tau is even larger than the inner confidence set; it is conservative, but provides more certainty.
These reported sets of tau are valid for the ordinal distribution comparison and can further be interpreted in terms of the underlying latent distributions as long as the reporting thresholds do not differ in a way that makes the smaller-{it:groupvar} values of {varname} higher; for details, see the formal assumptions and discussion in {help "latentcs##references":Kaplan and Zhao (2022)}.

{pstd}
The methodology was proposed by Kaplan and Zhao (2022); see {help "latentcs##references":References}.


{marker options}{...}
{title:Options}
{* dlgtab:Main}
{phang}
{opth "by(varname:groupvar)"} is required.
It specifies a binary variable that identifies the two groups whose distributions are compared.
Specifically, interest is in the set of tau values for which the tau-quantile of the larger-{it:groupvar} distribution of {varname} is higher than the tau-quantile of the smaller-{it:groupvar} distribution of {varname}.{p_end}

{phang}
{opt switch} reverses the direction of comparison, to find evidence of higher {varname} values given lower {it:{help varname:groupvar}}.{p_end}

{phang}
{opt level(#)}; see 
{helpb estimation options##level():[R] estimation options}.{p_end}

{phang}
{opt noplot} suppresses the plot.{p_end}

{phang}
{opt seed(#)} sets the random-number seed (as in {helpb bootstrap}).
The default is recommended to avoid the appearance of manipulation.
The current seed is restored upon completion of the command.
{p_end}

{phang}
{opt ndraw(#)} sets the number of random draws for simulating the critical value.
The default 10,000 should be sufficient for most purposes; smaller values can be used to speed computation for exploratory work, and larger values can be used for extra accuracy in results for publication.
{p_end}


{marker examples}{...}
{title:Examples}

{phang}{input:. webuse hanley , clear}{p_end}
{phang}{input:. latentcs rating , by(disease)}{p_end}
{pstd}
Additional examples are in the latentcs_examples.do file.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:latentcs} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(tau_bw_est)}}estimated ranges of tau for which the second group (or first group, with {opt switch} option) has larger tau-quantile (evidence of between-group inequality).
Each row is one such range (lower endpoint, upper endpoint).
Note: if none, then this is a 1-by-2 matrix with . as each entry.{p_end}
{synopt:{cmd:r(tau_bw_cs)}}similar to {cmd:r(tau_bw_est)} but for the more conservative inner confidence set of ranges, at confidence level {cmd:r(level)}.{p_end}
{synopt:{cmd:r(N)}}number of observations (overall, first group, second group){p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
R code, including code to replicate results from the paper, and a working paper version are available at {browse "https://kaplandm.github.io"}


{marker author}{...}
{title:Author}

{pstd}
David M. Kaplan{break}Department of Economics, University of Missouri{break}
kaplandm@missouri.edu{break}{browse "https://kaplandm.github.io"}

{marker references}
{title:References}

{pstd}
Kaplan, D. M., and Zhao, W. (2022).
Comparing latent inequality with ordinal data.
URL: {browse "https://kaplandm.github.io/#ordinal"}
{* browse "https://doi.org/10.1016/j.jeconom.2018.04.003"}
{* The Econometrics Journal, XXX(XXX):XXX-XXX.}
