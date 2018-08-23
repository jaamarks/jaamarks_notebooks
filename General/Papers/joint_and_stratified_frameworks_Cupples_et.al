An Empirical Comparison of Joint and Stratified Frameworks for Studying GxE Interactions: Systolic Blood Pressure and Smoking in the CHARGE Gene-Lifestyle Interactions Working Group

A link to this paper is https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911246/


#==========================================================================
This paper discusses joint 2df vs stratified framework.

The main statistical framework for studying a gene-interaction models is the
joint framework which couples the interaction effects, such as GxE (or GxSex), along with
the genetic main.

There is an alternative framework called the stratified framework. This 
combines results from genetic main-effect analyses carried out separately
within exposed and unexposed groups. 

This paper discussed how a joint framework, like the 2df test that
couples the genetic-main interaction with the GxE interaction. It
details that this framework is useful to identify variants with low
main effect and moderate interaction effects. 

There is also the stratified framework. This has emerged for dichotomous
(binary) exposure variables. Subjects are stratified to exposed and unexposed
groupls. Genetic main is performed separately in each stratum. The results
of the stratum-specific genetic effects are subsequently combined to perform a 1df test
or a 2df test. Note that the stratified framework approximates the joint framework.
There are main-effect models readily available in many software packages and
therefore it is easier to implement in large-scale consortium setting.

Note, that I have a project for which this stratified framework would be
applicable to. I have been working with the UHS123 cohort and I could
split this cohort up into subjects that are male and subjects that are
female because this is a dichotomous variable that I am interested in
in my model. In fact, Dana Hancock suggested this in GitHub issue 97 and
actually has a joint-to-stratified calculation.
