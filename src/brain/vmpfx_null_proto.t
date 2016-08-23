#
# Input Specification
#

data.root           ****DATA_ROOT****
data.path           ****DATA_PATH****
<arg> data.root     root
<arg> data.path     path


data.file           $<data.path>$<data.root>
data.format         pgh

<arg> data.file   data
<arg> data.format format
<require> data.file


#
# Experimental Paradigm
#

conditions               1
design.images            1
design.aligned           1
fixed                    0
IAI                      ****IAI****

<data> design
# condition  stimulus.start   stimulus.len 
0 0 ****NUM_IMAGES****
<end>


#
# Fixed Parameters
#

shape.params              4
drift.degree              3
paste.form                0
paste.scale               0.001

<arg> drift.degree degree

#
# Prior Specification
#

prior.baseline          2000.0      1.0e-6
prior.drift               0.001     1.0      4.0
prior.resp                0.005    50.0
prior.sigmasq             1.8     200.0

<data> prior.shape 
        2.0  4.30
        4.0  1.0
        2.0  4.30
        4.0  0.43
<end>

prior.drift.scaled      false


#
# Algorithm Parameters --- vmpfx
#

vmpfx.iterations     2500
vmpfx.evaluations    25000
vmpfx.scd.tol         1.0e-5
vmpfx.cov.method        analytic
vmpfx.knots             ****NUM_KNOTS****
vmpfx.fix.knots         true

vmpfx.init           auto   
vmpfx.produce        residuals
vmpfx.model          full


#
# Algorithm Parameters --- vfpfx
#

vfpfx.samples       10500
vfpfx.burnin          500  
vfpfx.increment         1

vfpfx.init          auto

<data> vfpfx.jump
#         name              default      [value]
          baseline          0.5
          shape             0.05
          responsiveness    0.001
<end>

vfpfx.se.expand         2.5

vfpfx.knots.max         6
vfpfx.knots.init        $<vmpfx.knots>
vfpfx.fix.knots         $<vmpfx.fix.knots>

<data> vfpfx.products
#      name               class          which
       baseline           cumulants      2
       responsiveness     samples        1
       drift              samples        25
       shape              quantiles      5
       noise.precision    quantiles      5
<end>


#
# Output Format
#
    
output.root           vmpfx_results
output.path           ""
<arg> output.root     oroot
<arg> output.path     opath

output.file           $<output.path>$<output.root>.out

<cond> output.append $<0>
      vmpfx     false
      vfpfx     true
      <default> true
<end>

<arg> output.file    output
<arg> output.logfile logfile

