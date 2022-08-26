from os import system

#tap_dir  = '/Users/ocots/Boulot/recherche/Logiciels/dev/hampath/hampath308/install/tapenade3.8'
tap_dir  = '/net/home/c/ct705827/Documents/tapenade/tapenade3.8'
TAP      = tap_dir+'/bin/tapenade'
tap_home = tap_dir

build_dir='./'
#OPTT     = '-tangent -fixinterface -inputlanguage fortran90 -outputlanguage fortran90 -nooptim activity -O ' + build_dir
OPTT     = '-tangent -fixinterface -inputlanguage fortran90 -outputlanguage fortran90 -O ' + build_dir

system('TAPENADE_HOME=' + tap_home + ' ' + TAP + ' ' + OPTT + ' ' + '-tgtfuncname _da' + ' ' + '-head "hfunalpha(M, deltaI, pdeltaI, alpha)>(h)" -o hfunalpha hfunalpha.f90 ')
system('TAPENADE_HOME=' + tap_home + ' ' + TAP + ' ' + OPTT + ' ' + '-tgtfuncname _dc' + ' ' + '-head "hfuncone(M, pdeltaI)>(hcone)" -o hfuncone hfunalpha.f90')
system('TAPENADE_HOME=' + tap_home + ' ' + TAP + ' ' + OPTT + ' ' + '-tgtfuncname _da' + ' ' + '-head "hfunalpha_da(alpha)>(hd)" -o hfunalpha_da hfunalpha_da.f90')
system('TAPENADE_HOME=' + tap_home + ' ' + TAP + ' ' + OPTT + ' ' + '-tgtfuncname _di' + ' ' + '-head "hfunimplicit(M, q, p, cont)>(h)" -o hfunimplicit hfunimplicit.f90 hfunalpha.f90 hfunalpha_da.f90 hfunalpha_da_da.f90 hfuncone_dc.f90')
system('TAPENADE_HOME=' + tap_home + ' ' + TAP + ' ' + OPTT + ' ' + '-tgtfuncname _di' + ' ' + '-head "hfunimplicit_di(M, q, p, cont)>(hd)" -o hfunimplicit_di hfunimplicit_di.f90')

