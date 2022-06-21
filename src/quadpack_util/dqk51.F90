!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine dqk51(f,a,b,reslt,abserr,resabs,resasc)
!***begin prologue  dqk51
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  51-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
!                           estimate
!                       j = integral of abs(f) over (a,b)
!***description
!
!           integration rules
!           standard fortran subroutine
!           real*8 version
!
!           parameters
!            on entry
!              f      - real*8
!                       function subprogram defining the integrand
!                       function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real*8
!                       lower limit of integration
!
!              b      - real*8
!                       upper limit of integration
!
!            on return
!              reslt  - real*8
!                       approximation to the integral i
!                       result is computed by applying the 51-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 25-point gauss rule (resg).
!
!              abserr - real*8
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-reslt)
!
!              resabs - real*8
!                       approximation to the integral j
!
!              resasc - real*8
!                       approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  dqk51

use Constants, only: Zero, One, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
interface
  function f(x)
    import :: wp
    real(kind=wp) :: f
    real(kind=wp), intent(in) :: x
  end function f
end interface
real(kind=wp), intent(in) :: a, b
real(kind=wp), intent(out) :: reslt, abserr, resabs, resasc
real(kind=wp) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(25), fv2(25), hlgth, resg, resk, reskh, uflow
integer(kind=iwp) :: j, jtw, jtwm1
real(kind=wp), external :: d1mach

! the abscissae and weights are given for the interval (-1,1).
! because of symmetry only the positive abscissae and their
! corresponding weights are given.
!
! xgk    - abscissae of the 51-point kronrod rule
!          xgk(2), xgk(4), ...  abscissae of the 25-point
!          gauss rule
!          xgk(1), xgk(3), ...  abscissae which are optimally
!          added to the 25-point gauss rule
!
! wgk    - weights of the 51-point kronrod rule
!
! wg     - weights of the 25-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

real(kind=wp), parameter :: wg(13) = [0.011393798501026287947902964113235_wp, &
                                      0.026354986615032137261901815295299_wp, &
                                      0.040939156701306312655623487711646_wp, &
                                      0.054904695975835191925936891540473_wp, &
                                      0.068038333812356917207187185656708_wp, &
                                      0.080140700335001018013234959669111_wp, &
                                      0.091028261982963649811497220702892_wp, &
                                      0.100535949067050644202206890392686_wp, &
                                      0.108519624474263653116093957050117_wp, &
                                      0.114858259145711648339325545869556_wp, &
                                      0.119455763535784772228178126512901_wp, &
                                      0.122242442990310041688959518945852_wp, &
                                      0.123176053726715451203902873079050_wp]
real(kind=wp), parameter :: xgk(26) = [0.999262104992609834193457486540341_wp, &
                                       0.995556969790498097908784946893902_wp, &
                                       0.988035794534077247637331014577406_wp, &
                                       0.976663921459517511498315386479594_wp, &
                                       0.961614986425842512418130033660167_wp, &
                                       0.942974571228974339414011169658471_wp, &
                                       0.920747115281701561746346084546331_wp, &
                                       0.894991997878275368851042006782805_wp, &
                                       0.865847065293275595448996969588340_wp, &
                                       0.833442628760834001421021108693570_wp, &
                                       0.797873797998500059410410904994307_wp, &
                                       0.759259263037357630577282865204361_wp, &
                                       0.717766406813084388186654079773298_wp, &
                                       0.673566368473468364485120633247622_wp, &
                                       0.626810099010317412788122681624518_wp, &
                                       0.577662930241222967723689841612654_wp, &
                                       0.526325284334719182599623778158010_wp, &
                                       0.473002731445714960522182115009192_wp, &
                                       0.417885382193037748851814394594572_wp, &
                                       0.361172305809387837735821730127641_wp, &
                                       0.303089538931107830167478909980339_wp, &
                                       0.243866883720988432045190362797452_wp, &
                                       0.183718939421048892015969888759528_wp, &
                                       0.122864692610710396387359818808037_wp, &
                                       0.061544483005685078886546392366797_wp, &
                                       0.000000000000000000000000000000000_wp]
real(kind=wp), parameter :: wgk(26) = [0.001987383892330315926507851882843_wp, &
                                       0.005561932135356713758040236901066_wp, &
                                       0.009473973386174151607207710523655_wp, &
                                       0.013236229195571674813656405846976_wp, &
                                       0.016847817709128298231516667536336_wp, &
                                       0.020435371145882835456568292235939_wp, &
                                       0.024009945606953216220092489164881_wp, &
                                       0.027475317587851737802948455517811_wp, &
                                       0.030792300167387488891109020215229_wp, &
                                       0.034002130274329337836748795229551_wp, &
                                       0.037116271483415543560330625367620_wp, &
                                       0.040083825504032382074839284467076_wp, &
                                       0.042872845020170049476895792439495_wp, &
                                       0.045502913049921788909870584752660_wp, &
                                       0.047982537138836713906392255756915_wp, &
                                       0.050277679080715671963325259433440_wp, &
                                       0.052362885806407475864366712137873_wp, &
                                       0.054251129888545490144543370459876_wp, &
                                       0.055950811220412317308240686382747_wp, &
                                       0.057437116361567832853582693939506_wp, &
                                       0.058689680022394207961974175856788_wp, &
                                       0.059720340324174059979099291932562_wp, &
                                       0.060539455376045862945360267517565_wp, &
                                       0.061128509717053048305859030416293_wp, &
                                       0.061471189871425316661544131965264_wp, &
                                       ! note: wgk (26) was calculated from the values of wgk(1..25)
                                       0.061580818067832935078759824240066_wp]

! list of major variables
! -----------------------
!
! centr  - mid point of the interval
! hlgth  - half-length of the interval
! absc   - abscissa
! fval*  - function value
! resg   - result of the 25-point gauss formula
! resk   - result of the 51-point kronrod formula
! reskh  - approximation to the mean value of f over (a,b),
!          i.e. to i/(b-a)
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! uflow is the smallest positive magnitude.
!
!***first executable statement  dqk51

epmach = d1mach(4)
uflow = d1mach(1)

centr = Half*(a+b)
hlgth = Half*(b-a)
dhlgth = abs(hlgth)

! compute the 51-point kronrod approximation to
! the integral, and estimate the absolute error.

fc = f(centr)
resg = fc*wg(13)
resk = fc*wgk(26)
resabs = abs(resk)
do j=1,12
  jtw = j*2
  absc = hlgth*xgk(jtw)
  fval1 = f(centr-absc)
  fval2 = f(centr+absc)
  fv1(jtw) = fval1
  fv2(jtw) = fval2
  fsum = fval1+fval2
  resg = resg+wg(j)*fsum
  resk = resk+wgk(jtw)*fsum
  resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
end do
do j=1,13
  jtwm1 = j*2-1
  absc = hlgth*xgk(jtwm1)
  fval1 = f(centr-absc)
  fval2 = f(centr+absc)
  fv1(jtwm1) = fval1
  fv2(jtwm1) = fval2
  fsum = fval1+fval2
  resk = resk+wgk(jtwm1)*fsum
  resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
end do
reskh = resk*Half
resasc = wgk(26)*abs(fc-reskh)
do j=1,25
  resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
end do
reslt = resk*hlgth
resabs = resabs*dhlgth
resasc = resasc*dhlgth
abserr = abs((resk-resg)*hlgth)
if ((resasc /= Zero) .and. (abserr /= Zero)) abserr = resasc*min(One,(200.0_wp*abserr/resasc)**OneHalf)
if (resabs > uflow/(50.0_wp*epmach)) abserr = max((epmach*50.0_wp)*resabs,abserr)

return

end subroutine dqk51
