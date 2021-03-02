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

subroutine dqk31(f,a,b,reslt,abserr,resabs,resasc)
!***begin prologue  dqk31
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  31-point gauss-kronrod rules
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
!                       result is computed by applying the 31-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 15-point gauss rule (resg).
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
!***end prologue  dqk31

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
real(kind=wp) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(15), fv2(15), hlgth, resg, resk, reskh, uflow
integer(kind=iwp) :: j, jtw, jtwm1
real(kind=wp), external :: d1mach

! the abscissae and weights are given for the interval (-1,1).
! because of symmetry only the positive abscissae and their
! corresponding weights are given.
!
! xgk    - abscissae of the 31-point kronrod rule
!          xgk(2), xgk(4), ...  abscissae of the 15-point
!          gauss rule
!          xgk(1), xgk(3), ...  abscissae which are optimally
!          added to the 15-point gauss rule
!
! wgk    - weights of the 31-point kronrod rule
!
! wg     - weights of the 15-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

real(kind=wp), parameter :: wg(8) = [0.030753241996117268354628393577204_wp, &
                                     0.070366047488108124709267416450667_wp, &
                                     0.107159220467171935011869546685869_wp, &
                                     0.139570677926154314447804794511028_wp, &
                                     0.166269205816993933553200860481209_wp, &
                                     0.186161000015562211026800561866423_wp, &
                                     0.198431485327111576456118326443839_wp, &
                                     0.202578241925561272880620199967519_wp]
real(kind=wp), parameter :: xgk(16) = [0.998002298693397060285172840152271_wp, &
                                       0.987992518020485428489565718586613_wp, &
                                       0.967739075679139134257347978784337_wp, &
                                       0.937273392400705904307758947710209_wp, &
                                       0.897264532344081900882509656454496_wp, &
                                       0.848206583410427216200648320774217_wp, &
                                       0.790418501442465932967649294817947_wp, &
                                       0.724417731360170047416186054613938_wp, &
                                       0.650996741297416970533735895313275_wp, &
                                       0.570972172608538847537226737253911_wp, &
                                       0.485081863640239680693655740232351_wp, &
                                       0.394151347077563369897207370981045_wp, &
                                       0.299180007153168812166780024266389_wp, &
                                       0.201194093997434522300628303394596_wp, &
                                       0.101142066918717499027074231447392_wp, &
                                       0.000000000000000000000000000000000_wp]
real(kind=wp), parameter :: wgk(16) = [0.005377479872923348987792051430128_wp, &
                                       0.015007947329316122538374763075807_wp, &
                                       0.025460847326715320186874001019653_wp, &
                                       0.035346360791375846222037948478360_wp, &
                                       0.044589751324764876608227299373280_wp, &
                                       0.053481524690928087265343147239430_wp, &
                                       0.062009567800670640285139230960803_wp, &
                                       0.069854121318728258709520077099147_wp, &
                                       0.076849680757720378894432777482659_wp, &
                                       0.083080502823133021038289247286104_wp, &
                                       0.088564443056211770647275443693774_wp, &
                                       0.093126598170825321225486872747346_wp, &
                                       0.096642726983623678505179907627589_wp, &
                                       0.099173598721791959332393173484603_wp, &
                                       0.100769845523875595044946662617570_wp, &
                                       0.101330007014791549017374792767493_wp]

! list of major variables
! -----------------------
!
! centr  - mid point of the interval
! hlgth  - half-length of the interval
! absc   - abscissa
! fval*  - function value
! resg   - result of the 15-point gauss formula
! resk   - result of the 31-point kronrod formula
! reskh  - approximation to the mean value of f over (a,b),
!          i.e. to i/(b-a)
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! uflow is the smallest positive magnitude.
!
!***first executable statement  dqk31

epmach = d1mach(4)
uflow = d1mach(1)

centr = Half*(a+b)
hlgth = Half*(b-a)
dhlgth = abs(hlgth)

! compute the 31-point kronrod approximation to
! the integral, and estimate the absolute error.

fc = f(centr)
resg = fc*wg(8)
resk = fc*wgk(16)
resabs = abs(resk)
do j=1,7
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
do j=1,8
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
resasc = wgk(16)*abs(fc-reskh)
do j=1,15
  resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
end do
reslt = resk*hlgth
resabs = resabs*dhlgth
resasc = resasc*dhlgth
abserr = abs((resk-resg)*hlgth)
if ((resasc /= Zero) .and. (abserr /= Zero)) abserr = resasc*min(One,(200.0_wp*abserr/resasc)**OneHalf)
if (resabs > uflow/(50.0_wp*epmach)) abserr = max((epmach*50.0_wp)*resabs,abserr)

return

end subroutine dqk31
