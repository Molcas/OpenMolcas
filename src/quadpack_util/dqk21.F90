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

subroutine dqk21(f,a,b,reslt,abserr,resabs,resasc)
!***begin prologue  dqk21
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  21-point gauss-kronrod rules
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
!                       result is computed by applying the 21-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 10-point gauss rule (resg).
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
!***end prologue  dqk21

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
real(kind=wp) :: absc, centr, dhlgth, epmach, fc, fsum, fval1, fval2, fv1(10), fv2(10), hlgth, resg, resk, reskh, uflow
integer(kind=iwp) :: j, jtw, jtwm1
real(kind=wp), external :: d1mach

! the abscissae and weights are given for the interval (-1,1).
! because of symmetry only the positive abscissae and their
! corresponding weights are given.
!
! xgk    - abscissae of the 21-point kronrod rule
!          xgk(2), xgk(4), ...  abscissae of the 10-point
!          gauss rule
!          xgk(1), xgk(3), ...  abscissae which are optimally
!          added to the 10-point gauss rule
!
! wgk    - weights of the 21-point kronrod rule
!
! wg     - weights of the 10-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.

real(kind=wp), parameter :: wg(5) = [0.066671344308688137593568809893332_wp, &
                                     0.149451349150580593145776339657697_wp, &
                                     0.219086362515982043995534934228163_wp, &
                                     0.269266719309996355091226921569469_wp, &
                                     0.295524224714752870173892994651338_wp]
real(kind=wp), parameter :: xgk(11) = [0.995657163025808080735527280689003_wp, &
                                       0.973906528517171720077964012084452_wp, &
                                       0.930157491355708226001207180059508_wp, &
                                       0.865063366688984510732096688423493_wp, &
                                       0.780817726586416897063717578345042_wp, &
                                       0.679409568299024406234327365114874_wp, &
                                       0.562757134668604683339000099272694_wp, &
                                       0.433395394129247190799265943165784_wp, &
                                       0.294392862701460198131126603103866_wp, &
                                       0.148874338981631210884826001129720_wp, &
                                       0.000000000000000000000000000000000_wp]
real(kind=wp), parameter :: wgk(11) = [0.011694638867371874278064396062192_wp, &
                                       0.032558162307964727478818972459390_wp, &
                                       0.054755896574351996031381300244580_wp, &
                                       0.075039674810919952767043140916190_wp, &
                                       0.093125454583697605535065465083366_wp, &
                                       0.109387158802297641899210590325805_wp, &
                                       0.123491976262065851077958109831074_wp, &
                                       0.134709217311473325928054001771707_wp, &
                                       0.142775938577060080797094273138717_wp, &
                                       0.147739104901338491374841515972068_wp, &
                                       0.149445554002916905664936468389821_wp]

! list of major variables
! -----------------------
!
! centr  - mid point of the interval
! hlgth  - half-length of the interval
! absc   - abscissa
! fval*  - function value
! resg   - result of the 10-point gauss formula
! resk   - result of the 21-point kronrod formula
! reskh  - approximation to the mean value of f over (a,b),
!          i.e. to i/(b-a)
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! uflow is the smallest positive magnitude.
!
!***first executable statement  dqk21

epmach = d1mach(4)
uflow = d1mach(1)

centr = Half*(a+b)
hlgth = Half*(b-a)
dhlgth = abs(hlgth)

! compute the 21-point kronrod approximation to
! the integral, and estimate the absolute error.

fc = f(centr)
resg = Zero
resk = fc*wgk(11)
resabs = abs(resk)
do j=1,5
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
do j=1,5
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
resasc = wgk(11)*abs(fc-reskh)
do j=1,10
  resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
end do
reslt = resk*hlgth
resabs = resabs*dhlgth
resasc = resasc*dhlgth
abserr = abs((resk-resg)*hlgth)
if ((resasc /= Zero) .and. (abserr /= Zero)) abserr = resasc*min(One,(200.0_wp*abserr/resasc)**OneHalf)
if (resabs > uflow/(50.0_wp*epmach)) abserr = max((epmach*50.0_wp)*resabs,abserr)

return

end subroutine dqk21
