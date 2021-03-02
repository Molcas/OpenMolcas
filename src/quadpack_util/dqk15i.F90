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

subroutine dqk15i(f,boun,inf,a,b,reslt,abserr,resabs,resasc)
!***begin prologue  dqk15i
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a2,h2a4a2
!***keywords  15-point transformed gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the original (infinite integration range is mapped
!            onto the interval (0,1) and (a,b) is a part of (0,1).
!            it is the purpose to compute
!            i = integral of transformed integrand over (a,b),
!            j = integral of abs(transformed integrand) over (a,b).
!***description
!
!           integration rule
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
!              boun   - real*8
!                       finite bound of original integration
!                       range (set to zero if inf = +2)
!
!              inf    - integer
!                       if inf = -1, the original interval is
!                                   (-infinity,bound),
!                       if inf = +1, the original interval is
!                                   (bound,+infinity),
!                       if inf = +2, the original interval is
!                                   (-infinity,+infinity) and
!                       the integral is computed as the sum of two
!                       integrals, one over (-infinity,0) and one over
!                       (0,+infinity).
!
!              a      - real*8
!                       lower limit for integration over subrange
!                       of (0,1)
!
!              b      - real*8
!                       upper limit for integration over subrange
!                       of (0,1)
!
!            on return
!              reslt  - real*8
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the 7-point gauss rule (resg).
!
!              abserr - real*8
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-reslt)
!
!              resabs - real*8
!                       approximation to the integral j
!
!              resasc - real*8
!                       approximation to the integral of
!                       abs((transformed integrand)-i/(b-a)) over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  dqk15i

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
real(kind=wp), intent(in) :: boun, a, b
real(kind=wp), intent(out) :: reslt, abserr, resabs, resasc
integer(kind=iwp), intent(in) :: inf
integer(kind=iwp) :: j
real(kind=wp) :: absc, absc1, absc2, centr, dinf, epmach, fc, fsum, fval1, fval2, fv1(7), fv2(7), hlgth, resg, resk, reskh, &
                 tabsc1, tabsc2, uflow
real(kind=wp), external :: d1mach

! the abscissae and weights are supplied for the interval
! (-1,1).  because of symmetry only the positive abscissae and
! their corresponding weights are given.
!
! xgk    - abscissae of the 15-point kronrod rule
!          xgk(2), xgk(4), ...  abscissae of the 7-point
!          gauss rule
!          xgk(1), xgk(3), ...  abscissae which are optimally
!          added to the 7-point gauss rule
!
! wgk    - weights of the 15-point kronrod rule
!
! wg     - weights of the 7-point gauss rule, corresponding
!          to the abscissae xgk(2), xgk(4), ...
!          wg(1), wg(3), ... are set to zero.

real(kind=wp), parameter :: wg(8) = [0.000000000000000000000000000000000_wp, &
                                     0.129484966168869693270611432679082_wp, &
                                     0.000000000000000000000000000000000_wp, &
                                     0.279705391489276667901467771423780_wp, &
                                     0.000000000000000000000000000000000_wp, &
                                     0.381830050505118944950369775488975_wp, &
                                     0.000000000000000000000000000000000_wp, &
                                     0.417959183673469387755102040816327_wp]
real(kind=wp), parameter :: xgk(8) = [0.991455371120812639206854697526329_wp, &
                                      0.949107912342758524526189684047851_wp, &
                                      0.864864423359769072789712788640926_wp, &
                                      0.741531185599394439863864773280788_wp, &
                                      0.586087235467691130294144838258730_wp, &
                                      0.405845151377397166906606412076961_wp, &
                                      0.207784955007898467600689403773245_wp, &
                                      0.000000000000000000000000000000000_wp]
real(kind=wp), parameter :: wgk(8) = [0.022935322010529224963732008058970_wp, &
                                      0.063092092629978553290700663189204_wp, &
                                      0.104790010322250183839876322541518_wp, &
                                      0.140653259715525918745189590510238_wp, &
                                      0.169004726639267902826583426598550_wp, &
                                      0.190350578064785409913256402421014_wp, &
                                      0.204432940075298892414161999234649_wp, &
                                      0.209482141084727828012999174891714_wp]

! list of major variables
! -----------------------
!
! centr  - mid point of the interval
! hlgth  - half-length of the interval
! absc*  - abscissa
! tabsc* - transformed abscissa
! fval*  - function value
! resg   - result of the 7-point gauss formula
! resk   - result of the 15-point kronrod formula
! reskh  - approximation to the mean value of the transformed
!          integrand over (a,b), i.e. to i/(b-a)
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! uflow is the smallest positive magnitude.
!
!***first executable statement  dqk15i

epmach = d1mach(4)
uflow = d1mach(1)
dinf = min(1,inf)

centr = Half*(a+b)
hlgth = Half*(b-a)
tabsc1 = boun+dinf*(One-centr)/centr
fval1 = f(tabsc1)
if (inf == 2) fval1 = fval1+f(-tabsc1)
fc = (fval1/centr)/centr

! compute the 15-point kronrod approximation to
! the integral, and estimate the error.

resg = wg(8)*fc
resk = wgk(8)*fc
resabs = abs(resk)
do j=1,7
  absc = hlgth*xgk(j)
  absc1 = centr-absc
  absc2 = centr+absc
  tabsc1 = boun+dinf*(One-absc1)/absc1
  tabsc2 = boun+dinf*(One-absc2)/absc2
  fval1 = f(tabsc1)
  fval2 = f(tabsc2)
  if (inf == 2) fval1 = fval1+f(-tabsc1)
  if (inf == 2) fval2 = fval2+f(-tabsc2)
  fval1 = (fval1/absc1)/absc1
  fval2 = (fval2/absc2)/absc2
  fv1(j) = fval1
  fv2(j) = fval2
  fsum = fval1+fval2
  resg = resg+wg(j)*fsum
  resk = resk+wgk(j)*fsum
  resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
end do
reskh = resk*Half
resasc = wgk(8)*abs(fc-reskh)
do j=1,7
  resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
end do
reslt = resk*hlgth
resabs = resabs*hlgth
resasc = resasc*hlgth
abserr = abs((resk-resg)*hlgth)
if ((resasc /= Zero) .and. (abserr /= Zero)) abserr = resasc*min(One,(200.0_wp*abserr/resasc)**OneHalf)
if (resabs > uflow/(50.0_wp*epmach)) abserr = max((epmach*50.0_wp)*resabs,abserr)

return

end subroutine dqk15i
