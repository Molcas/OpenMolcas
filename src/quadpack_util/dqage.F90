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

subroutine dqage(f,a,b,epsabs,epsrel,key,limit,reslt,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
!***begin prologue  dqage
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
!***keywords  automatic integrator, general-purpose,
!             integrand examinator, globally adaptive,
!             gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            definite integral   i = integral of f over (a,b),
!            hopefully satisfying following claim for accuracy
!            abs(i-reslt) <= max(epsabs,epsrel*abs(i)).
!***description
!
!        computation of a definite integral
!        standard fortran subroutine
!        real*8 version
!
!        parameters
!         on entry
!            f      - real*8
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            a      - real*8
!                     lower limit of integration
!
!            b      - real*8
!                     upper limit of integration
!
!            epsabs - real*8
!                     absolute accuracy requested
!            epsrel - real*8
!                     relative accuracy requested
!                     if  epsabs <= 0
!                     and epsrel < max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            key    - integer
!                     key for choice of local integration rule
!                     a gauss-kronrod pair is used with
!                          7 - 15 points if key < 2,
!                         10 - 21 points if key = 2,
!                         15 - 31 points if key = 3,
!                         20 - 41 points if key = 4,
!                         25 - 51 points if key = 5,
!                         30 - 61 points if key > 5.
!
!            limit  - integer
!                     gives an upperbound on the number of subintervals
!                     in the partition of (a,b), limit >= 1.
!
!         on return
!            reslt  - real*8
!                     approximation to the integral
!
!            abserr - real*8
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-reslt)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for result and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value
!                             of limit.
!                             however, if this yields no improvement it
!                             is rather advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a local
!                             difficulty can be determined(e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. if possible, an appropriate
!                             special-purpose integrator should be used
!                             which is designed for handling the type of
!                             difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 the input is invalid, because
!                             (epsabs <= 0 and
!                              epsrel < max(50*rel.mach.acc.,0.5d-28),
!                             reslt, abserr, neval, last, rlist(1) ,
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to a and b
!                             respectively.
!
!            alist   - real*8
!                      vector of dimension at least limit, the first
!                       last  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            blist   - real*8
!                      vector of dimension at least limit, the first
!                       last  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (a,b)
!
!            rlist   - real*8
!                      vector of dimension at least limit, the first
!                       last  elements of which are the
!                      integral approximations on the subintervals
!
!            elist   - real*8
!                      vector of dimension at least limit, the first
!                       last  elements of which are the moduli of the
!                      absolute error estimates on the subintervals
!
!            iord    - integer
!                      vector of dimension at least limit, the first k
!                      elements of which are pointers to the
!                      error estimates over the subintervals,
!                      such that elist(iord(1)), ...,
!                      elist(iord(k)) form a decreasing sequence,
!                      with k = last if last <= (limit/2+2), and
!                      k = limit+1-last otherwise
!
!            last    - integer
!                      number of subintervals actually produced in the
!                      subdivision process
!
!***references  (none)
!***routines called  d1mach,dqk15,dqk21,dqk31,
!                    dqk41,dqk51,dqk61,dqpsrt
!***end prologue  dqage

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
interface
  function f(x)
    import :: wp
    real(kind=wp) :: f
    real(kind=wp), intent(in) :: x
  end function f
end interface
real(kind=wp), intent(in) :: a, b, epsabs, epsrel
integer(kind=iwp), intent(in) :: key, limit
real(kind=wp), intent(out) :: reslt, abserr, alist(limit), blist(limit), rlist(limit), elist(limit)
integer(kind=iwp), intent(out) :: neval, ier, iord(limit), last
real(kind=wp) :: area, area1, area12, area2, a1, a2, b1, b2, defabs, defab1, defab2, epmach, errbnd, errmax, error1, error2, &
                 erro12, errsum, resabs, uflow
integer(kind=iwp) :: iroff1, iroff2, k, keyf, maxerr, nrmax
real(kind=wp), external :: d1mach

!  list of major variables
!  -----------------------
!
! alist     - list of left end points of all subintervals
!             considered up to now
! blist     - list of right end points of all subintervals
!             considered up to now
! rlist(i)  - approximation to the integral over
!            (alist(i),blist(i))
! elist(i)  - error estimate applying to rlist(i)
! maxerr    - pointer to the interval with largest
!             error estimate
! errmax    - elist(maxerr)
! area      - sum of the integrals over the subintervals
! errsum    - sum of the errors over the subintervals
! errbnd    - requested accuracy max(epsabs,epsrel*
!             abs(reslt))
! *****1    - variable for the left subinterval
! *****2    - variable for the right subinterval
! last      - index for subdivision
!
!
! machine dependent constants
! ---------------------------
!
! epmach  is the largest relative spacing.
! uflow  is the smallest positive magnitude.
!
!***first executable statement  dqage

epmach = d1mach(4)
uflow = d1mach(1)

! test on validity of parameters
! ------------------------------

ier = 0
neval = 0
last = 0
reslt = Zero
abserr = Zero
alist(1) = a
blist(1) = b
rlist(1) = Zero
elist(1) = Zero
iord(1) = 0
if ((epsabs <= Zero) .and. (epsrel < max(50.0_wp*epmach,5.0e-29_wp))) ier = 6
if (ier == 6) return

! first approximation to the integral
! -----------------------------------

keyf = key
if (key <= 0) keyf = 1
if (key >= 7) keyf = 6
neval = 0
if (keyf == 1) call dqk15(f,a,b,reslt,abserr,defabs,resabs)
if (keyf == 2) call dqk21(f,a,b,reslt,abserr,defabs,resabs)
if (keyf == 3) call dqk31(f,a,b,reslt,abserr,defabs,resabs)
if (keyf == 4) call dqk41(f,a,b,reslt,abserr,defabs,resabs)
if (keyf == 5) call dqk51(f,a,b,reslt,abserr,defabs,resabs)
if (keyf == 6) call dqk61(f,a,b,reslt,abserr,defabs,resabs)
last = 1
rlist(1) = reslt
elist(1) = abserr
iord(1) = 1

! test on accuracy.

errbnd = max(epsabs,epsrel*abs(reslt))
if ((abserr <= 50.0_wp*epmach*defabs) .and. (abserr > errbnd)) ier = 2
if (limit == 1) ier = 1
if ((ier /= 0) .or. ((abserr <= errbnd) .and. (abserr /= resabs)) .or. (abserr == Zero)) then
  call finish()
  return
end if

! initialization
! --------------

errmax = abserr
maxerr = 1
area = reslt
errsum = abserr
nrmax = 1
iroff1 = 0
iroff2 = 0

! main do-loop
! ------------

do last=2,limit

  ! bisect the subinterval with the largest error estimate.

  a1 = alist(maxerr)
  b1 = Half*(alist(maxerr)+blist(maxerr))
  a2 = b1
  b2 = blist(maxerr)
  if (keyf == 1) call dqk15(f,a1,b1,area1,error1,resabs,defab1)
  if (keyf == 2) call dqk21(f,a1,b1,area1,error1,resabs,defab1)
  if (keyf == 3) call dqk31(f,a1,b1,area1,error1,resabs,defab1)
  if (keyf == 4) call dqk41(f,a1,b1,area1,error1,resabs,defab1)
  if (keyf == 5) call dqk51(f,a1,b1,area1,error1,resabs,defab1)
  if (keyf == 6) call dqk61(f,a1,b1,area1,error1,resabs,defab1)
  if (keyf == 1) call dqk15(f,a2,b2,area2,error2,resabs,defab2)
  if (keyf == 2) call dqk21(f,a2,b2,area2,error2,resabs,defab2)
  if (keyf == 3) call dqk31(f,a2,b2,area2,error2,resabs,defab2)
  if (keyf == 4) call dqk41(f,a2,b2,area2,error2,resabs,defab2)
  if (keyf == 5) call dqk51(f,a2,b2,area2,error2,resabs,defab2)
  if (keyf == 6) call dqk61(f,a2,b2,area2,error2,resabs,defab2)

  ! improve previous approximations to integral
  ! and error and test for accuracy.

  neval = neval+1
  area12 = area1+area2
  erro12 = error1+error2
  errsum = errsum+erro12-errmax
  area = area+area12-rlist(maxerr)
  if ((defab1 /= error1) .and. (defab2 /= error2)) then
    if ((abs(rlist(maxerr)-area12) <= 1.0e-5_wp*abs(area12)) .and. (erro12 >= 0.99_wp*errmax)) iroff1 = iroff1+1
    if ((last > 10) .and. (erro12 > errmax)) iroff2 = iroff2+1
  end if
  rlist(maxerr) = area1
  rlist(last) = area2
  errbnd = max(epsabs,epsrel*abs(area))
  if (errsum > errbnd) then

    ! test for roundoff error and eventually set error flag.

    if ((iroff1 >= 6) .or. (iroff2 >= 20)) ier = 2

    ! set error flag in the case that the number of subintervals
    ! equals limit.

    if (last == limit) ier = 1

    ! set error flag in the case of bad integrand behaviour
    ! at a point of the integration range.

    if (max(abs(a1),abs(b2)) <= (One+100.0_wp*epmach)*(abs(a2)+1000.0_wp*uflow)) ier = 3
  end if

  ! append the newly-created intervals to the list.

  if (error2 <= error1) then
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
  else
    alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
  end if

  ! call subroutine dqpsrt to maintain the descending ordering
  ! in the list of error estimates and select the subinterval
  ! with the largest error estimate (to be bisected next).

  call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
  ! ***jump out of do-loop
  if ((ier /= 0) .or. (errsum <= errbnd)) exit
end do

! compute final result.
! ---------------------

reslt = Zero
do k=1,last
  reslt = reslt+rlist(k)
end do
abserr = errsum

call finish()

return

contains

subroutine finish()
  if (keyf /= 1) neval = (10*keyf+1)*(2*neval+1)
  if (keyf == 1) neval = 30*neval+15
end subroutine finish

end subroutine dqage
