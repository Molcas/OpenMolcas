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

subroutine dqagie(f,bound,inf,epsabs,epsrel,limit,reslt,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
!***begin prologue  dqagie
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a3a1,h2a4a1
!***keywords  automatic integrator, infinite intervals,
!             general-purpose, transformation, extrapolation,
!             globally adaptive
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!            integral   i = integral of f over (bound,+infinity)
!            or i = integral of f over (-infinity,bound)
!            or i = integral of f over (-infinity,+infinity),
!            hopefully satisfying following claim for accuracy
!            abs(i-reslt) <= max(epsabs,epsrel*abs(i))
!***description
!
! integration over infinite intervals
! standard fortran subroutine
!
!            f      - real*8
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            bound  - real*8
!                     finite bound of integration range
!                     (has no meaning if interval is doubly-infinite)
!
!            inf    - real*8
!                     indicating the kind of integration range involved
!                     inf = 1 corresponds to  (bound,+infinity),
!                     inf = -1            to  (-infinity,bound),
!                     inf = 2             to (-infinity,+infinity).
!
!            epsabs - real*8
!                     absolute accuracy requested
!            epsrel - real*8
!                     relative accuracy requested
!                     if  epsabs <= 0
!                     and epsrel < max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit >= 1
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
!                   - ier > 0 abnormal termination of the routine. the
!                             estimates for result and error are less
!                             reliable. it is assumed that the requested
!                             accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however,if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties.
!                             if the position of a local difficulty can
!                             be determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is assumed that the requested tolerance
!                             cannot be achieved, and that the returned
!                             result is the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             (epsabs <= 0 and
!                              epsrel < max(50*rel.mach.acc.,0.5d-28),
!                             reslt, abserr, neval, last, rlist(1),
!                             elist(1) and iord(1) are set to zero.
!                             alist(1) and blist(1) are set to 0
!                             and 1 respectively.
!
!            alist  - real*8
!                     vector of dimension at least limit, the first
!                      last  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the transformed integration range (0,1).
!
!            blist  - real*8
!                     vector of dimension at least limit, the first
!                      last  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the transformed integration range (0,1).
!
!            rlist  - real*8
!                     vector of dimension at least limit, the first
!                      last  elements of which are the integral
!                     approximations on the subintervals
!
!            elist  - real*8
!                     vector of dimension at least limit,  the first
!                     last elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            iord   - integer
!                     vector of dimension limit, the first k
!                     elements of which are pointers to the
!                     error estimates over the subintervals,
!                     such that elist(iord(1)), ..., elist(iord(k))
!                     form a decreasing sequence, with k = last
!                     if last <= (limit/2+2), and k = limit+1-last
!                     otherwise
!
!            last   - integer
!                     number of subintervals actually produced
!                     in the subdivision process
!
!***references  (none)
!***routines called  d1mach,dqelg,dqk15i,dqpsrt
!***end prologue  dqagie

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
real(kind=wp), intent(in) :: bound, epsabs, epsrel
integer(kind=iwp), intent(in) :: inf, limit
real(kind=wp), intent(out) :: reslt, abserr, alist(limit), blist(limit), rlist(limit), elist(limit)
integer(kind=iwp), intent(out) :: neval, ier, iord(limit), last
real(kind=wp) :: abseps, area, area1, area12, area2, a1, a2, boun, b1, b2, correc, defabs, defab1, defab2, dres, epmach, erlarg, &
                 erlast, errbnd, errmax, error1, error2, erro12, errsum, ertest, oflow, resabs, reseps, res3la(3), rlist2(52), &
                 small, uflow
integer(kind=iwp) :: id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, maxerr, nres, nrmax, numrl2
logical(kind=iwp) :: extrap, noext
real(kind=wp), external :: d1mach

!  the dimension of rlist2 is determined by the value of
!  limexp in subroutine dqelg.
!
!
!  list of major variables
!  -----------------------
!
! alist     - list of left end points of all subintervals
!             considered up to now
! blist     - list of right end points of all subintervals
!             considered up to now
! rlist(i)  - approximation to the integral over
!             (alist(i),blist(i))
! rlist2    - array of dimension at least (limexp+2),
!             containing the part of the epsilon table
!             wich is still needed for further computations
! elist(i)  - error estimate applying to rlist(i)
! maxerr    - pointer to the interval with largest error
!             estimate
! errmax    - elist(maxerr)
! erlast    - error on the interval currently subdivided
!             (before that subdivision has taken place)
! area      - sum of the integrals over the subintervals
! errsum    - sum of the errors over the subintervals
! errbnd    - requested accuracy max(epsabs,epsrel*
!             abs(reslt))
! *****1    - variable for the left subinterval
! *****2    - variable for the right subinterval
! last      - index for subdivision
! nres      - number of calls to the extrapolation routine
! numrl2    - number of elements currently in rlist2. if an
!             appropriate approximation to the compounded
!             integral has been obtained, it is put in
!             rlist2(numrl2) after numrl2 has been increased
!             by one.
! small     - length of the smallest interval considered up
!             to now, multiplied by 1.5
! erlarg    - sum of the errors over the intervals larger
!             than the smallest interval considered up to now
! extrap    - logical variable denoting that the routine
!             is attempting to perform extrapolation. i.e.
!             before subdividing the smallest interval we
!             try to decrease the value of erlarg.
! noext     - logical variable denoting that extrapolation
!             is no longer allowed (true-value)
!
!  machine dependent constants
!  ---------------------------
!
! epmach is the largest relative spacing.
! uflow is the smallest positive magnitude.
! oflow is the largest positive magnitude.

!***first executable statement  dqagie

epmach = d1mach(4)

! test on validity of parameters
! -----------------------------

ier = 0
neval = 0
last = 0
reslt = Zero
abserr = Zero
alist(1) = Zero
blist(1) = One
rlist(1) = Zero
elist(1) = Zero
iord(1) = 0

correc = Zero ! dummy initialize
erlarg = Zero ! dummy initialize
ertest = Zero ! dummy initialize
small = Zero ! dummy initialize

if ((epsabs <= Zero) .and. (epsrel < max(50.0_wp*epmach,0.5e-28_wp))) ier = 6
if (ier == 6) return

! first approximation to the integral
! -----------------------------------
!
! determine the interval to be mapped onto (0,1).
! if inf = 2 the integral is computed as i = i1+i2, where
! i1 = integral of f over (-infinity,0),
! i2 = integral of f over (0,+infinity).

boun = bound
if (inf == 2) boun = Zero
call dqk15i(f,boun,inf,Zero,One,reslt,abserr,defabs,resabs)

! test on accuracy

last = 1
rlist(1) = reslt
elist(1) = abserr
iord(1) = 1
dres = abs(reslt)
errbnd = max(epsabs,epsrel*dres)
if ((abserr <= 100.0_wp*epmach*defabs) .and. (abserr > errbnd)) ier = 2
if (limit == 1) ier = 1
if ((ier /= 0) .or. ((abserr <= errbnd) .and. (abserr /= resabs)) .or. (abserr == Zero)) then
  call finish(compute=.false.)
  return
end if

! initialization
! --------------

uflow = d1mach(1)
oflow = d1mach(2)
rlist2(1) = reslt
errmax = abserr
maxerr = 1
area = reslt
errsum = abserr
abserr = oflow
nrmax = 1
nres = 0
ktmin = 0
numrl2 = 2
extrap = .false.
noext = .false.
ierro = 0
iroff1 = 0
iroff2 = 0
iroff3 = 0
ksgn = -1
if (dres >= (One-50.0_wp*epmach)*defabs) ksgn = 1

! main do-loop
! ------------

main: do last=2,limit

  ! bisect the subinterval with nrmax-th largest error estimate.

  a1 = alist(maxerr)
  b1 = Half*(alist(maxerr)+blist(maxerr))
  a2 = b1
  b2 = blist(maxerr)
  erlast = errmax
  call dqk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1)
  call dqk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2)

  ! improve previous approximations to integral
  ! and error and test for accuracy.

  area12 = area1+area2
  erro12 = error1+error2
  errsum = errsum+erro12-errmax
  area = area+area12-rlist(maxerr)
  if ((defab1 /= error1) .and. (defab2 /= error2)) then
    if ((abs(rlist(maxerr)-area12) <= 0.1e-4_wp*abs(area12)) .and. (erro12 >= 0.99_wp*errmax)) then
      if (extrap) iroff2 = iroff2+1
      if (.not. extrap) iroff1 = iroff1+1
    end if
    if ((last > 10) .and. (erro12 > errmax)) iroff3 = iroff3+1
  end if
  rlist(maxerr) = area1
  rlist(last) = area2
  errbnd = max(epsabs,epsrel*abs(area))

  ! test for roundoff error and eventually set error flag.

  if ((iroff1+iroff2 >= 10) .or. (iroff3 >= 20)) ier = 2
  if (iroff2 >= 5) ierro = 3

  ! set error flag in the case that the number of
  ! subintervals equals limit.

  if (last == limit) ier = 1

  ! set error flag in the case of bad integrand behaviour
  ! at some points of the integration range.

  if (max(abs(a1),abs(b2)) <= (One+100.0_wp*epmach)*(abs(a2)+1000.0_wp*uflow)) ier = 4

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
  ! with nrmax-th largest error estimate (to be bisected next).

  call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
  if (errsum <= errbnd) then
    call finish(compute=.true.)
    return
  end if
  if (ier /= 0) exit main
  if (last == 2) then
    small = 0.375_wp
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area
    cycle main
  end if
  if (noext) cycle main
  erlarg = erlarg-erlast
  if (abs(b1-a1) > small) erlarg = erlarg+erro12
  if (.not. extrap) then

    ! test whether the interval to be bisected next is the
    ! smallest interval.

    if (abs(blist(maxerr)-alist(maxerr)) > small) cycle main
    extrap = .true.
    nrmax = 2
  end if
  if ((ierro /= 3) .and. (erlarg > ertest)) then

    ! the smallest interval has the largest error.
    ! before bisecting decrease the sum of the errors over the
    ! larger intervals (erlarg) and perform extrapolation.

    id = nrmax
    jupbnd = last
    if (last > (2+limit/2)) jupbnd = limit+3-last
    do k=id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if (abs(blist(maxerr)-alist(maxerr)) > small) cycle main
      nrmax = nrmax+1
    end do
  end if

  ! perform extrapolation.

  numrl2 = numrl2+1
  rlist2(numrl2) = area
  call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
  ktmin = ktmin+1
  if ((ktmin > 5) .and. (abserr < 1.0e-3_wp*errsum)) ier = 5
  if (abseps < abserr) then
    ktmin = 0
    abserr = abseps
    reslt = reseps
    correc = erlarg
    ertest = max(epsabs,epsrel*abs(reseps))
    if (abserr <= ertest) exit main
  end if

  ! prepare bisection of the smallest interval.

  if (numrl2 == 1) noext = .true.
  if (ier == 5) exit main
  maxerr = iord(1)
  errmax = elist(maxerr)
  nrmax = 1
  extrap = .false.
  small = small*Half
  erlarg = errsum
end do main

! set final result and error estimate.
! ------------------------------------

if (abserr == oflow) then
  call finish(compute=.true.)
  return
end if
if ((ier+ierro) /= 0) then
  if (ierro == 3) abserr = abserr+correc
  if (ier == 0) ier = 3
  if ((reslt == Zero) .or. (area == Zero)) then
    if (abserr > errsum) then
      call finish(compute=.true.)
      return
    end if
    if (area == Zero) then
      call finish(compute=.false.)
      return
    end if
  else
    if (abserr/abs(reslt) > errsum/abs(area)) then
      call finish(compute=.true.)
      return
    end if
  end if
end if

! test on divergence

if ((ksgn /= -1) .or. (max(abs(reslt),abs(area)) > defabs*1.0e-2_wp)) then
  if ((1.0e-2_wp > reslt/area) .or. (reslt/area > 100.0_wp) .or. (errsum > abs(area))) ier = 6
end if
call finish(compute=.false.)

return

contains

subroutine finish(compute)

  logical(kind=iwp), intent(in) :: compute
  integer(kind=iwp) :: i

  if (compute) then

    ! compute global integral sum.

    reslt = Zero
    do i=1,last
      reslt = reslt+rlist(i)
    end do
    abserr = errsum

  end if

  neval = 30*last-15
  if (inf == 2) neval = 2*neval
  if (ier > 2) ier = ier-1

end subroutine finish

end subroutine dqagie
