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

subroutine dqagi(f,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,lenw,last,iwork,work)
!***begin prologue  dqagi
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
!            or i = integral of f over (-infinity,+infinity)
!            hopefully satisfying following claim for accuracy
!            abs(i-result).le.max(epsabs,epsrel*abs(i)).
!***description
!
!        integration over infinite intervals
!        standard fortran subroutine
!
!        parameters
!         on entry
!            f      - real*8
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            bound  - real*8
!                     finite bound of integration range
!                     (has no meaning if interval is doubly-infinite)
!
!            inf    - integer
!                     indicating the kind of integration range involved
!                     inf = 1 corresponds to  (bound,+infinity),
!                     inf = -1            to  (-infinity,bound),
!                     inf = 2             to (-infinity,+infinity).
!
!            epsabs - real*8
!                     absolute accuracy requested
!            epsrel - real*8
!                     relative accuracy requested
!                     if  epsabs.le.0
!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!
!         on return
!            result - real*8
!                     approximation to the integral
!
!            abserr - real*8
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                   - ier.gt.0 abnormal termination of the routine. the
!                             estimates for result and error are less
!                             reliable. it is assumed that the requested
!                             accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (e.g. singularity,
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
!                             (epsabs.le.0 and
!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                              or limit.lt.1 or leniw.lt.limit*4.
!                             result, abserr, neval, last are set to
!                             zero. exept when limit or leniw is
!                             invalid, iwork(1), work(limit*2+1) and
!                             work(limit*3+1) are set to zero, work(1)
!                             is set to a and work(limit+1) to b.
!
!         dimensioning parameters
!            limit - integer
!                    dimensioning parameter for iwork
!                    limit determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (a,b), limit.ge.1.
!                    if limit.lt.1, the routine will end with ier = 6.
!
!            lenw  - integer
!                    dimensioning parameter for work
!                    lenw must be at least limit*4.
!                    if lenw.lt.limit*4, the routine will end
!                    with ier = 6.
!
!            last  - integer
!                    on return, last equals the number of subintervals
!                    produced in the subdivision process, which
!                    determines the number of significant elements
!                    actually in the work arrays.
!
!         work arrays
!            iwork - integer
!                    vector of dimension at least limit, the first
!                    k elements of which contain pointers
!                    to the error estimates over the subintervals,
!                    such that work(limit*3+iwork(1)),... ,
!                    work(limit*3+iwork(k)) form a decreasing
!                    sequence, with k = last if last.le.(limit/2+2), and
!                    k = limit+1-last otherwise
!
!            work  - real*8
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left
!                     end points of the subintervals in the
!                     partition of (a,b),
!                    work(limit+1), ..., work(limit+last) contain
!                     the right end points,
!                    work(limit*2+1), ...,work(limit*2+last) contain the
!                     integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3)
!                     contain the error estimates.
!***references  (none)
!***routines called  dqagie,xerror
!***end prologue  dqagi

real*8 abserr, bound, epsabs, epsrel, f, result, work
integer ier, inf, iwork, last, lenw, limit, lvl, l1, l2, l3, neval

dimension iwork(limit), work(lenw)

external f

!***first executable statement  dqagi

! check validity of limit and lenw.

ier = 6
neval = 0
last = 0
result = 0.0d+00
abserr = 0.0d+00
if (limit < 1 .or. lenw < limit*4) go to 10

! prepare call for dqagie.

l1 = limit+1
l2 = limit+l1
l3 = limit+l2

call dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)

! call error handler if necessary.

lvl = 0
10 continue
if (ier == 6) lvl = 1
if (ier /= 0) call xerror('abnormal return from dqagi',26,ier,lvl)

return

end subroutine dqagi

subroutine dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
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
!            abs(i-result).le.max(epsabs,epsrel*abs(i))
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
!                     if  epsabs.le.0
!                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!            limit  - integer
!                     gives an upper bound on the number of subintervals
!                     in the partition of (a,b), limit.ge.1
!
!         on return
!            result - real*8
!                     approximation to the integral
!
!            abserr - real*8
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                   - ier.gt.0 abnormal termination of the routine. the
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
!                             (epsabs.le.0 and
!                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                             result, abserr, neval, last, rlist(1),
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
!                     if last.le.(limit/2+2), and k = limit+1-last
!                     otherwise
!
!            last   - integer
!                     number of subintervals actually produced
!                     in the subdivision process
!
!***references  (none)
!***routines called  d1mach,dqelg,dqk15i,dqpsrt
!***end prologue  dqagie

real*8 abseps, abserr, alist, area, area1, area12, area2, a1, a2, blist, boun, bound, b1, b2, correc, defabs, defab1, defab2, &
       dres, d1mach, elist, epmach, epsabs, epsrel, erlarg, erlast, errbnd, errmax, error1, error2, erro12, errsum, ertest, f, &
       oflow, resabs, reseps, result, res3la, rlist, rlist2, small, uflow
integer id, ier, ierro, inf, iord, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, last, limit, maxerr, neval, nres, nrmax, numrl2
logical extrap, noext

dimension alist(limit), blist(limit), elist(limit), iord(limit), res3la(3), rlist(limit), rlist2(52)

external f

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
!             abs(result))
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
result = 0.0d+00
abserr = 0.0d+00
alist(1) = 0.0d+00
blist(1) = 0.1d+01
rlist(1) = 0.0d+00
elist(1) = 0.0d+00
iord(1) = 0

correc = 0.0d+00 ! dummy initialize
erlarg = 0.0d+00 ! dummy initialize
ertest = 0.0d+00 ! dummy initialize
small = 0.0d+00 ! dummy initialize

if (epsabs <= 0.0d+00 .and. epsrel < max(0.5d+02*epmach,0.5d-28)) ier = 6
if (ier == 6) go to 999

! first approximation to the integral
! -----------------------------------
!
! determine the interval to be mapped onto (0,1).
! if inf = 2 the integral is computed as i = i1+i2, where
! i1 = integral of f over (-infinity,0),
! i2 = integral of f over (0,+infinity).

boun = bound
if (inf == 2) boun = 0.0d+00
call dqk15i(f,boun,inf,0.0d+00,0.1d+01,result,abserr,defabs,resabs)

! test on accuracy

last = 1
rlist(1) = result
elist(1) = abserr
iord(1) = 1
dres = abs(result)
errbnd = max(epsabs,epsrel*dres)
if (abserr <= 1.0d+02*epmach*defabs .and. abserr > errbnd) ier = 2
if (limit == 1) ier = 1
if (ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs) .or. abserr == 0.0d+00) go to 130

! initialization
! --------------

uflow = d1mach(1)
oflow = d1mach(2)
rlist2(1) = result
errmax = abserr
maxerr = 1
area = result
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
if (dres >= (0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1

! main do-loop
! ------------

do 90 last=2,limit

  ! bisect the subinterval with nrmax-th largest error estimate.

  a1 = alist(maxerr)
  b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
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
  if (defab1 == error1 .or. defab2 == error2) go to 15
  if (abs(rlist(maxerr)-area12) > 0.1d-04*abs(area12) .or. erro12 < 0.99d+00*errmax) go to 10
  if (extrap) iroff2 = iroff2+1
  if (.not. extrap) iroff1 = iroff1+1
  10 continue
  if (last > 10 .and. erro12 > errmax) iroff3 = iroff3+1
  15 continue
  rlist(maxerr) = area1
  rlist(last) = area2
  errbnd = max(epsabs,epsrel*abs(area))

  ! test for roundoff error and eventually set error flag.

  if (iroff1+iroff2 >= 10 .or. iroff3 >= 20) ier = 2
  if (iroff2 >= 5) ierro = 3

  ! set error flag in the case that the number of
  ! subintervals equals limit.

  if (last == limit) ier = 1

  ! set error flag in the case of bad integrand behaviour
  ! at some points of the integration range.

  if (max(abs(a1),abs(b2)) <= (0.1d+01+0.1d+03*epmach)*(abs(a2)+0.1d+04*uflow)) ier = 4

  ! append the newly-created intervals to the list.

  if (error2 > error1) go to 20
  alist(last) = a2
  blist(maxerr) = b1
  blist(last) = b2
  elist(maxerr) = error1
  elist(last) = error2
  go to 30
  20 continue
  alist(maxerr) = a2
  alist(last) = a1
  blist(last) = b1
  rlist(maxerr) = area2
  rlist(last) = area1
  elist(maxerr) = error2
  elist(last) = error1

  ! call subroutine dqpsrt to maintain the descending ordering
  ! in the list of error estimates and select the subinterval
  ! with nrmax-th largest error estimate (to be bisected next).

  30 continue
  call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
  if (errsum <= errbnd) go to 115
  if (ier /= 0) go to 100
  if (last == 2) go to 80
  if (noext) go to 90
  erlarg = erlarg-erlast
  if (abs(b1-a1) > small) erlarg = erlarg+erro12
  if (extrap) go to 40

  ! test whether the interval to be bisected next is the
  ! smallest interval.

  if (abs(blist(maxerr)-alist(maxerr)) > small) go to 90
  extrap = .true.
  nrmax = 2
  40 continue
  if (ierro == 3 .or. erlarg <= ertest) go to 60

  ! the smallest interval has the largest error.
  ! before bisecting decrease the sum of the errors over the
  ! larger intervals (erlarg) and perform extrapolation.

  id = nrmax
  jupbnd = last
  if (last > (2+limit/2)) jupbnd = limit+3-last
  do 50 k=id,jupbnd
    maxerr = iord(nrmax)
    errmax = elist(maxerr)
    if (abs(blist(maxerr)-alist(maxerr)) > small) go to 90
    nrmax = nrmax+1
  50 continue

  ! perform extrapolation.

  60 continue
  numrl2 = numrl2+1
  rlist2(numrl2) = area
  call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
  ktmin = ktmin+1
  if (ktmin > 5 .and. abserr < 0.1d-02*errsum) ier = 5
  if (abseps >= abserr) go to 70
  ktmin = 0
  abserr = abseps
  result = reseps
  correc = erlarg
  ertest = max(epsabs,epsrel*abs(reseps))
  if (abserr <= ertest) go to 100

  ! prepare bisection of the smallest interval.

  70 continue
  if (numrl2 == 1) noext = .true.
  if (ier == 5) go to 100
  maxerr = iord(1)
  errmax = elist(maxerr)
  nrmax = 1
  extrap = .false.
  small = small*0.5d+00
  erlarg = errsum
  go to 90
  80 continue
   small = 0.375d+00
  erlarg = errsum
  ertest = errbnd
  rlist2(2) = area
90 continue

! set final result and error estimate.
! ------------------------------------

100 continue
if (abserr == oflow) go to 115
if ((ier+ierro) == 0) go to 110
if (ierro == 3) abserr = abserr+correc
if (ier == 0) ier = 3
if (result /= 0.0d+00 .and. area /= 0.0d+00) go to 105
if (abserr > errsum) go to 115
if (area == 0.0d+00) go to 130
go to 110
105 continue
if (abserr/abs(result) > errsum/abs(area)) go to 115

! test on divergence

110 continue
if (ksgn == (-1) .and. max(abs(result),abs(area)) <= defabs*0.1d-01) go to 130
if (0.1d-01 > (result/area) .or. (result/area) > 0.1d+03 .or. errsum > abs(area)) ier = 6
go to 130

! compute global integral sum.

115 continue
result = 0.0d+00
do 120 k=1,last
  result = result+rlist(k)
120 continue
abserr = errsum
130 continue
neval = 30*last-15
if (inf == 2) neval = 2*neval
if (ier > 2) ier = ier-1
999 continue

return

end subroutine dqagie

subroutine dqelg(n,epstab,result,abserr,res3la,nres)
!***begin prologue  dqelg
!***refer to  dqagie,dqagoe,dqagpe,dqagse
!***routines called  d1mach
!***revision date  830518   (yymmdd)
!***keywords  epsilon algorithm, convergence acceleration,
!             extrapolation
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!           de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine determines the limit of a given sequence of
!            approximations, by means of the epsilon algorithm of
!            p.wynn. an estimate of the absolute error is also given.
!            the condensed epsilon table is computed. only those
!            elements needed for the computation of the next diagonal
!            are preserved.
!***description
!
!           epsilon algorithm
!           standard fortran subroutine
!           real*8 version
!
!           parameters
!              n      - integer
!                       epstab(n) contains the new element in the
!                       first column of the epsilon table.
!
!              epstab - real*8
!                       vector of dimension 52 containing the elements
!                       of the two lower diagonals of the triangular
!                       epsilon table. the elements are numbered
!                       starting at the right-hand corner of the
!                       triangle.
!
!              result - real*8
!                       resulting approximation to the integral
!
!              abserr - real*8
!                       estimate of the absolute error computed from
!                       result and the 3 previous results
!
!              res3la - real*8
!                       vector of dimension 3 containing the last 3
!                       results
!
!              nres   - integer
!                       number of calls to the routine
!                       (should be zero at first call)
!
!***end prologue  dqelg

real*8 abserr, delta1, delta2, delta3, d1mach, epmach, epsinf, epstab, error, err1, err2, err3, e0, e1, e1abs, e2, e3, oflow, res, &
       result, res3la, ss, tol1, tol2, tol3
integer i, ib, ib2, ie, indx, k1, k2, k3, limexp, n, newelm, nres, num
dimension epstab(52), res3la(3)

! list of major variables
! -----------------------
!
! e0     - the 4 elements on which the computation of a new
! e1       element in the epsilon table is based
! e2
! e3                 e0
!              e3    e1    new
!                    e2
! newelm - number of elements to be computed in the new
!          diagonal
! error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
! result - the element in the new diagonal with least value
!          of error
!
! machine dependent constants
! ---------------------------
!
! epmach is the largest relative spacing.
! oflow is the largest positive magnitude.
! limexp is the maximum number of elements the epsilon
! table can contain. if this number is reached, the upper
! diagonal of the epsilon table is deleted.
!
!***first executable statement  dqelg

epmach = d1mach(4)
oflow = d1mach(2)
nres = nres+1
abserr = oflow
result = epstab(n)
if (n < 3) go to 100
limexp = 50
epstab(n+2) = epstab(n)
newelm = (n-1)/2
epstab(n) = oflow
num = n
k1 = n
do 40 i=1,newelm
  k2 = k1-1
  k3 = k1-2
  res = epstab(k1+2)
  e0 = epstab(k3)
  e1 = epstab(k2)
  e2 = res
  e1abs = abs(e1)
  delta2 = e2-e1
  err2 = abs(delta2)
  tol2 = max(abs(e2),e1abs)*epmach
  delta3 = e1-e0
  err3 = abs(delta3)
  tol3 = max(e1abs,abs(e0))*epmach
  if (err2 > tol2 .or. err3 > tol3) go to 10

  ! if e0, e1 and e2 are equal to within machine
  ! accuracy, convergence is assumed.
  ! result = e2
  ! abserr = abs(e1-e0)+abs(e2-e1)

  result = res
  abserr = err2+err3
  ! ***jump out of do-loop
  go to 100
  10 continue
  e3 = epstab(k1)
  epstab(k1) = e1
  delta1 = e1-e3
  err1 = abs(delta1)
  tol1 = max(e1abs,abs(e3))*epmach

  ! if two elements are very close to each other, omit
  ! a part of the table by adjusting the value of n

  if (err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3) go to 20
  ss = 0.1d+01/delta1+0.1d+01/delta2-0.1d+01/delta3
  epsinf = abs(ss*e1)

  ! test to detect irregular behaviour in the table, and
  ! eventually omit a part of the table adjusting the value
  ! of n.

  if (epsinf > 0.1d-03) go to 30
  20 continue
  n = i+i-1
  ! ***jump out of do-loop
  go to 50

  ! compute a new element and eventually adjust
  ! the value of result.

  30 continue
  res = e1+0.1d+01/ss
  epstab(k1) = res
  k1 = k1-2
  error = err2+abs(res-e2)+err3
  if (error > abserr) go to 40
  abserr = error
  result = res
40 continue

! shift the table.

50 continue
if (n == limexp) n = 2*(limexp/2)-1
ib = 1
if ((num/2)*2 == num) ib = 2
ie = newelm+1
do 60 i=1,ie
  ib2 = ib+2
  epstab(ib) = epstab(ib2)
  ib = ib2
60 continue
if (num == n) go to 80
indx = num-n+1
do 70 i=1,n
  epstab(i) = epstab(indx)
  indx = indx+1
70 continue
80 continue
if (nres >= 4) go to 90
res3la(nres) = result
abserr = oflow
go to 100

! compute error estimate

90 continue
abserr = abs(result-res3la(3))+abs(result-res3la(2))+abs(result-res3la(1))
res3la(1) = res3la(2)
res3la(2) = res3la(3)
res3la(3) = result
100 continue
abserr = max(abserr,0.5d+01*epmach*abs(result))

return

end subroutine dqelg

subroutine dqk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
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
!                       fuction subprogram defining the integrand
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
!              result - real*8
!                       approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule(resk) obtained by optimal addition
!                       of abscissae to the 7-point gauss rule(resg).
!
!              abserr - real*8
!                       estimate of the modulus of the absolute error,
!                       which should equal or exceed abs(i-result)
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

real*8 a, absc, absc1, absc2, abserr, b, boun, centr, dinf, d1mach, epmach, f, fc, fsum, fval1, fval2, fv1, fv2, hlgth, resabs, &
       resasc, resg, resk, reskh, result, tabsc1, tabsc2, uflow, wg, wgk, xgk
integer inf, j
external f

dimension fv1(7), fv2(7), xgk(8), wgk(8), wg(8)

! the abscissae and weights are supplied for the interval
! (-1,1).  because of symmetry only the positive abscissae and
! their corresponding weights are given.
!
! xgk    - abscissae of the 15-point kronrod rule
!          xgk(2), xgk(4), ... abscissae of the 7-point
!          gauss rule
!          xgk(1), xgk(3), ...  abscissae which are optimally
!          added to the 7-point gauss rule
!
! wgk    - weights of the 15-point kronrod rule
!
! wg     - weights of the 7-point gauss rule, corresponding
!          to the abscissae xgk(2), xgk(4), ...
!          wg(1), wg(3), ... are set to zero.

data wg(1)/0.0d0/
data wg(2)/0.129484966168869693270611432679082d0/
data wg(3)/0.0d0/
data wg(4)/0.279705391489276667901467771423780d0/
data wg(5)/0.0d0/
data wg(6)/0.381830050505118944950369775488975d0/
data wg(7)/0.0d0/
data wg(8)/0.417959183673469387755102040816327d0/

data xgk(1)/0.991455371120812639206854697526329d0/
data xgk(2)/0.949107912342758524526189684047851d0/
data xgk(3)/0.864864423359769072789712788640926d0/
data xgk(4)/0.741531185599394439863864773280788d0/
data xgk(5)/0.586087235467691130294144838258730d0/
data xgk(6)/0.405845151377397166906606412076961d0/
data xgk(7)/0.207784955007898467600689403773245d0/
data xgk(8)/0.000000000000000000000000000000000d0/

data wgk(1)/0.022935322010529224963732008058970d0/
data wgk(2)/0.063092092629978553290700663189204d0/
data wgk(3)/0.104790010322250183839876322541518d0/
data wgk(4)/0.140653259715525918745189590510238d0/
data wgk(5)/0.169004726639267902826583426598550d0/
data wgk(6)/0.190350578064785409913256402421014d0/
data wgk(7)/0.204432940075298892414161999234649d0/
data wgk(8)/0.209482141084727828012999174891714d0/

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

centr = 0.5d+00*(a+b)
hlgth = 0.5d+00*(b-a)
tabsc1 = boun+dinf*(0.1d+01-centr)/centr
fval1 = f(tabsc1)
if (inf == 2) fval1 = fval1+f(-tabsc1)
fc = (fval1/centr)/centr

! compute the 15-point kronrod approximation to
! the integral, and estimate the error.

resg = wg(8)*fc
resk = wgk(8)*fc
resabs = abs(resk)
do 10 j=1,7
  absc = hlgth*xgk(j)
  absc1 = centr-absc
  absc2 = centr+absc
  tabsc1 = boun+dinf*(0.1d+01-absc1)/absc1
  tabsc2 = boun+dinf*(0.1d+01-absc2)/absc2
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
10 continue
reskh = resk*0.5d+00
resasc = wgk(8)*abs(fc-reskh)
do 20 j=1,7
  resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
20 continue
result = resk*hlgth
resasc = resasc*hlgth
resabs = resabs*hlgth
abserr = abs((resk-resg)*hlgth)
if (resasc /= 0.0d+00 .and. abserr /= 0.d0) abserr = resasc*min(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
if (resabs > uflow/(0.5d+02*epmach)) abserr = max((epmach*0.5d+02)*resabs,abserr)

return

end subroutine dqk15i
