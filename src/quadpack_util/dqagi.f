************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine dqagi(f,bound,inf,epsabs,epsrel,result,abserr,neval,
     *   ier,limit,lenw,last,iwork,work)
c***begin prologue  dqagi
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1,h2a4a1
c***keywords  automatic integrator, infinite intervals,
c             general-purpose, transformation, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. -k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            integral   i = integral of f over (bound,+infinity)
c            or i = integral of f over (-infinity,bound)
c            or i = integral of f over (-infinity,+infinity)
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c        integration over infinite intervals
c        standard fortran subroutine
c
c        parameters
c         on entry
c            f      - real*8
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            bound  - real*8
c                     finite bound of integration range
c                     (has no meaning if interval is doubly-infinite)
c
c            inf    - integer
c                     indicating the kind of integration range involved
c                     inf = 1 corresponds to  (bound,+infinity),
c                     inf = -1            to  (-infinity,bound),
c                     inf = 2             to (-infinity,+infinity).
c
c            epsabs - real*8
c                     absolute accuracy requested
c            epsrel - real*8
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c
c         on return
c            result - real*8
c                     approximation to the integral
c
c            abserr - real*8
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine. the
c                             estimates for result and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however, if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties. if
c                             the position of a local difficulty can be
c                             determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is assumed that the requested tolerance
c                             cannot be achieved, and that the returned
c                             result is the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
c                              or limit.lt.1 or leniw.lt.limit*4.
c                             result, abserr, neval, last are set to
c                             zero. exept when limit or leniw is
c                             invalid, iwork(1), work(limit*2+1) and
c                             work(limit*3+1) are set to zero, work(1)
c                             is set to a and work(limit+1) to b.
c
c         dimensioning parameters
c            limit - integer
c                    dimensioning parameter for iwork
c                    limit determines the maximum number of subintervals
c                    in the partition of the given integration interval
c                    (a,b), limit.ge.1.
c                    if limit.lt.1, the routine will end with ier = 6.
c
c            lenw  - integer
c                    dimensioning parameter for work
c                    lenw must be at least limit*4.
c                    if lenw.lt.limit*4, the routine will end
c                    with ier = 6.
c
c            last  - integer
c                    on return, last equals the number of subintervals
c                    produced in the subdivision process, which
c                    determines the number of significant elements
c                    actually in the work arrays.
c
c         work arrays
c            iwork - integer
c                    vector of dimension at least limit, the first
c                    k elements of which contain pointers
c                    to the error estimates over the subintervals,
c                    such that work(limit*3+iwork(1)),... ,
c                    work(limit*3+iwork(k)) form a decreasing
c                    sequence, with k = last if last.le.(limit/2+2), and
c                    k = limit+1-last otherwise
c
c            work  - real*8
c                    vector of dimension at least lenw
c                    on return
c                    work(1), ..., work(last) contain the left
c                     end points of the subintervals in the
c                     partition of (a,b),
c                    work(limit+1), ..., work(limit+last) contain
c                     the right end points,
c                    work(limit*2+1), ...,work(limit*2+last) contain the
c                     integral approximations over the subintervals,
c                    work(limit*3+1), ..., work(limit*3)
c                     contain the error estimates.
c***references  (none)
c***routines called  dqagie,xerror
c***end prologue  dqagi
c
      real*8 abserr,bound,epsabs,epsrel,f,result,work
      integer ier,inf,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
c
      dimension iwork(limit),work(lenw)
c
      external f
c
c         check validity of limit and lenw.
c
c***first executable statement  dqagi
      ier = 6
      neval = 0
      last = 0
      result = 0.0d+00
      abserr = 0.0d+00
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
c
c         prepare call for dqagie.
c
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
c
      call dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,
     *  neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
c
c         call error handler if necessary.
c
       lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror('abnormal return from dqagi',26,ier,lvl)
      return
      end
      subroutine dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr,
     *   neval,ier,alist,blist,rlist,elist,iord,last)
c***begin prologue  dqagie
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a1,h2a4a1
c***keywords  automatic integrator, infinite intervals,
c             general-purpose, transformation, extrapolation,
c             globally adaptive
c***author  piessens,robert,appl. math & progr. div - k.u.leuven
c           de doncker,elise,appl. math & progr. div - k.u.leuven
c***purpose  the routine calculates an approximation result to a given
c            integral   i = integral of f over (bound,+infinity)
c            or i = integral of f over (-infinity,bound)
c            or i = integral of f over (-infinity,+infinity),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i))
c***description
c
c integration over infinite intervals
c standard fortran subroutine
c
c            f      - real*8
c                     function subprogram defining the integrand
c                     function f(x). the actual name for f needs to be
c                     declared e x t e r n a l in the driver program.
c
c            bound  - real*8
c                     finite bound of integration range
c                     (has no meaning if interval is doubly-infinite)
c
c            inf    - real*8
c                     indicating the kind of integration range involved
c                     inf = 1 corresponds to  (bound,+infinity),
c                     inf = -1            to  (-infinity,bound),
c                     inf = 2             to (-infinity,+infinity).
c
c            epsabs - real*8
c                     absolute accuracy requested
c            epsrel - real*8
c                     relative accuracy requested
c                     if  epsabs.le.0
c                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                     the routine will end with ier = 6.
c
c            limit  - integer
c                     gives an upper bound on the number of subintervals
c                     in the partition of (a,b), limit.ge.1
c
c         on return
c            result - real*8
c                     approximation to the integral
c
c            abserr - real*8
c                     estimate of the modulus of the absolute error,
c                     which should equal or exceed abs(i-result)
c
c            neval  - integer
c                     number of integrand evaluations
c
c            ier    - integer
c                     ier = 0 normal and reliable termination of the
c                             routine. it is assumed that the requested
c                             accuracy has been achieved.
c                   - ier.gt.0 abnormal termination of the routine. the
c                             estimates for result and error are less
c                             reliable. it is assumed that the requested
c                             accuracy has not been achieved.
c            error messages
c                     ier = 1 maximum number of subdivisions allowed
c                             has been achieved. one can allow more
c                             subdivisions by increasing the value of
c                             limit (and taking the according dimension
c                             adjustments into account). however,if
c                             this yields no improvement it is advised
c                             to analyze the integrand in order to
c                             determine the integration difficulties.
c                             if the position of a local difficulty can
c                             be determined (e.g. singularity,
c                             discontinuity within the interval) one
c                             will probably gain from splitting up the
c                             interval at this point and calling the
c                             integrator on the subranges. if possible,
c                             an appropriate special-purpose integrator
c                             should be used, which is designed for
c                             handling the type of difficulty involved.
c                         = 2 the occurrence of roundoff error is
c                             detected, which prevents the requested
c                             tolerance from being achieved.
c                             the error may be under-estimated.
c                         = 3 extremely bad integrand behaviour occurs
c                             at some points of the integration
c                             interval.
c                         = 4 the algorithm does not converge.
c                             roundoff error is detected in the
c                             extrapolation table.
c                             it is assumed that the requested tolerance
c                             cannot be achieved, and that the returned
c                             result is the best which can be obtained.
c                         = 5 the integral is probably divergent, or
c                             slowly convergent. it must be noted that
c                             divergence can occur with any other value
c                             of ier.
c                         = 6 the input is invalid, because
c                             (epsabs.le.0 and
c                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                             result, abserr, neval, last, rlist(1),
c                             elist(1) and iord(1) are set to zero.
c                             alist(1) and blist(1) are set to 0
c                             and 1 respectively.
c
c            alist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the left
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            blist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the right
c                     end points of the subintervals in the partition
c                     of the transformed integration range (0,1).
c
c            rlist  - real*8
c                     vector of dimension at least limit, the first
c                      last  elements of which are the integral
c                     approximations on the subintervals
c
c            elist  - real*8
c                     vector of dimension at least limit,  the first
c                     last elements of which are the moduli of the
c                     absolute error estimates on the subintervals
c
c            iord   - integer
c                     vector of dimension limit, the first k
c                     elements of which are pointers to the
c                     error estimates over the subintervals,
c                     such that elist(iord(1)), ..., elist(iord(k))
c                     form a decreasing sequence, with k = last
c                     if last.le.(limit/2+2), and k = limit+1-last
c                     otherwise
c
c            last   - integer
c                     number of subintervals actually produced
c                     in the subdivision process
c
c***references  (none)
c***routines called  d1mach,dqelg,dqk15i,dqpsrt
c***end prologue  dqagie
      real*8 abseps,abserr,alist,area,area1,area12,area2,a1,
     *  a2,blist,boun,bound,b1,b2,correc,defabs,defab1,defab2,
     *  dres,d1mach,elist,epmach,epsabs,epsrel,erlarg,erlast,
     *  errbnd,errmax,error1,error2,erro12,errsum,ertest,f,oflow,resabs,
     *  reseps,result,res3la,rlist,rlist2,small,uflow
      integer id,ier,ierro,inf,iord,iroff1,iroff2,iroff3,jupbnd,k,ksgn,
     *  ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
      logical extrap,noext
c
      dimension alist(limit),blist(limit),elist(limit),iord(limit),
     *  res3la(3),rlist(limit),rlist2(52)
c
      external f
c
c            the dimension of rlist2 is determined by the value of
c            limexp in subroutine dqelg.
c
c
c            list of major variables
c            -----------------------
c
c           alist     - list of left end points of all subintervals
c                       considered up to now
c           blist     - list of right end points of all subintervals
c                       considered up to now
c           rlist(i)  - approximation to the integral over
c                       (alist(i),blist(i))
c           rlist2    - array of dimension at least (limexp+2),
c                       containing the part of the epsilon table
c                       wich is still needed for further computations
c           elist(i)  - error estimate applying to rlist(i)
c           maxerr    - pointer to the interval with largest error
c                       estimate
c           errmax    - elist(maxerr)
c           erlast    - error on the interval currently subdivided
c                       (before that subdivision has taken place)
c           area      - sum of the integrals over the subintervals
c           errsum    - sum of the errors over the subintervals
c           errbnd    - requested accuracy max(epsabs,epsrel*
c                       abs(result))
c           *****1    - variable for the left subinterval
c           *****2    - variable for the right subinterval
c           last      - index for subdivision
c           nres      - number of calls to the extrapolation routine
c           numrl2    - number of elements currently in rlist2. if an
c                       appropriate approximation to the compounded
c                       integral has been obtained, it is put in
c                       rlist2(numrl2) after numrl2 has been increased
c                       by one.
c           small     - length of the smallest interval considered up
c                       to now, multiplied by 1.5
c           erlarg    - sum of the errors over the intervals larger
c                       than the smallest interval considered up to now
c           extrap    - logical variable denoting that the routine
c                       is attempting to perform extrapolation. i.e.
c                       before subdividing the smallest interval we
c                       try to decrease the value of erlarg.
c           noext     - logical variable denoting that extrapolation
c                       is no longer allowed (true-value)
c
c            machine dependent constants
c            ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c           oflow is the largest positive magnitude.
c
c***first executable statement  dqagie
       epmach = d1mach(4)
c
c           test on validity of parameters
c           -----------------------------
c
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
*
      correc = 0.0d+00 ! dummy initialize
      erlarg = 0.0d+00 ! dummy initialize
      ertest = 0.0d+00 ! dummy initialize
      small  = 0.0d+00 ! dummy initialize
*
      if(epsabs.le.0.0d+00.and.epsrel.lt.max(0.5d+02*epmach,0.5d-28))
     *  ier = 6
       if(ier.eq.6) go to 999
c
c
c           first approximation to the integral
c           -----------------------------------
c
c           determine the interval to be mapped onto (0,1).
c           if inf = 2 the integral is computed as i = i1+i2, where
c           i1 = integral of f over (-infinity,0),
c           i2 = integral of f over (0,+infinity).
c
      boun = bound
      if(inf.eq.2) boun = 0.0d+00
      call dqk15i(f,boun,inf,0.0d+00,0.1d+01,result,abserr,
     *  defabs,resabs)
c
c           test on accuracy
c
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      dres = abs(result)
      errbnd = max(epsabs,epsrel*dres)
      if(abserr.le.1.0d+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or.
     *  abserr.eq.0.0d+00) go to 130
c
c           initialization
c           --------------
c
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
      if(dres.ge.(0.1d+01-0.5d+02*epmach)*defabs) ksgn = 1
c
c           main do-loop
c           ------------
c
      do 90 last = 2,limit
c
c           bisect the subinterval with nrmax-th largest error estimate.
c
        a1 = alist(maxerr)
        b1 = 0.5d+00*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        erlast = errmax
        call dqk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1)
        call dqk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2)
c
c           improve previous approximations to integral
c           and error and test for accuracy.
c
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2)go to 15
        if(abs(rlist(maxerr)-area12).gt.0.1d-04*abs(area12)
     *  .or.erro12.lt.0.99d+00*errmax) go to 10
        if(extrap) iroff2 = iroff2+1
        if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = max(epsabs,epsrel*abs(area))
c
c           test for roundoff error and eventually set error flag.
c
        if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
        if(iroff2.ge.5) ierro = 3
c
c           set error flag in the case that the number of
c           subintervals equals limit.
c
        if(last.eq.limit) ier = 1
c
c           set error flag in the case of bad integrand behaviour
c           at some points of the integration range.
c
        if(max(abs(a1),abs(b2)).le.(0.1d+01+0.1d+03*epmach)*
     *  (abs(a2)+0.1d+04*uflow)) ier = 4
c
c           append the newly-created intervals to the list.
c
        if(error2.gt.error1) go to 20
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 30
   20   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
c
c           call subroutine dqpsrt to maintain the descending ordering
c           in the list of error estimates and select the subinterval
c           with nrmax-th largest error estimate (to be bisected next).
c
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
        if(errsum.le.errbnd) go to 115
        if(ier.ne.0) go to 100
        if(last.eq.2) go to 80
        if(noext) go to 90
        erlarg = erlarg-erlast
        if(abs(b1-a1).gt.small) erlarg = erlarg+erro12
        if(extrap) go to 40
c
c           test whether the interval to be bisected next is the
c           smallest interval.
c
        if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
        extrap = .true.
        nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
c
c           the smallest interval has the largest error.
c           before bisecting decrease the sum of the errors over the
c           larger intervals (erlarg) and perform extrapolation.
c
        id = nrmax
        jupbnd = last
        if(last.gt.(2+limit/2)) jupbnd = limit+3-last
        do 50 k = id,jupbnd
          maxerr = iord(nrmax)
          errmax = elist(maxerr)
          if(abs(blist(maxerr)-alist(maxerr)).gt.small) go to 90
          nrmax = nrmax+1
   50   continue
c
c           perform extrapolation.
c
   60   numrl2 = numrl2+1
        rlist2(numrl2) = area
        call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
        ktmin = ktmin+1
        if(ktmin.gt.5.and.abserr.lt.0.1d-02*errsum) ier = 5
        if(abseps.ge.abserr) go to 70
        ktmin = 0
        abserr = abseps
        result = reseps
        correc = erlarg
        ertest = max(epsabs,epsrel*abs(reseps))
        if(abserr.le.ertest) go to 100
c
c            prepare bisection of the smallest interval.
c
   70   if(numrl2.eq.1) noext = .true.
        if(ier.eq.5) go to 100
        maxerr = iord(1)
        errmax = elist(maxerr)
        nrmax = 1
        extrap = .false.
        small = small*0.5d+00
        erlarg = errsum
        go to 90
   80   small = 0.375d+00
        erlarg = errsum
        ertest = errbnd
        rlist2(2) = area
   90 continue
c
c           set final result and error estimate.
c           ------------------------------------
c
  100 if(abserr.eq.oflow) go to 115
      if((ier+ierro).eq.0) go to 110
      if(ierro.eq.3) abserr = abserr+correc
      if(ier.eq.0) ier = 3
      if(result.ne.0.0d+00.and.area.ne.0.0d+00)go to 105
      if(abserr.gt.errsum)go to 115
      if(area.eq.0.0d+00) go to 130
      go to 110
  105 if(abserr/abs(result).gt.errsum/abs(area))go to 115
c
c           test on divergence
c
  110 if(ksgn.eq.(-1).and.max(abs(result),abs(area)).le.
     * defabs*0.1d-01) go to 130
      if(0.1d-01.gt.(result/area).or.(result/area).gt.0.1d+03.
     *or.errsum.gt.abs(area)) ier = 6
      go to 130
c
c           compute global integral sum.
c
  115 result = 0.0d+00
      do 120 k = 1,last
        result = result+rlist(k)
  120 continue
      abserr = errsum
  130 neval = 30*last-15
      if(inf.eq.2) neval = 2*neval
      if(ier.gt.2) ier=ier-1
  999 return
      end
      subroutine dqelg(n,epstab,result,abserr,res3la,nres)
c***begin prologue  dqelg
c***refer to  dqagie,dqagoe,dqagpe,dqagse
c***routines called  d1mach
c***revision date  830518   (yymmdd)
c***keywords  epsilon algorithm, convergence acceleration,
c             extrapolation
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math & progr. div. - k.u.leuven
c***purpose  the routine determines the limit of a given sequence of
c            approximations, by means of the epsilon algorithm of
c            p.wynn. an estimate of the absolute error is also given.
c            the condensed epsilon table is computed. only those
c            elements needed for the computation of the next diagonal
c            are preserved.
c***description
c
c           epsilon algorithm
c           standard fortran subroutine
c           real*8 version
c
c           parameters
c              n      - integer
c                       epstab(n) contains the new element in the
c                       first column of the epsilon table.
c
c              epstab - real*8
c                       vector of dimension 52 containing the elements
c                       of the two lower diagonals of the triangular
c                       epsilon table. the elements are numbered
c                       starting at the right-hand corner of the
c                       triangle.
c
c              result - real*8
c                       resulting approximation to the integral
c
c              abserr - real*8
c                       estimate of the absolute error computed from
c                       result and the 3 previous results
c
c              res3la - real*8
c                       vector of dimension 3 containing the last 3
c                       results
c
c              nres   - integer
c                       number of calls to the routine
c                       (should be zero at first call)
c
c***end prologue  dqelg
c
      real*8 abserr,delta1,delta2,delta3,d1mach,
     *  epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3,
     *  oflow,res,result,res3la,ss,tol1,tol2,tol3
      integer i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm,nres,num
      dimension epstab(52),res3la(3)
c
c           list of major variables
c           -----------------------
c
c           e0     - the 4 elements on which the computation of a new
c           e1       element in the epsilon table is based
c           e2
c           e3                 e0
c                        e3    e1    new
c                              e2
c           newelm - number of elements to be computed in the new
c                    diagonal
c           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
c           result - the element in the new diagonal with least value
c                    of error
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           oflow is the largest positive magnitude.
c           limexp is the maximum number of elements the epsilon
c           table can contain. if this number is reached, the upper
c           diagonal of the epsilon table is deleted.
c
c***first executable statement  dqelg
      epmach = d1mach(4)
      oflow = d1mach(2)
      nres = nres+1
      abserr = oflow
      result = epstab(n)
      if(n.lt.3) go to 100
      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = oflow
      num = n
      k1 = n
      do 40 i = 1,newelm
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
        if(err2.gt.tol2.or.err3.gt.tol3) go to 10
c
c           if e0, e1 and e2 are equal to within machine
c           accuracy, convergence is assumed.
c           result = e2
c           abserr = abs(e1-e0)+abs(e2-e1)
c
        result = res
        abserr = err2+err3
c ***jump out of do-loop
        go to 100
   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = abs(delta1)
        tol1 = max(e1abs,abs(e3))*epmach
c
c           if two elements are very close to each other, omit
c           a part of the table by adjusting the value of n
c
        if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
        ss = 0.1d+01/delta1+0.1d+01/delta2-0.1d+01/delta3
        epsinf = abs(ss*e1)
c
c           test to detect irregular behaviour in the table, and
c           eventually omit a part of the table adjusting the value
c           of n.
c
        if(epsinf.gt.0.1d-03) go to 30
   20   n = i+i-1
c ***jump out of do-loop
        go to 50
c
c           compute a new element and eventually adjust
c           the value of result.
c
   30   res = e1+0.1d+01/ss
        epstab(k1) = res
        k1 = k1-2
        error = err2+abs(res-e2)+err3
        if(error.gt.abserr) go to 40
        abserr = error
        result = res
   40 continue
c
c           shift the table.
c
   50 if(n.eq.limexp) n = 2*(limexp/2)-1
      ib = 1
      if((num/2)*2.eq.num) ib = 2
      ie = newelm+1
      do 60 i=1,ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
   60 continue
      if(num.eq.n) go to 80
      indx = num-n+1
      do 70 i = 1,n
        epstab(i)= epstab(indx)
        indx = indx+1
   70 continue
   80 if(nres.ge.4) go to 90
      res3la(nres) = result
      abserr = oflow
      go to 100
c
c           compute error estimate
c
   90 abserr = abs(result-res3la(3))+abs(result-res3la(2))
     *  +abs(result-res3la(1))
      res3la(1) = res3la(2)
      res3la(2) = res3la(3)
      res3la(3) = result
  100 abserr = max(abserr,0.5d+01*epmach*abs(result))
      return
      end
      subroutine dqk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
c***begin prologue  dqk15i
c***date written   800101   (yymmdd)
c***revision date  830518   (yymmdd)
c***category no.  h2a3a2,h2a4a2
c***keywords  15-point transformed gauss-kronrod rules
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  the original (infinite integration range is mapped
c            onto the interval (0,1) and (a,b) is a part of (0,1).
c            it is the purpose to compute
c            i = integral of transformed integrand over (a,b),
c            j = integral of abs(transformed integrand) over (a,b).
c***description
c
c           integration rule
c           standard fortran subroutine
c           real*8 version
c
c           parameters
c            on entry
c              f      - real*8
c                       fuction subprogram defining the integrand
c                       function f(x). the actual name for f needs to be
c                       declared e x t e r n a l in the calling program.
c
c              boun   - real*8
c                       finite bound of original integration
c                       range (set to zero if inf = +2)
c
c              inf    - integer
c                       if inf = -1, the original interval is
c                                   (-infinity,bound),
c                       if inf = +1, the original interval is
c                                   (bound,+infinity),
c                       if inf = +2, the original interval is
c                                   (-infinity,+infinity) and
c                       the integral is computed as the sum of two
c                       integrals, one over (-infinity,0) and one over
c                       (0,+infinity).
c
c              a      - real*8
c                       lower limit for integration over subrange
c                       of (0,1)
c
c              b      - real*8
c                       upper limit for integration over subrange
c                       of (0,1)
c
c            on return
c              result - real*8
c                       approximation to the integral i
c                       result is computed by applying the 15-point
c                       kronrod rule(resk) obtained by optimal addition
c                       of abscissae to the 7-point gauss rule(resg).
c
c              abserr - real*8
c                       estimate of the modulus of the absolute error,
c                       which should equal or exceed abs(i-result)
c
c              resabs - real*8
c                       approximation to the integral j
c
c              resasc - real*8
c                       approximation to the integral of
c                       abs((transformed integrand)-i/(b-a)) over (a,b)
c
c***references  (none)
c***routines called  d1mach
c***end prologue  dqk15i
c
      real*8 a,absc,absc1,absc2,abserr,b,boun,centr,dinf,
     *  d1mach,epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,
     *  resabs,resasc,resg,resk,reskh,result,tabsc1,tabsc2,uflow,wg,wgk,
     *  xgk
      integer inf,j
      external f
c
      dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(8)
c
c           the abscissae and weights are supplied for the interval
c           (-1,1).  because of symmetry only the positive abscissae and
c           their corresponding weights are given.
c
c           xgk    - abscissae of the 15-point kronrod rule
c                    xgk(2), xgk(4), ... abscissae of the 7-point
c                    gauss rule
c                    xgk(1), xgk(3), ...  abscissae which are optimally
c                    added to the 7-point gauss rule
c
c           wgk    - weights of the 15-point kronrod rule
c
c           wg     - weights of the 7-point gauss rule, corresponding
c                    to the abscissae xgk(2), xgk(4), ...
c                    wg(1), wg(3), ... are set to zero.
c
      data wg(1) / 0.0d0 /
      data wg(2) / 0.1294849661 6886969327 0611432679 082d0 /
      data wg(3) / 0.0d0 /
      data wg(4) / 0.2797053914 8927666790 1467771423 780d0 /
      data wg(5) / 0.0d0 /
      data wg(6) / 0.3818300505 0511894495 0369775488 975d0 /
      data wg(7) / 0.0d0 /
      data wg(8) / 0.4179591836 7346938775 5102040816 327d0 /
c
      data xgk(1) / 0.9914553711 2081263920 6854697526 329d0 /
      data xgk(2) / 0.9491079123 4275852452 6189684047 851d0 /
      data xgk(3) / 0.8648644233 5976907278 9712788640 926d0 /
      data xgk(4) / 0.7415311855 9939443986 3864773280 788d0 /
      data xgk(5) / 0.5860872354 6769113029 4144838258 730d0 /
      data xgk(6) / 0.4058451513 7739716690 6606412076 961d0 /
      data xgk(7) / 0.2077849550 0789846760 0689403773 245d0 /
      data xgk(8) / 0.0000000000 0000000000 0000000000 000d0 /
c
      data wgk(1) / 0.0229353220 1052922496 3732008058 970d0 /
      data wgk(2) / 0.0630920926 2997855329 0700663189 204d0 /
      data wgk(3) / 0.1047900103 2225018383 9876322541 518d0 /
      data wgk(4) / 0.1406532597 1552591874 5189590510 238d0 /
      data wgk(5) / 0.1690047266 3926790282 6583426598 550d0 /
      data wgk(6) / 0.1903505780 6478540991 3256402421 014d0 /
      data wgk(7) / 0.2044329400 7529889241 4161999234 649d0 /
      data wgk(8) / 0.2094821410 8472782801 2999174891 714d0 /
c
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the interval
c           hlgth  - half-length of the interval
c           absc*  - abscissa
c           tabsc* - transformed abscissa
c           fval*  - function value
c           resg   - result of the 7-point gauss formula
c           resk   - result of the 15-point kronrod formula
c           reskh  - approximation to the mean value of the transformed
c                    integrand over (a,b), i.e. to i/(b-a)
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqk15i
      epmach = d1mach(4)
      uflow = d1mach(1)
      dinf = min(1,inf)
c
      centr = 0.5d+00*(a+b)
      hlgth = 0.5d+00*(b-a)
      tabsc1 = boun+dinf*(0.1d+01-centr)/centr
      fval1 = f(tabsc1)
      if(inf.eq.2) fval1 = fval1+f(-tabsc1)
      fc = (fval1/centr)/centr
c
c           compute the 15-point kronrod approximation to
c           the integral, and estimate the error.
c
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
        if(inf.eq.2) fval1 = fval1+f(-tabsc1)
        if(inf.eq.2) fval2 = fval2+f(-tabsc2)
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
      if(resasc.ne.0.0d+00.and.abserr.ne.0.d0) abserr = resasc*
     * min(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = max
     * ((epmach*0.5d+02)*resabs,abserr)
      return
      end
