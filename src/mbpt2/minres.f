************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1995,1999,2003, Michael A. Saunders                    *
************************************************************************

************************************************************************
* This file from MINRES:                                               *
*   http://web.stanford.edu/group/SOL/software/minres/                 *
*                                                                      *
* The software for MINRES (f77 version) is provided by SOL, Stanford   *
* University under the terms of the OSI Common Public License (CPL):   *
* http://www.opensource.org/licenses/cpl1.0.php                        *
************************************************************************
      subroutine MINRES( n, b, r1, r2, v, w, w1, w2, x, y,
     $                   Aprod, Msolve, checkA, precon, shift,
     $                   nout , itnlim, rtol,
     $                   istop, itn, Anorm, Acond, rnorm, ynorm )

      implicit           none
      external           Aprod, Msolve
      integer            n, nout, itnlim, istop, itn
      logical            checkA, precon
      REAL*8   shift, rtol, Anorm, Acond, rnorm, ynorm,
     $                   b(n), r1(n), r2(n),
     $                   v(n), w(n), w1(n), w2(n), x(n), y(n)
!     ------------------------------------------------------------------
!
!     MINRES  is designed to solve the system of linear equations
!
!                Ax = b
!
!     or the least-squares problem
!
!         min || Ax - b ||_2,
!
!     where A is an n by n symmetric matrix and b is a given vector.
!     The matrix A may be indefinite and/or singular.
!
!     1. If A is known to be positive definite, the Conjugate Gradient
!        Method might be preferred, since it requires the same number
!        of iterations as MINRES but less work per iteration.
!
!     2. If A is indefinite but Ax = b is known to have a solution
!        (e.g. if A is nonsingular), SYMMLQ might be preferred,
!        since it requires the same number of iterations as MINRES
!        but slightly less work per iteration.
!
!     The matrix A is intended to be large and sparse.  It is accessed
!     by means of a subroutine call of the form
!     SYMMLQ development:
!
!                call Aprod ( n, x, y )
!
!     which must return the product y = Ax for any given vector x.
!
!
!     More generally, MINRES is designed to solve the system
!
!                (A - shift*I) x = b
!     or
!         min || (A - shift*I) x - b ||_2,
!
!     where  shift  is a specified scalar value.  Again, the matrix
!     (A - shift*I) may be indefinite and/or singular.
!     The work per iteration is very slightly less if  shift = 0.
!
!     Note: If  shift  is an approximate eigenvalue of  A
!     and  b  is an approximate eigenvector,  x  might prove to be
!     a better approximate eigenvector, as in the methods of
!     inverse iteration and/or Rayleigh-quotient iteration.
!     However, we're not yet sure on that -- it may be better
!     to use SYMMLQ.
!
!     A further option is that of preconditioning, which may reduce
!     the number of iterations required.  If M = C C' is a positive
!     definite matrix that is known to approximate  (A - shift*I)
!     in some sense, and if systems of the form  My = x  can be
!     solved efficiently, the parameters precon and Msolve may be
!     used (see below).  When  precon = .true., MINRES will
!     implicitly solve the system of equations
!
!             P (A - shift*I) P' xbar  =  P b,
!
!     i.e.                  Abar xbar  =  bbar
!     where                         P  =  C**(-1),
!                                Abar  =  P (A - shift*I) P',
!                                bbar  =  P b,
!
!     and return the solution       x  =  P' xbar.
!     The associated residual is rbar  =  bbar - Abar xbar
!                                      =  P (b - (A - shift*I)x)
!                                      =  P r.
!
!     In the discussion below, eps refers to the machine precision.
!     eps is computed by MINRES.  A typical value is eps = 2.22d-16
!     for IEEE double-precision arithmetic.
!
!     Parameters
!     ----------
!
!     n       input      The dimension of the matrix A.
!
!     b(n)    input      The rhs vector b.
!
!     r1(n)   workspace
!     r2(n)   workspace
!     v(n)    workspace
!     w(n)    workspace
!     w1(n)   workspace
!     w2(n)   workspace
!
!     x(n)    output     Returns the computed solution  x.
!
!     y(n)    workspace
!
!     Aprod   external   A subroutine defining the matrix A.
!                        For a given vector x, the statement
!
!                              call Aprod ( n, x, y )
!
!                        must return the product y = Ax
!                        without altering the vector x.
!
!     Msolve  external   An optional subroutine defining a
!                        preconditioning matrix M, which should
!                        approximate (A - shift*I) in some sense.
!                        M must be positive definite.
!                        For a given vector x, the statement
!
!                              call Msolve( n, x, y )
!
!                        must solve the linear system My = x
!                        without altering the vector x.
!
!                        In general, M should be chosen so that Abar has
!                        clustered eigenvalues.  For example,
!                        if A is positive definite, Abar would ideally
!                        be close to a multiple of I.
!                        If A or A - shift*I is indefinite, Abar might
!                        be close to a multiple of diag( I  -I ).
!
!                        NOTE.  The program calling MINRES must declare
!                        Aprod and Msolve to be external.
!
!     checkA  input      If checkA = .true., an extra call of Aprod will
!                        be used to check if A is symmetric.  Also,
!                        if precon = .true., an extra call of Msolve
!                        will be used to check if M is symmetric.
!
!     precon  input      If precon = .true., preconditioning will
!                        be invoked.  Otherwise, subroutine Msolve
!                        will not be referenced; in this case the
!                        actual parameter corresponding to Msolve may
!                        be the same as that corresponding to Aprod.
!
!     shift   input      Should be zero if the system Ax = b is to be
!                        solved.  Otherwise, it could be an
!                        approximation to an eigenvalue of A, such as
!                        the Rayleigh quotient b'Ab / (b'b)
!                        corresponding to the vector b.
!                        If b is sufficiently like an eigenvector
!                        corresponding to an eigenvalue near shift,
!                        then the computed x may have very large
!                        components.  When normalized, x may be
!                        closer to an eigenvector than b.
!
!     nout    input      A file number.
!                        If nout .gt. 0, a summary of the iterations
!                        will be printed on unit nout.
!
!     itnlim  input      An upper limit on the number of iterations.
!
!     rtol    input      A user-specified tolerance.  MINRES terminates
!                        if it appears that norm(rbar) is smaller than
!                              rtol * norm(Abar) * norm(xbar),
!                        where rbar is the transformed residual vector,
!                              rbar = bbar - Abar xbar.
!
!                        If shift = 0 and precon = .false., MINRES
!                        terminates if norm(b - A*x) is smaller than
!                              rtol * norm(A) * norm(x).
!
!     istop   output     An integer giving the reason for termination...
!
!              -1        beta2 = 0 in the Lanczos iteration; i.e. the
!                        second Lanczos vector is zero.  This means the
!                        rhs is very special.
!                        If there is no preconditioner, b is an
!                        eigenvector of A.
!                        Otherwise (if precon is true), let My = b.
!                        If shift is zero, y is a solution of the
!                        generalized eigenvalue problem Ay = lambda My,
!                        with lambda = alpha1 from the Lanczos vectors.
!
!                        In general, (A - shift*I)x = b
!                        has the solution         x = (1/alpha1) y
!                        where My = b.
!
!               0        b = 0, so the exact solution is x = 0.
!                        No iterations were performed.
!
!               1        Norm(rbar) appears to be less than
!                        the value  rtol * norm(Abar) * norm(xbar).
!                        The solution in  x  should be acceptable.
!
!               2        Norm(rbar) appears to be less than
!                        the value  eps * norm(Abar) * norm(xbar).
!                        This means that the residual is as small as
!                        seems reasonable on this machine.
!
!               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!                        which should indicate that x has essentially
!                        converged to an eigenvector of A
!                        corresponding to the eigenvalue shift.
!
!               4        Acond (see below) has exceeded 0.1/eps, so
!                        the matrix Abar must be very ill-conditioned.
!                        x may not contain an acceptable solution.
!
!               5        The iteration limit was reached before any of
!                        the previous criteria were satisfied.
!
!               6        The matrix defined by Aprod does not appear
!                        to be symmetric.
!                        For certain vectors y = Av and r = Ay, the
!                        products y'y and r'v differ significantly.
!
!               7        The matrix defined by Msolve does not appear
!                        to be symmetric.
!                        For vectors satisfying My = v and Mr = y, the
!                        products y'y and r'v differ significantly.
!
!               8        An inner product of the form  x' M**(-1) x
!                        was not positive, so the preconditioning matrix
!                        M does not appear to be positive definite.
!
!                        If istop .ge. 5, the final x may not be an
!                        acceptable solution.
!
!     itn     output     The number of iterations performed.
!
!     Anorm   output     An estimate of the norm of the matrix operator
!                        Abar = P (A - shift*I) P',   where P = C**(-1).
!
!     Acond   output     An estimate of the condition of Abar above.
!                        This will usually be a substantial
!                        under-estimate of the true condition.
!
!     rnorm   output     An estimate of the norm of the final
!                        transformed residual vector,
!                           P (b  -  (A - shift*I) x).
!
!     ynorm   output     An estimate of the norm of xbar.
!                        This is sqrt( x'Mx ).  If precon is false,
!                        ynorm is an estimate of norm(x).
!     ------------------------------------------------------------------
!
!
!     MINRES is an implementation of the algorithm described in
!     the following reference:
!
!     C. C. Paige and M. A. Saunders (1975),
!     Solution of sparse indefinite systems of linear equations,
!     SIAM J. Numer. Anal. 12(4), pp. 617-629.
!     ------------------------------------------------------------------
!
!
!     MINRES development:
!            1972: First version, similar to original SYMMLQ.
!                  Later lost @#%*!
!        Oct 1995: Tried to reconstruct MINRES from
!                  1995 version of SYMMLQ.
!     30 May 1999: Need to make it more like LSQR.
!                  In middle of major overhaul.
!     19 Jul 2003: Next attempt to reconstruct MINRES.
!                  Seems to need two vectors more than SYMMLQ.  (w1, w2)
!                  Lanczos is now at the top of the loop,
!                  so the operator Aprod is called in just one place
!                  (not counting the initial check for symmetry).
!     22 Jul 2003: Success at last.  Preconditioning also works.
!                  minres.f added to http://www.stanford.edu/group/SOL/.
!
!     FUTURE WORK: A stopping rule is needed for singular systems.
!                  We need to estimate ||Ar|| as in LSQR.  This will be
!                  joint work with Sou Cheng Choi, SCCM, Stanford.
!                  Note that ||Ar|| small => r is a null vector for A.
!
!
!     Michael A. Saunders           na.msaunders@na-net.ornl.gov
!     Department of MS&E            saunders@stanford.edu
!     Stanford University
!     Stanford, CA 94305-4026       (650) 723-1875
!     ------------------------------------------------------------------
!
!
!     Subroutines and functions
!
!     USER       Aprod , Msolve
!     BLAS1      daxpy , dcopy , ddot  , dnrm2  } These are all in
!     Utilities  daxpy2, dload2, dscal2         } the file minresblas.f

!     Functions

      real*8, external :: ddot_

!     Local variables

      REAL*8   alfa  , beta  , beta1 , cs    ,
     $                   dbar  , delta , denom , diag  ,
     $                   eps   , epsa  , epsln , epsr  , epsx  ,
     $                   gamma , gbar  , gmax  , gmin  ,
     $                   oldb  , oldeps, qrnorm, phi   , phibar,
     $                   rhs1  , rhs2  , s     , sn    , t     ,
     $                   tnorm2, ynorm2, z
      integer            i
      logical            debug, prnt

      REAL*8   zero, one, two, ten
      parameter        ( zero = 0.0d+0,  one =  1.0d+0,
     $                   two  = 2.0d+0,  ten = 10.0d+0 )

      character*16       enter, exit
      character*52       msg(-1:8)

      data               enter /' Enter MINRES.  '/,
     $                   exit  /' Exit  MINRES.  '/

      data               msg
     $ / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',
     $   'beta1 = 0.  The exact solution is  x = 0',
     $   'Requested accuracy achieved, as determined by rtol',
     $   'Reasonable accuracy achieved, given eps',
     $   'x has converged to an eigenvector',
     $   'Acond has exceeded 0.1/eps',
     $   'The iteration limit was reached',
     $   'Aprod  does not define a symmetric matrix',
     $   'Msolve does not define a symmetric matrix',
     $   'Msolve does not define a pos-def preconditioner' /
!     ------------------------------------------------------------------

      debug = .false.
      gmax= 99999999
      gmin=-99999999

!     ------------------------------------------------------------------
!     Compute eps, the machine precision.  The call to daxpy is
!     intended to fool compilers that use extra-length registers.
!     31 May 1999: Hardwire eps so the debugger can step thru easily.
!     ------------------------------------------------------------------
      eps    = 2.22d-16    ! Set eps = zero here if you want it computed

      if (eps .le. zero) then
         eps    = two**(-12)
   10       eps    = eps / two
            x(1)   = eps
            y(1)   = one
            call daxpy_( 1, one, x, 1, y, 1 )
            if (y(1) .gt. one) go to 10
         eps    = eps * two
      end if

!     ------------------------------------------------------------------
!     Print heading and initialize.
!     ------------------------------------------------------------------
      if (nout .gt. 0) then
         write(nout, 1000) enter, n, checkA, precon,
     $                     itnlim, rtol, shift
      end if
      istop  = 0
      itn    = 0
      Anorm  = zero
      Acond  = zero
      rnorm  = zero
      ynorm  = zero
      call dload2( n, zero, x )

!     ------------------------------------------------------------------
!     Set up y and v for the first Lanczos vector v1.
!     y  =  beta1 P' v1,  where  P = C**(-1).
!     v is really P' v1.
!     ------------------------------------------------------------------
      call dcopy_( n, b, 1, y , 1 )         ! y  = b
      call dcopy_( n, b, 1, r1, 1 )         ! r1 = b
      if ( precon ) call Msolve( n, b, y )
      beta1  = ddot_( n, b, 1, y, 1 )

      if (beta1 .lt. zero) then    ! M must be indefinite.
         istop = 8
         go to 900
      end if

      if (beta1 .eq. zero) then    ! b = 0 exactly.  Stop with x = 0.
         istop = 0
         go to 900
      end if

      beta1  = sqrt( beta1 )       ! Normalize y to get v1 later.

!     ------------------------------------------------------------------
!     See if Msolve is symmetric.
!     ------------------------------------------------------------------
      if (checkA  .and.  precon) then
         call Msolve( n, y, r2 )
         s      = ddot_( n, y, 1, y, 1 )
         t      = ddot_( n,r1, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333d+0
         if (z .gt. epsa) then
            istop = 7
            go to 900
         end if
      end if

!     ------------------------------------------------------------------
!     See if Aprod  is symmetric.
!     ------------------------------------------------------------------
      if (checkA) then
         call Aprod ( n, y, w )
         call Aprod ( n, w, r2 )
         s      = ddot_( n, w, 1, w, 1 )
         t      = ddot_( n, y, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333d+0
         if (z .gt. epsa) then
            istop = 6
            go to 900
         end if
      end if

!     ------------------------------------------------------------------
!     Initialize other quantities.
!     ------------------------------------------------------------------
      oldb   = zero
      beta   = beta1
      dbar   = zero
      epsln  = zero
      qrnorm = beta1
      phibar = beta1
      rhs1   = beta1
      rhs2   = zero
      tnorm2 = zero
      ynorm2 = zero
      cs     = - one
      sn     = zero
      call dload2( n, zero, w  )        ! w  = 0
      call dload2( n, zero, w2 )        ! w2 = 0
      call dcopy_( n, r1, 1, r2, 1 )    ! r2 = r1

      if (debug) then
         write(6,*) ' '
         write(6,*) 'b    ', b
         write(6,*) 'beta ', beta
         write(6,*) ' '
      end if

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------
  100 itn    = itn   +  1               ! k = itn = 1 first time through
      if (istop .ne. 0) go to 900

      !-----------------------------------------------------------------
      ! Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
      ! The general iteration is similar to the case k = 1 with v0 = 0:
      !
      !   p1      = Operator * v1  -  beta1 * v0,
      !   alpha1  = v1'p1,  '
      !   q2      = p2  -  alpha1 * v1,
      !   beta2^2 = q2'q2,  '
      !   v2      = (1/beta2) q2.
      !
      ! Again, y = betak P vk,  where  P = C**(-1).
      ! .... more description needed.
      !-----------------------------------------------------------------
      s      = one / beta            ! Normalize previous vector (in y).
      call dscal2( n, s, y, v )      ! v = vk if P = I

      call Aprod ( n, v, y )
      call daxpy_( n, (- shift), v, 1, y, 1 )
      if (itn .ge. 2) then
         call daxpy_( n, (- beta/oldb), r1, 1, y, 1 )
      end if

      alfa   = ddot_( n, v, 1, y, 1 )     ! alphak

      call daxpy_( n, (- alfa/beta), r2, 1, y, 1 )
      call dcopy_( n, r2, 1, r1, 1 )
      call dcopy_( n,  y, 1, r2, 1 )
      if ( precon ) call Msolve( n, r2, y )

      oldb   = beta                        ! oldb = betak
      beta   = ddot_( n, r2, 1, y, 1 )    ! beta = betak+1^2
      if (beta .lt. zero) then
         istop = 6
         go to 900
      end if

      beta   = sqrt( beta )                ! beta = betak+1
      tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

      if (itn .eq. 1) then                 ! Initialize a few things.
         if (beta/beta1 .le. ten*eps) then ! beta2 = 0 or ~ 0.
            istop = -1                     ! Terminate later.
         end if
        !tnorm2 = alfa**2
         gmax   = abs( alfa )              ! alpha1
         gmin   = gmax                     ! alpha1
      end if

      ! Apply previous rotation Qk-1 to get
      !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
      !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

      oldeps = epsln
      delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
      gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
      epsln  =               sn * beta ! epsln2 = 0         epslnk+1
      dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

      ! Compute the next plane rotation Qk

      gamma  = sqrt( gbar**2 + beta**2 )   ! gammak
      cs     = gbar / gamma                ! ck
      sn     = beta / gamma                ! sk
      phi    = cs * phibar                 ! phik
      phibar = sn * phibar                 ! phibark+1

      if (debug) then
         write(6,*) ' '
         write(6,*) 'v    ', v
         write(6,*) 'alfa ', alfa
         write(6,*) 'beta ', beta
         write(6,*) 'gamma', gamma
         write(6,*) 'delta', delta
         write(6,*) 'gbar ', gbar
         write(6,*) 'epsln', epsln
         write(6,*) 'dbar ', dbar
         write(6,*) 'phi  ', phi
         write(6,*) 'phiba', phibar
         write(6,*) ' '
      end if

      ! Update  x.

      denom = one/gamma

      do i = 1, n
         w1(i) = w2(i)
         w2(i) = w(i)
         w(i)  = ( v(i) - oldeps*w1(i) - delta*w2(i) ) * denom
         x(i)  =   x(i) +   phi * w(i)
      end do

      ! Go round again.

      gmax   = max( gmax, gamma )
      gmin   = min( gmin, gamma )
      z      = rhs1 / gamma
      ynorm2 = z**2  +  ynorm2
      rhs1   = rhs2  -  delta * z
      rhs2   =       -  epsln * z

      ! Estimate various norms and test for convergence.

      Anorm  = sqrt( tnorm2 )
      ynorm  = sqrt( ynorm2 )
      epsa   = Anorm * eps
      epsx   = Anorm * ynorm * eps
      epsr   = Anorm * ynorm * rtol
      diag   = gbar
      if (diag .eq. zero) diag = epsa

      qrnorm = phibar
      rnorm  = qrnorm

      ! Estimate  cond(A).
      ! In this version we look at the diagonals of  R  in the
      ! factorization of the lower Hessenberg matrix,  Q * H = R,
      ! where H is the tridiagonal matrix from Lanczos with one
      ! extra row, beta(k+1) e_k^T.

      Acond  = gmax / gmin

      ! See if any of the stopping criteria are satisfied.
      ! In rare cases, istop is already -1 from above (Abar = const*I).

      if (istop .eq. 0) then
         if (itn    .ge. itnlim    ) istop = 5
         if (Acond  .ge. 0.1d+0/eps) istop = 4
         if (epsx   .ge. beta1     ) istop = 3
         if (qrnorm .le. epsx      ) istop = 2
         if (qrnorm .le. epsr      ) istop = 1
      end if


      ! See if it is time to print something.

      if (nout .gt. 0) then
         prnt   = .false.
         if (n      .le. 40         ) prnt = .true.
         if (itn    .le. 10         ) prnt = .true.
         if (itn    .ge. itnlim - 10) prnt = .true.
         if (mod(itn,10)  .eq.     0) prnt = .true.
         if (qrnorm .le.  ten * epsx) prnt = .true.
         if (qrnorm .le.  ten * epsr) prnt = .true.
         if (Acond  .ge. 1.0d-2/eps ) prnt = .true.
         if (istop  .ne.  0         ) prnt = .true.

         if ( prnt ) then
            if (    itn     .eq. 1) write(nout, 1200)
            write(nout, 1300) itn, x(1), qrnorm, Anorm, Acond
            if (mod(itn,10) .eq. 0) write(nout, 1500)
         end if
      end if

      go to 100
!     ------------------------------------------------------------------
!     End of main iteration loop.
!     ------------------------------------------------------------------

      ! Display final status.

  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, istop, itn,
     $                     exit, Anorm, Acond,
     $                     exit, rnorm, ynorm
         write(nout, 3000) exit, msg(istop)
      end if

      return


 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b'
     $       / ' n      =', i7, 5x, 'checkA =', l4, 12x,
     $          'precon =', l4
     $       / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,
     $          'shift  =', e23.14)
 1200 format(// 5x, 'itn', 8x, 'x(1)', 10x,
     $         'norm(r)', 3x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, 3e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 5x, 'istop =', i3,   14x, 'itn   =', i8
     $       /     a, 5x, 'Anorm =', e12.4, 5x, 'Acond =', e12.4
     $       /     a, 5x, 'rnorm =', e12.4, 5x, 'ynorm =', e12.4)
 3000 format(      a, 5x, a )

      end ! subroutine MINRES

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dload2( n, const, x )

      implicit           none
      integer            n
      REAL*8   const, x(n)

*     ------------------------------------------------------------------
*     dload2 loads all elements of x with const.
*     ------------------------------------------------------------------

      integer            i

      do i = 1, n
         x(i) = const
      end do

      end ! subroutine dload2

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dscal2( n, a, x, y )

      implicit           none
      integer            n
      REAL*8   a, x(n), y(n)

*     ------------------------------------------------------------------
*     dscal2 sets y = a*x.
*     ------------------------------------------------------------------

      integer            i

      do i = 1, n
         y(i) = a*x(i)
      end do

      end ! subroutine dscal2
