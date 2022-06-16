!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Walter Gautschi                                  *
!               1994, Gene H. Golub                                    *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

module RysScratch
! Compute Rys roots and weights from scratch

use Definitions, only: wp, iwp

implicit none
private

! WhichQuad: which shifted Legendre quadrature to use
!            (number of points would be 24,27,30,34,37,39,42,45,46,50,51,54,56,59,61,63,66,68,70,73,300)
! TAsymp: asymptotic limit for the t parameter, for values larger than
!         this, the scaled Hermite quadrature is accurate enough.
!         These values are quite conservative, with estimated errors
!         below 1e-16.

integer(kind=iwp), parameter :: naux(11) = [30,35,40,45,50,55,60,65,70,75,300], &
                                WhichQuad(21) = [1,1,1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,9,10,size(naux)]
real(kind=wp), parameter :: TAsymp(20) = [39.0_wp,47.0_wp,54.0_wp,60.0_wp,66.0_wp,72.0_wp,78.0_wp,83.0_wp,89.0_wp,94.0_wp,99.0_wp, &
                                          104.0_wp,109.0_wp,115.0_wp,120.0_wp,125.0_wp,130.0_wp,134.0_wp,139.0_wp,144.0_wp]
real(kind=wp), allocatable :: Leg_r(:,:), Leg_w(:,:)

public :: RysRtsWgh, SetAux, UnSetAux

contains

! Compute and store roots and weights for a shifted Legendre quadrature,
! used for boot-strapping Rys roots and weights.
! Different sets of roots and weights are computed

subroutine SetAux(eps)

  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One, Four, Half, Quart
  use Definitions, only: u6

  real(kind=wp), intent(in) :: eps
  integer(kind=iwp) :: Err, i, j, maux
  real(kind=wp), allocatable :: a(:), b(:)
  integer(kind=iwp), parameter :: nquad = size(naux)

  if (allocated(Leg_r)) return
  maux = maxval(naux)
  call mma_allocate(Leg_r,maux,nquad,label='Leg_r')
  call mma_allocate(Leg_w,maux,nquad,label='Leg_w')
  call mma_allocate(a,maux)
  call mma_allocate(b,maux)
  do j=1,nquad
    a(1:naux(j)) = Half
    b(1) = One
    do i=2,naux(j)
      b(i) = Quart/(Four-One/(i-1)**2)
    end do
    call GaussQuad(naux(j),a,b,eps,Leg_r(1,j),Leg_w(1,j),Err)
    if (Err /= 0) then
      write(u6,*) Err
      call WarningMessage(2,'Error in GaussQuad')
      call AbEnd()
    end if
    Leg_r(1:naux(j),j) = Leg_r(1:naux(j),j)**2
  end do
  call mma_deallocate(a)
  call mma_deallocate(b)

end subroutine SetAux

subroutine UnSetAux()

  use stdalloc, only: mma_deallocate

  if (allocated(Leg_r)) call mma_deallocate(Leg_r)
  if (allocated(Leg_w)) call mma_deallocate(Leg_w)

end subroutine UnSetAux

! Compute Rys roots and weights from scratch:
! - For high t values use the asymptotic scaled Hermite quadrature
! - For lower t, get first alpha and beta from the auxiliary Legendre
!   quadrature, then compute roots and weights

subroutine RysRtsWgh(TValues,nT,Roots,Weights,Order)

  use vRys_RW, only: HerR2, HerW2, iHerR2, iHerW2
  use Gateway_global, only: asymptotic_Rys
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Five
  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: nT, Order
  real(kind=wp), intent(in) :: TValues(nT)
  real(kind=wp), intent(out) :: Roots(Order,nT), Weights(Order,nT)
  integer(kind=iwp) :: i, iquad, Err
  real(kind=wp) :: Alpha(Order), Beta(Order), TA
  real(kind=wp), allocatable :: a(:), b(:)
  real(kind=wp), parameter :: eps = 1.0e-16_wp

  if (Order > size(TAsymp)) then
    ! Rough fit for asymptotic T
    TA = 50.0_wp+Five*Order
  else
    TA = TAsymp(Order)
  end if

  do i=1,nT
    if ((TValues(i) > TA) .or. asymptotic_Rys) then
      Roots(:,i) = HerR2(iHerR2(Order):iHerR2(Order)+Order-1)/TValues(i)
      Weights(:,i) = HerW2(iHerW2(Order):iHerW2(Order)+Order-1)/sqrt(TValues(i))
    else
      iquad = WhichQuad(min(Order,size(WhichQuad)))
      call mma_allocate(a,naux(iquad))
      call mma_allocate(b,naux(iquad))
      a(1:naux(iquad)) = Leg_r(1:naux(iquad),iquad)
      b(1:naux(iquad)) = Leg_w(1:naux(iquad),iquad)*exp(-TValues(i)*a(1:naux(iquad)))
      call Lanczos(Order,naux(iquad),a,b,Alpha,Beta,Err)
      if (Err /= 0) then
        write(u6,*) Err
        call WarningMessage(2,'Error in Lanczos')
        call AbEnd()
      end if
      call GaussQuad(Order,Alpha,Beta,eps,Roots(1,i),Weights(1,i),Err)
      if (Err /= 0) then
        write(u6,*) Err
        call WarningMessage(2,'Error in GaussQuad 2')
        call AbEnd()
      end if
      call mma_deallocate(a)
      call mma_deallocate(b)
    end if
  end do

end subroutine RysRtsWgh

!***********************************************************************
! Routines GaussQuad and Lanczos adapted from:
!
! Algorithm 726: ORTHPOL -- A package of Routines for Generating
! Orthogonal Polynomials and Gauss-Type Quadrature Rules
!   Walter Gautschi. ACM Trans. Math. Softw. 20 (1994) 21-62
!   doi:10.1145/174603.174605
!***********************************************************************

subroutine GaussQuad(n,alpha,beta,eps,roots,weights,ierr)
! Given  n  and a measure  dlambda, this routine generates the n-point
! Gaussian quadrature formula
!
!     integral over supp(dlambda) of f(x)dlambda(x)
!
!        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
!
! The nodes are returned as  roots(k)=x(k) and the weights as
! weights(k)=w(k), k=1,2,...,n. The user has to supply the recursion
! coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
! dlambda. The routine computes the nodes as eigenvalues, and the
! weights in term of the first component of the respective normalized
! eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
! It uses a translation and adaptation of the algol procedure  imtql2,
! Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
! by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
! Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
! routine  imtql2.
!
!        Input:  n - - the number of points in the Gaussian quadrature
!                      formula; type integer
!                alpha,beta - arrays of dimension  n  to be filled
!                      with the values of  alpha(k-1), beta(k-1), k=1,2,
!                      ...,n
!                eps - the relative accuracy desired in the nodes
!                      and weights
!
!        Output: roots - array of dimension  n  containing the Gaussian
!                      nodes (in increasing order)  roots(k)=x(k), k=1,2,
!                      ...,n
!                weights - array of dimension  n  containing the
!                      Gaussian weights  weights(k)=w(k), k=1,2,...,n
!                ierr - an error flag equal to  0  on normal return,
!                      equal to  i  if the QR algorithm does not
!                      converge within 30 iterations on evaluating the
!                      i-th eigenvalue, equal to  -1  if  n  is not in
!                      range, and equal to  -2  if one of the beta's is
!                      negative.

  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Two

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(in) :: alpha(n), beta(n), eps
  real(kind=wp), intent(out) :: roots(n), weights(n)
  integer(kind=iwp), intent(out) :: ierr
  integer(kind=iwp) :: i, ii, j, k, l, m, mml
  real(kind=wp) :: b, c, f, g, p, r, s
  real(kind=wp), allocatable :: e(:)
  integer(kind=iwp), parameter :: maxcyc = 30

  ierr = 0
  if (n < 1) then
    ierr = -1
    return
  end if

  ! Initialization

  call mma_allocate(e,n,label='e')
  roots(:) = alpha
  weights(:) = Zero
  if (any(beta < Zero)) then
    ierr = -2
    return
  end if
  e(1:n-1) = sqrt(beta(2:n))
  if (n == 1) then
    weights(1) = beta(1)
    return
  else
    weights(1) = One
    e(n) = Zero
  end if

  ! Loop over roots

  do l=1,n
    do j=1,maxcyc

      ! Look for a small subdiagonal element.

      do m=l,n-1
        if (abs(e(m)) <= eps*(abs(roots(m))+abs(roots(m+1)))) exit
      end do
      p = roots(l)
      if (m == l) exit

      ! Form shift.

      g = (roots(l+1)-p)/(Two*e(l))
      r = sqrt(g*g+One)
      g = roots(m)-p+e(l)/(g+sign(r,g))
      s = One
      c = One
      p = Zero
      mml = m-l

      ! For i=m-1 step -1 until l do ...

      do ii=1,mml
        i = m-ii
        f = s*e(i)
        b = c*e(i)
        if (abs(f) < abs(g)) then
          s = f/g
          r = sqrt(s*s+One)
          e(i+1) = g*r
          c = One/r
          s = s*c
        else
          c = g/f
          r = sqrt(c*c+One)
          e(i+1) = f*r
          s = One/r
          c = c*s
        end if
        g = roots(i+1)-p
        r = (roots(i)-g)*s+Two*c*b
        p = s*r
        roots(i+1) = g+p
        g = c*r-b

        ! Form first component of vector.

        f = weights(i+1)
        weights(i+1) = s*weights(i)+c*f
        weights(i) = c*weights(i)-s*f
      end do
      roots(l) = roots(l)-p
      e(l) = g
      e(m) = Zero
    end do

    ! Set error - no convergence to an eigenvalue after maxcyc iterations.

    if (j > maxcyc) then
      ierr = l
      return
    end if
  end do
  call mma_deallocate(e)

  ! Order eigenvalues and eigenvectors.

  do ii=2,n
    i = ii-1
    k = i
    p = roots(i)
    do j=ii,n
      if (roots(j) < p) then
        k = j
        p = roots(j)
      end if
    end do
    if (k /= i) then
      roots(k) = roots(i)
      roots(i) = p
      p = weights(i)
      weights(i) = weights(k)
      weights(k) = p
    end if
  end do
  weights(:) = beta(1)*weights**2

  return

end subroutine GaussQuad

subroutine Lanczos(n,ncap,x,w,alpha,beta,ierr)
! This routine carries out the same task as the routine  sti, but
! uses the more stable Lanczos method. The meaning of the input
! and output parameters is the same as in the routine  sti. (This
! routine is adapted from the routine RKPW in W.B. Gragg and
! W.J. Harrod, "The numerically stable reconstruction of Jacobi
! matrices from spectral data", Numer. Math. 44, 1984, 317-335.)
!
! Routine sti:
!
! This routine applies "Stieltjes's procedure" (cf. Section 2.1 of
! W. Gautschi, "On generating orthogonal polynomials", SIAM J. Sci.
! Statist. Comput. 3, 1982, 289-317) to generate the recursion
! coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
! (monic) orthogonal polynomials associated with the inner product
!
!     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
!
! The integer  n  must be between  1  and  ncap, inclusive; otherwise,
! there is an error exit with  ierr=1. The results are stored in the
! arrays  alpha, beta.
!
! If there is a threat of underflow or overflow in the calculation
! of the coefficients  alpha(k)  and  beta(k), the routine exits with
! the error flag  ierr  set equal to  -k  (in the case of underflow)
! or  +k  (in the case of overflow), where  k  is the recursion index
! for which the problem occurs. The former [latter] can often be avoided
! by multiplying all weights  w(k)  by a sufficiently large [small]
! scaling factor prior to entering the routine, and, upon exit, divide
! the coefficient  beta(0)  by the same factor.
!
! This routine should be used with caution if  n  is relatively close
! to  ncap, since there is a distinct possibility of numerical
! instability developing. (See W. Gautschi, "Is the recurrence relation
! for orthogonal polynomials always stable?", BIT, 1993, to appear.)
! In that case, the routine  lancz  should be used.

  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One

  integer(kind=iwp), intent(in) :: n, ncap
  real(kind=wp), intent(in) :: x(ncap), w(ncap)
  real(kind=wp), intent(out) :: alpha(n), beta(n)
  integer(kind=iwp), intent(out) :: ierr
  integer(kind=iwp) :: i, k
  real(kind=wp) :: gam, pj, rho, sig, t, tk, tmp, tsig, xlam
  real(kind=wp), allocatable :: p0(:), p1(:)

  ierr = 0
  if ((n <= 0) .or. (n > ncap)) then
    ierr = 1
    return
  end if
  call mma_allocate(p0,ncap,label='p0')
  call mma_allocate(p1,ncap,label='p1')
  p0(:) = x
  p1(:) = Zero
  p1(1) = w(1)
  do i=1,ncap-1
    pj = w(i+1)
    gam = One
    sig = Zero
    t = Zero
    xlam = x(i+1)
    do k=1,i+1
      rho = p1(k)+pj
      tmp = gam*rho
      tsig = sig
      if (rho <= Zero) then
        gam = One
        sig = Zero
      else
        gam = p1(k)/rho
        sig = pj/rho
      end if
      tk = sig*(p0(k)-xlam)-gam*t
      p0(k) = p0(k)-(tk-t)
      t = tk
      if (sig <= Zero) then
        pj = tsig*p1(k)
      else
        pj = (t**2)/sig
      end if
      tsig = sig
      p1(k) = tmp
    end do
  end do
  alpha(:) = p0(1:n)
  beta(:) = p1(1:n)
  call mma_deallocate(p0)
  call mma_deallocate(p1)

  return

end subroutine Lanczos

end module RysScratch
