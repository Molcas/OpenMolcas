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

module fmm_W_worker

use fmm_global_paras, only: INTK, REALK, LUPRI, ZERO_VECT_TOL, Zero, One, Half
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_get_ltsqr_W_matrix, &
          fmm_get_boundary_W_matrix, &
          fmm_contract_Wq

contains

!-------------------------------------------------------------------------------
! Routine to build the lower-triangular elements of a square W_matrix.
! Only the non-zero elements of the full square matrix are written.
! LMAX is the order of translation:
!  W_lm,jk is built for all l and j <= LMAX
! r_ab is the translation vector (final-initial).
! (note to translate a potential -r_ab should be passed, and W' used)

subroutine fmm_get_ltsqr_W_matrix(LMAX,JMAX,r_ab,W_matrix)

  implicit none
  integer(INTK), intent(in) :: LMAX, JMAX
  real(REALK), intent(in)   :: r_ab(3)
  real(REALK), intent(out)  :: W_matrix(:,:)
  integer(INTK) :: n, m

  n = max(LMAX,JMAX)
  m = min(LMAX,JMAX)
  call fmm_build_ltsqr_W_matrix()

contains

subroutine fmm_build_ltsqr_W_matrix()
  implicit none
  real(REALK) :: R_sh(-n:n,0:n)

  call fmm_generate_R(n,-r_ab,R_sh)
  ! note minus sign above - see TUH eqn.(9.13.58)
  call fmm_generate_W(n,m,R_sh,W_matrix)

end subroutine fmm_build_ltsqr_W_matrix

end subroutine fmm_get_ltsqr_W_matrix

!-------------------------------------------------------------------------------
! Get single-column W-matrix for getting boundary potential only

subroutine fmm_get_boundary_W_matrix(LMAX,r_ab,W_matrix)

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: r_ab(3)
  real(REALK), intent(out)  :: W_matrix(:,:)

  real(REALK) :: R_sh(-LMAX:LMAX,0:LMAX)

  call fmm_generate_R(LMAX,-r_ab,R_sh)
  ! note minus sign above - see TUH eqn.(9.13.58)
  call fmm_generate_boundary_W(LMAX,R_sh,W_matrix)

end subroutine fmm_get_boundary_W_matrix

!-------------------------------------------------------------------------------

subroutine fmm_generate_R(LMAX,point,R_sh)

  ! Subroutine to generate scaled regular solid harmonics
  ! See page 416, "Molecular Electronic Structure Theory"
  !-------------------------------------------------------
  ! R_sh(M,L)  == COS component of R_sh(M,L)
  ! R_sh(-M,L) == SIN component of R_sh(-M,L)
  ! where L and M are positive
  ! Note that R_sh(0,L) == COS component of R_sh(0,L) !
  !-------------------------------------------------------

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: point(3)
  real(REALK), intent(out)  :: R_sh(-LMAX:LMAX,0:LMAX)

  real(REALK) :: x, y, z, r_2, r_minus2
  real(REALK) :: A(6)
  integer(INTK) :: signl, l, m

  A(:) = zero
  x = point(1)
  y = point(2)
  z = point(3)
  r_2 = dot_product(point,point)
  r_minus2 = one/(r_2)

  if (r_2 < ZERO_VECT_TOL*ZERO_VECT_TOL) then
    write(LUPRI,'(3E25.15)') point
    call fmm_quit('ERROR: Why do we try to do zero W-vector.')
  end if

  if (LMAX == 0) then
    R_sh(0,0) = one
    return
  end if

  ! Start with initialization of lowest order terms

  R_sh(0,0) = one       ! cannot assign "sin" terms for m==0
  R_sh(-1,1) = -half*y
  R_sh(0,1) = z         ! cannot assign "sin" terms for m==0
  R_sh(1,1) = -half*x

  ! Now iterate higher order terms

  signl = -1
  lloop: do l=2,LMAX
    signl = -signl
    A(1) = one/(2*l)
    A(2) = A(1)*x
    A(3) = A(1)*y*signl

    ! cos terms:
    R_sh(l,l) = A(3)*R_sh(1-l,l-1)-A(2)*R_sh(l-1,l-1)
    ! sin terms:
    R_sh(-l,l) = A(3)*R_sh(l-1,l-1)+A(2)*R_sh(1-l,l-1)

    A(5) = (2*l-1)*(z)*(r_minus2)
    do m=0,(l-2)
      A(6) = (r_2)/(l*l-m*m)
      ! cos terms:
      R_sh(m,l) = A(6)*(A(5)*R_sh(m,l-1)-R_sh(m,l-2))
      ! sin terms:
      R_sh(-m,l) = A(6)*(A(5)*R_sh(-m,l-1)-R_sh(-m,l-2))
    end do

    ! Now do (l,l-1) terms to avoid calling elements that do not exist
    m = l-1
    ! cos terms:
    R_sh(m,l) = z*R_sh(m,l-1)
    ! sin terms:
    R_sh(-m,l) = z*R_sh(-m,l-1)

  end do lloop

end subroutine fmm_generate_R

!-------------------------------------------------------------------------------

subroutine fmm_generate_W(LMAX,JMAX,R_sh,W_matrix)

  ! Subroutine to generate Real Translation Matrix for
  ! given point from the Regular Solid Harmonics R_sh(m,l)
  ! cf Helgaker et. al. pp 413
  ! Note that W-matrix is lower triangular wrt W(lm,jk) and that
  ! only NON-ZERO elements are added here.
  ! We must therefore zero W_matrix initially, or use a special contractor
  !------------------------------------------------------------------------

  implicit none
  integer(INTK), intent(in) :: LMAX, JMAX
  real(REALK), intent(in)   :: R_sh(-LMAX:LMAX,0:LMAX)
  real(REALK), intent(out)  :: W_matrix(:,:)

  integer(INTK) :: l, m, j, k, a, b, q, pp, qq, lm_min, lm_max
  integer(INTK) :: signk, phase

  ! l==j case first:
  ! diagonal = one
  ! all other off-diagonal elements identically zero
  do j=0,JMAX
    qq = j*(j+1)+1
    do k=-j,j
      q = qq+k
      W_matrix(q,q) = one
    end do
  end do

  if (LMAX == 0) return
  if (JMAX > LMAX) call fmm_quit('error in fmm_generate_W')

  ! remaining l:
  jloop: do j=0,JMAX
    qq = j*(j+1)+1

    ! loop over positive k
    signk = -1
    kloop1: do k=0,j
      signk = -signk
      q = qq+k

      lloop1: do l=(j+1),LMAX
        a = (l-j)
        pp = l*(l+1)+1

        ! cos-cos terms
        !--------------

        ! (m-k) terms
        phase = 1
        do m=(k-1),max((k-a),0),-1     ! reverse order
          ! (m-k)<0
          phase = -phase
          b = -(m-k)
          W_matrix(pp+m,q) = phase*R_sh(b,a)
        end do
        do m=k,(a+k)      ! (a+k)<l
          ! (m-k)>=0
          b = (m-k)
          W_matrix(pp+m,q) = R_sh(b,a)
        end do
        ! (m+k) terms
        do m=0,(a-k)   ! always overlaps (m-k) terms
          ! (m+k)>=0
          b = (m+k)
          W_matrix(pp+m,q) = W_matrix(pp+m,q)+signk*R_sh(b,a)
        end do

        ! sin-cos terms
        !--------------

        ! m==-k first
        m = min(-k,-1)
        W_matrix(pp+m,q) = zero
        ! (m+k) terms
        do m=-(k+a),-(k+1)
          ! (m+k)<0
          b = (m+k)
          W_matrix(pp+m,q) = signk*R_sh(b,a)
        end do
        phase = -signk
        do m=-k+1,min(-1,(a-k))
          ! (m+k)>0
          phase = -phase
          b = -(m+k)
          W_matrix(pp+m,q) = phase*R_sh(b,a)
        end do
        ! (m-k) terms
        do m=(k-a),-1   ! always overlap (m+k) terms
          ! (m-k)<0
          b = (m-k)
          W_matrix(pp+m,q) = W_matrix(pp+m,q)+R_sh(b,a)
        end do

      end do lloop1
    end do kloop1

    ! (half)**(delta(0,k)) correction:
    lm_min = j*(j+2)+2
    lm_max = LMAX*(LMAX+2)+1
    W_matrix(lm_min:lm_max,qq) = half*W_matrix(lm_min:lm_max,qq)   !k==0

    ! loop over negative k
    signk = 1
    kloop2: do k=-1,-j,-1     ! reverse order
      signk = -signk
      q = qq+k

      lloop2: do l=(j+1),LMAX
        a = (l-j)
        pp = l*(l+1)+1

        ! cos-sin terms
        !--------------
        ! (m+k) terms
        do m=max(0,-(a+k)),(-k-1)
          ! (m+k)<0
          b = (m+k)
          W_matrix(pp+m,q) = signk*R_sh(b,a)
        end do
        m = -k
        W_matrix(pp+m,q) = zero
        phase = -signk
        do m=-k+1,(a-k)  !<=l
          ! (m+k)>=0
          phase = -phase
          b = -(m+k)
          W_matrix(pp+m,q) = phase*R_sh(b,a)
        end do
        ! (m-k) terms
        phase = signk
        do m=0,(a+k)    ! always overlap (m+k) terms
          ! (m-k)>0
          phase = -phase
          b = -(m-k)
          W_matrix(pp+m,q) = W_matrix(pp+m,q)-phase*R_sh(b,a)
        end do

        ! sin-sin terms
        !--------------

        ! (m-k) terms
        phase = 1
        do m=(k-1),(k-a),-1    ! reverse order
          ! (m-k)<0
          phase = -phase
          b = -(m-k)
          W_matrix(pp+m,q) = phase*R_sh(b,a)
        end do
        do m=k,min(-1,(a+k))
          ! (m-k)>=0
          b = (m-k)
          W_matrix(pp+m,q) = R_sh(b,a)
        end do
        ! (m+k) terms
        phase = 1
        do m=-1,-(a+k),-1     ! reverse order
          ! (m+k)<0              ! always overlap (m-k) terms
          phase = -phase
          b = -(m+k)
          W_matrix(pp+m,q) = W_matrix(pp+m,q)-phase*R_sh(b,a)
        end do

      end do lloop2
    end do kloop2
  end do jloop

end subroutine fmm_generate_W

!-------------------------------------------------------------------------------

subroutine fmm_generate_boundary_W(LMAX,R_sh,W_matrix)

  ! This is a simplified version of fmm_generate_W_matrix when building
  ! only the boundary potential, when only one column of the W-matrix
  ! is required.

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: R_sh(-LMAX:LMAX,0:LMAX)
  real(REALK), intent(out)  :: W_matrix(:,:)

  integer(INTK) :: l, m, pp

  W_matrix(1,1) = one
  do l=1,LMAX
    pp = l*(l+1)+1
    do m=-l,l
      W_matrix(pp+m,1) = R_sh(m,l)
    end do
  end do

end subroutine fmm_generate_boundary_W

!-------------------------------------------------------------------------------

subroutine fmm_contract_Wq(N_or_T,W_mat,WLDA,vec_in,n,vec_out,m)

  implicit none
  character, intent(in)      :: N_or_T
  integer(INTK), intent(in)  :: n, m, WLDA
  real(REALK), intent(in)    :: W_mat(WLDA,WLDA)
  real(REALK), intent(in)    :: vec_in(n)
  real(REALK), intent(inout) :: vec_out(m)

  integer(INTK) :: q

  if (N_or_T == 'N') then
    do q=1,n
      vec_out(q:m) = vec_out(q:m)+W_mat(q:m,q)*vec_in(q)
    end do
  else
    ! Perform transpose contraction
    do q=1,m
      vec_out(q) = vec_out(q)+dot_product(W_mat(q:n,q),vec_in(q:n))
    end do
  end if

end subroutine fmm_contract_Wq

!-------------------------------------------------------------------------------

end module fmm_W_worker
