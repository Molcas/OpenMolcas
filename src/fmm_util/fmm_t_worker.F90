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

module fmm_T_worker

use fmm_global_paras, only: INTK, REALK, LUPRI, ZERO_VECT_TOL, Zero, One, Two, Half
use fmm_utils, only: fmm_quit

implicit none
private
! Public procedures
public :: fmm_get_SPLTSQ_T_matrix, &
          fmm_get_FLTSQ_T_matrix, &
          fmm_get_boundary_T_matrix, &
          fmm_contract_Tq, &
          fmm_postfac_Vff, &
          fmm_scale_vec

contains

!------------------------------------------------------------------------------
! Get sparse (l+j <= LMAX), lower-triangular, square-based T-matrix
!FIXME: would be nice to allow possibility of building rectangular-based too

subroutine fmm_get_SPLTSQ_T_matrix(LMAX,r_ab,T_matrix)

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: r_ab(3)
  real(REALK), intent(out)  :: T_matrix(:,:)

  real(REALK) :: I_sh((1+LMAX)**2)

  call fmm_generate_I(LMAX,r_ab,I_sh)
  call fmm_generate_T(LMAX,fmm_SPLTSQ_JMAX,I_sh,T_matrix)

end subroutine fmm_get_SPLTSQ_T_matrix

function fmm_SPLTSQ_JMAX(L,LMAX)
  implicit none
  integer(INTK), intent(in) :: L, LMAX
  integer(INTK) :: fmm_SPLTSQ_JMAX
  fmm_SPLTSQ_JMAX = LMAX-L
end function fmm_SPLTSQ_JMAX

!------------------------------------------------------------------------------
! Get full (all elements built), lower-triangular, square T-matrix

subroutine fmm_get_FLTSQ_T_matrix(LMAX,r_ab,T_matrix)

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: r_ab(3)
  real(REALK), intent(out)  :: T_matrix(:,:)

  real(REALK) :: I_sh((1+2*LMAX)**2)

  call fmm_generate_I(2*LMAX,r_ab,I_sh)
  call fmm_generate_T(LMAX,fmm_FLTSQ_JMAX,I_sh,T_matrix)

end subroutine fmm_get_FLTSQ_T_matrix

function fmm_FLTSQ_JMAX(L,LMAX)
# include "macros.fh"
  implicit none
  integer(INTK), intent(in) :: L, LMAX
  integer(INTK) :: fmm_FLTSQ_JMAX
  unused_var(L)
  fmm_FLTSQ_JMAX = LMAX
end function fmm_FLTSQ_JMAX

!------------------------------------------------------------------------------
! Get single-column T-matrix for building boundary potential only

subroutine fmm_get_boundary_T_matrix(LMAX,r_ab,T_matrix)

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: r_ab(3)
  real(REALK), intent(out)  :: T_matrix(:,:)

  real(REALK) :: I_sh((1+2*LMAX)**2)

  call fmm_generate_I(2*LMAX,r_ab,I_sh)
  call fmm_generate_boundary_T(LMAX,I_sh,T_matrix)

end subroutine fmm_get_boundary_T_matrix

!------------------------------------------------------------------------------

subroutine fmm_generate_I(LMAX,vector,I_sh)

  ! Subroutine to generate Scaled Irregular solid harmonics
  ! See page 416, "Molecular Electronic Structure Theory".
  !--------------------------------------------------------
  ! I_sh(M,L)  == COS component of I_sh(M,L)
  ! I_sh(-M,L) == SIN component of I_sh(-M,L)
  ! where L and M are positive
  ! Note that I_sh(0,L) == COS component of I_sh(0,L) !
  ! In fact we combine (l,m) into a single index and
  ! we assume the order 0, -1,1,+1, -2,-1,0,+1,+2 ...
  !--------------------------------------------------------

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: vector(3)
  real(REALK), intent(out)  :: I_sh((LMAX+1)**2)

  real(REALK) :: tmp1, tmp2, tmp3
  real(REALK) :: x, y, z, r_2, r_minus2
  integer(INTK) :: L, m, i, j, k, p, q, u, sign_L

  x = vector(1)
  y = vector(2)
  z = vector(3)
  r_2 = x**2+y**2+z**2
  r_minus2 = one/(r_2)

  !FIXME;
  if (r_2 < ZERO_VECT_TOL) then
    write(LUPRI,'(3E25.15)') vector
    call fmm_quit('Why do we try to do a zero T_vector?')
  end if

  if (LMAX == 0) then
    I_sh(1) = sqrt(r_minus2)
    return
  end if

  ! Start with initialization of lowest order terms
  I_sh(1) = sqrt(r_minus2)
  I_sh(2) = -(y)*(r_minus2)*I_sh(1)
  I_sh(3) = (z)*(r_minus2)*I_sh(1)
  I_sh(4) = -(x)*(r_minus2)*I_sh(1)

  ! Now iterate higher order terms
  sign_L = -1
  do L=2,LMAX
    sign_L = -sign_L
    i = (L+1)*(L+1)
    j = L*L+1
    p = j-2*L+1
    q = j-1
    tmp1 = (2*L-1)*r_minus2
    tmp2 = tmp1*y*sign_L
    tmp3 = tmp1*x
    !
    I_sh(i) = tmp2*I_sh(p)-tmp3*I_sh(q)
    I_sh(j) = tmp2*I_sh(q)+tmp3*I_sh(p)
    !
    p = p+L-1
    q = L*(L-3)+3
    u = q+L-2
    k = i-L
    tmp2 = tmp1*z
    do m=0,(L-2)
      tmp3 = (u-m*m)*(r_minus2)
      !
      I_sh(k+m) = tmp2*I_sh(p+m)-tmp3*I_sh(q+m)
      I_sh(k-m) = tmp2*I_sh(p-m)-tmp3*I_sh(q-m)
      ! I_sh(L,0)=0
      !
    end do
    ! Now do (L, L-1) terms to avoid calling elements that do not exist
    m = L-1
    I_sh(k+m) = tmp2*I_sh(p+m)
    I_sh(k-m) = tmp2*I_sh(p-m)
  end do

end subroutine fmm_generate_I

!------------------------------------------------------------------------------

subroutine fmm_generate_T(LMAX,JMAX,I_sh,T_matrix)

  ! Subroutine to generate Real Interaction Matrix for
  ! given point from Irregular Solid Harmonics (I_sh)
  ! cf Helgaker et. al. pp 415
  !-----------------------------------------------------------------------
  !  T_cos_cos(lm,jk) = 2*{I_cos(l+j,m+k) + (-1)^k*I_cos(l+j,m-k)}  etc.
  !
  ! Note we only write the non-zero elements here, and the precise form
  ! depends on the function JMAX passed in
  !-----------------------------------------------------------------------
  !FIXME: can optimise for j=l by noting that if we build just lower
  ! triangular, we need only half the elements there

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: I_sh(:)
  real(REALK), intent(out)  :: T_matrix(:,:)
  integer(INTK), external   :: JMAX

  integer(INTK) :: sign_k, sign_m
  integer(INTK) :: L, J, k, m, q, pp, qq, u

  if (LMAX == 0) then
    T_matrix(1,1) = two*I_sh(1)     ! NB scaling to make T symmetric
    return
  end if

  L_loop: do L=0,LMAX
    qq = L*(L+1)+1
    J_loop: do J=L,JMAX(L,LMAX)
      pp = J*(J+1)+1
      u = 1+(L+J)*(L+J+1)

      sign_m = -1
      positive_m: do m=0,L
        q = qq+m
        sign_m = -sign_m

        ! cos-cos terms
        !---------------

        sign_k = -1
        do k=0,min(m,J)
          sign_k = -sign_k
          ! (m+k)>0 and (m-k)>0
          T_matrix(pp+k,q) = I_sh(u+(m+k))+sign_k*I_sh(u+(m-k))
        end do

        do k=(m+1),J
          ! (m+k)>0 and (m-k)<0
          T_matrix(pp+k,q) = I_sh(u+(m+k))+sign_m*I_sh(u+(-(m-k)))
        end do

        ! cos-sin terms
        !---------------

        ! only need to consider (m+k)=0 separately as cannot
        ! also have (m-k)=0 unless k=0

        do k=(-J),min((-m-1),-1)
          ! (m+k)<0 and (m-k)>0
          T_matrix(pp+k,q) = I_sh(u+(m+k))+sign_m*I_sh(u+(-(m-k)))
        end do

        do k=max(-J,-m),min(-m,-1)
          ! (m+k)=0 and (m-k)>0
          T_matrix(pp+k,q) = sign_m*I_sh(u+(-(m-k)))
        end do

        sign_k = 1
        do k=-1,max((-m+1),(-J)),-1
          sign_k = -sign_k
          ! (m+k)>0 and (m-k)>0
          T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))+sign_m*I_sh(u+(-(m-k)))
        end do

      end do positive_m

      sign_m = 1
      negative_m: do m=-1,-L,-1
        q = qq+m
        sign_m = -sign_m

        ! sin-cos terms
        !---------------

        k = min(J,(-m))
        sign_k = 1
        if (btest(k,0)) sign_k = -1
        ! (m+k)=0 and (m-k)<0
        T_matrix(pp+k,q) = sign_k*I_sh(u+(m-k))

        sign_k = -1
        do k=0,min((-m-1),J)
          sign_k = -sign_k
          ! (m+k)<0 and (m-k)<0
          T_matrix(pp+k,q) = I_sh(u+(m+k))+sign_k*I_sh(u+(m-k))
        end do

        sign_k = -1
        if (btest(J,0)) sign_k = 1
        do k=J,(-m+1),-1
          sign_k = -sign_k
          ! (m+k)>0 and (m-k)<0
          T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))+sign_k*I_sh(u+(m-k))
        end do

        ! sin-sin terms
        !---------------

        sign_k = -1
        if (btest(J,0)) sign_k = 1
        do k=(-J),m
          sign_k = -sign_k
          ! (m+k)<0 and (m-k)>0
          T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))+sign_k*I_sh(u+(m-k))
        end do

        sign_k = 1
        do k=-1,max((-J),(m+1)),-1
          sign_k = -sign_k
          ! (m+k)<0 and (m-k)<0
          T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))+sign_m*I_sh(u+(-(m-k)))
        end do

      end do negative_m
    end do J_loop
  end do L_loop

end subroutine fmm_generate_T

!------------------------------------------------------------------------------

subroutine fmm_generate_boundary_T(LMAX,I_sh,T_matrix)

!FIXME comments!
! This is a simplified version of fmm_generate_T_matrix when building
! only the boundary potential, when only one column of the T_matrix
! is required.

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: I_sh(:)
  real(REALK), intent(out)  :: T_matrix(:,:)

  integer(INTK) :: J, k, pp, u

  if (LMAX == 0) then
    T_matrix(1,1) = two*I_sh(1)     ! NB scaling to make T symmetric
    return
  end if

  do J=0,LMAX
    pp = J*(J+1)+1
    u = 1+(J)*(J+1)
    do k=-J,J
      T_matrix(pp+k,1) = two*I_sh(u+k)
    end do
  end do

end subroutine fmm_generate_boundary_T

!------------------------------------------------------------------------------

function fmm_contract_Tq(LMAX,vect,T_matrix)

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: vect(:)
  real(REALK), intent(in)   :: T_matrix(:,:)
  real(REALK) :: fmm_contract_Tq(size(vect))

  integer(INTK) :: L, p, q, s, u, qmin, qmax

  ! do L=0 terms first
  p = (LMAX+1)*(LMAX+1)
  fmm_contract_Tq(1) = half*dot_product(vect(1:p),T_matrix(1:p,1))
  do s=2,p
    fmm_contract_Tq(s) = vect(1)*T_matrix(s,1)
  end do

  contract: do L=1,LMAX

    u = L*(L+1)+1
    p = (LMAX-L+1)*(LMAX-L+1)

    qmin = u-L
    qmax = min(u+L,p)

    qloop: do q=qmin,qmax
      fmm_contract_Tq(q) = fmm_contract_Tq(q)+dot_product(vect(q:p),T_matrix(q:p,q))
      do s=q+1,p
        fmm_contract_Tq(s) = fmm_contract_Tq(s)+T_matrix(s,q)*vect(q)
      end do
    end do qloop

    ! add on extra post-factors "after" contraction
    fmm_contract_Tq(u) = half*fmm_contract_Tq(u)    ! m=0

  end do contract

end function fmm_contract_Tq

!------------------------------------------------------------------------------

subroutine fmm_scale_vec(LMAX,scl_fact,scale_vec,prefactor)

  implicit none
  integer(INTK), intent(in) :: LMAX
  real(REALK), intent(in)   :: scl_fact
  real(REALK), intent(out)  :: scale_vec(:), prefactor

  integer(INTK) :: L, L2, lo, hi, u
  real(REALK) :: prefactor1, prevec

  prefactor1 = one/scl_fact

  prevec = one
  scale_vec(1) = prevec
  do L=1,LMAX
    prevec = prefactor1*prevec
    L2 = L*L
    lo = L2+1
    hi = L2+2*L+1
    do u=lo,hi
      scale_vec(u) = prevec
    end do
  end do

  if (scl_fact < zero) then
    prefactor = -prefactor1
  else
    prefactor = prefactor1
  end if

end subroutine fmm_scale_vec

!------------------------------------------------------------------------------

subroutine fmm_postfac_Vff(LMAX,Vff_tmp)

  implicit none
  integer(INTK), intent(in)  :: LMAX
  real(REALK), intent(inout) :: Vff_tmp(:)

  integer(INTK) :: l, u
  ! add extra postfactors  (for half**(delta(0,m))
  do l=0,LMAX
    u = l*(l+1)+1
    Vff_tmp(u) = half*Vff_tmp(u)
  end do

end subroutine fmm_postfac_Vff

!------------------------------------------------------------------------------

end module fmm_T_worker
