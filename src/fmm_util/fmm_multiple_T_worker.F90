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
module fmm_multiple_T_worker

use fmm_global_paras, only: INTK, REALK, One, Two, Half

implicit none
private
! Public procedures
public :: fmm_get_SPLTSQ_T_matrices, &
          fmm_get_FLTSQ_T_matrices, &
          fmm_contract_multi_Tq

contains

!------------------------------------------------------------------------------
! Get sparse (l+j <= LMAX), lower-triangular, square-based T-matrix
!FIXME: would be nice to allow possibility of building rectangular-based too

subroutine fmm_get_SPLTSQ_T_matrices(ndim,LMAX,r_ab,T_matrix)

  implicit none
  integer(INTK), intent(in) :: ndim, LMAX
  real(REALK), intent(in)   :: r_ab(:,:)
  real(REALK), intent(out)  :: T_matrix(:,:,:)

  real(REALK) :: I_sh(ndim,(1+LMAX)**2)

  call fmm_generate_I(ndim,LMAX,r_ab,I_sh)
  call fmm_generate_T(LMAX,.false.,I_sh,T_matrix)

end subroutine fmm_get_SPLTSQ_T_matrices

!------------------------------------------------------------------------------
! Get full (all elements built), lower-triangular, square T-matrix

subroutine fmm_get_FLTSQ_T_matrices(ndim,LMAX,r_ab,T_matrix)

  implicit none
  integer(INTK), intent(in) :: ndim, LMAX
  real(REALK), intent(in)   :: r_ab(:,:)
  real(REALK), intent(out)  :: T_matrix(:,:,:)

  real(REALK) :: I_sh(ndim,(1+2*LMAX)**2)

  call fmm_generate_I(ndim,2*LMAX,r_ab,I_sh)
  call fmm_generate_T(LMAX,.true.,I_sh,T_matrix)

end subroutine fmm_get_FLTSQ_T_matrices

!------------------------------------------------------------------------------

subroutine fmm_generate_I(ndim,LMAX,r_ab,I_sh)

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
  integer(INTK), intent(in) :: LMAX, NDIM
  real(REALK), intent(in)   :: r_ab(:,:)
  real(REALK), intent(out)  :: I_sh(:,:)

  real(REALK) :: r_minus2(NDIM)
  integer(INTK) :: L, m, i, j, k, p, q, u, sign_L, n, p1, q1, m1
  real(REALK) :: var1, var2, var3, rmin2

  r_minus2(:) = one/(r_ab(:,1)**2+r_ab(:,2)**2+r_ab(:,3)**2)

  if (LMAX == 0) then
    I_sh(:,1) = sqrt(r_minus2(:))
    return
  end if

  I_sh(:,1) = sqrt(r_minus2(:))
  I_sh(:,2) = -r_ab(:,2)*r_minus2(:)*I_sh(:,1)
  I_sh(:,3) = r_ab(:,3)*r_minus2(:)*I_sh(:,1)
  I_sh(:,4) = -r_ab(:,1)*r_minus2(:)*I_sh(:,1)

  ! Now iterate higher order terms
  sign_L = -1
  do L=2,LMAX

    sign_L = -sign_L
    i = (L+1)*(L+1)
    j = L*L+1
    p = j-2*L+1
    q = j-1
    p1 = p
    q1 = q
    p = p+L-1
    q = L*(L-3)+3
    u = q+L-2
    k = i-L
    m1 = L-1

    do n=1,ndim

      rmin2 = r_minus2(n)
      var1 = (2*L-1)*rmin2
      var2 = var1*r_ab(n,2)*sign_L
      var3 = var1*r_ab(n,1)
      I_sh(n,i) = var2*I_sh(n,p1)-var3*I_sh(n,q1)
      I_sh(n,j) = var2*I_sh(n,q1)+var3*I_sh(n,p1)
      var2 = var1*r_ab(n,3)

      do m=0,L-2
        var3 = (u-m*m)*rmin2
        I_sh(n,k+m) = var2*I_sh(n,p+m)-var3*I_sh(n,q+m)
        I_sh(n,k-m) = var2*I_sh(n,p-m)-var3*I_sh(n,q-m)
      end do

      I_sh(n,k+m1) = var2*I_sh(n,p+m1)
      I_sh(n,k-m1) = var2*I_sh(n,p-m1)

    end do

  end do

end subroutine fmm_generate_I

!------------------------------------------------------------------------------

subroutine fmm_generate_T(LMAX,TOLMAX,I_sh,T_matrix)

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
  logical, intent(in)       :: TOLMAX
  real(REALK), intent(in)   :: I_sh(:,:)
  real(REALK), intent(out)  :: T_matrix(:,:,:)

  integer(INTK) :: sign_m
  integer(INTK) :: L, J, k, m, q, pp, qq, u, Jlast

  if (LMAX == 0) then
    T_matrix(:,1,1) = two*I_sh(:,1)  ! NB scaling to make T symmetric
    return
  end if

  L_loop: do L=0,LMAX
    qq = L*(L+1)+1

    if (tolmax) then
      Jlast = LMAX
    else
      Jlast = LMAX-L
    end if

    J_loop: do J=L,Jlast
      pp = J*(J+1)+1
      u = 1+(L+J)*(L+J+1)

      sign_m = -1
      positive_m: do m=0,L
        q = qq+m
        sign_m = -sign_m

        ! cos-cos terms
        !---------------

        do k=0,min(m,J)
          if (btest(k,0)) then
            T_matrix(:,pp+k,q) = I_sh(:,u+m+k)-I_sh(:,u+m-k)
          else
            T_matrix(:,pp+k,q) = I_sh(:,u+m+k)+I_sh(:,u+m-k)
          end if
        end do
        if (sign_m > 0) then
          do k=m+1,J
            T_matrix(:,pp+k,q) = I_sh(:,u+m+k)+I_sh(:,u-m+k)
          end do
        else
          do k=m+1,J
            T_matrix(:,pp+k,q) = I_sh(:,u+m+k)-I_sh(:,u-m+k)
          end do
        end if

        ! cos-sin terms
        !---------------

        if (sign_m > 0) then
          do k=-J,min(-m-1,-1)
            T_matrix(:,pp+k,q) = I_sh(:,u+(m+k))+I_sh(:,u+(-(m-k)))
          end do
          do k=max(-J,-m),min(-m,-1)
            T_matrix(:,pp+k,q) = I_sh(:,u+(-(m-k)))
          end do
          do k=-1,max(-m+1,-J),-1
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))+I_sh(:,u+(-(m-k)))
            else
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))+I_sh(:,u+(-(m-k)))
            end if
          end do
        else
          do k=-J,min(-m-1,-1)
            T_matrix(:,pp+k,q) = I_sh(:,u+(m+k))-I_sh(:,u+(-(m-k)))
          end do
          do k=max(-J,-m),min(-m,-1)
            T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m-k)))
          end do
          do k=-1,max(-m+1,-J),-1
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))-I_sh(:,u+(-(m-k)))
            else
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))-I_sh(:,u+(-(m-k)))
            end if
          end do
        end if

      end do positive_m

      sign_m = 1
      negative_m: do m=-1,-L,-1
        q = qq+m
        sign_m = -sign_m

        ! sin-cos terms
        !---------------

        k = min(J,-m)
        if (btest(k,0)) then
          T_matrix(:,pp+k,q) = -I_sh(:,u+(m-k))
        else
          T_matrix(:,pp+k,q) = I_sh(:,u+(m-k))
        end if

        do k=0,min((-m-1),J)
          if (btest(k,0)) then
            T_matrix(:,pp+k,q) = I_sh(:,u+(m+k))-I_sh(:,u+(m-k))
          else
            T_matrix(:,pp+k,q) = I_sh(:,u+(m+k))+I_sh(:,u+(m-k))
          end if
        end do

        if (sign_m > 0) then
          do k=J,(-m+1),-1
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))-I_sh(:,u+(m-k))
            else
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))+I_sh(:,u+(m-k))
            end if
          end do
        else
          do k=J,(-m+1),-1
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))-I_sh(:,u+(m-k))
            else
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))+I_sh(:,u+(m-k))
            end if
          end do
        end if

        ! sin-sin terms
        !---------------

        if (sign_m > 0) then
          do k=-J,m
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))-I_sh(:,u+(m-k))
            else
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))+I_sh(:,u+(m-k))
            end if
          end do
          do k=-1,max(-J,m+1),-1
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))+I_sh(:,u+(-(m-k)))
            else
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))+I_sh(:,u+(-(m-k)))
            end if
          end do
        else
          do k=-J,m
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))-I_sh(:,u+(m-k))
            else
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))+I_sh(:,u+(m-k))
            end if
          end do
          do k=-1,max(-J,m+1),-1
            if (btest(k,0)) then
              T_matrix(:,pp+k,q) = -I_sh(:,u+(-(m+k)))-I_sh(:,u+(-(m-k)))
            else
              T_matrix(:,pp+k,q) = I_sh(:,u+(-(m+k)))-I_sh(:,u+(-(m-k)))
            end if
          end do
        end if

      end do negative_m
    end do J_loop
  end do L_loop

end subroutine fmm_generate_T

!------------------------------------------------------------------------------

function fmm_contract_multi_Tq(LMAX,vect,T_mats,ndim)

  implicit none
  integer(INTK), intent(in) :: LMAX, ndim
  real(REALK), intent(in)   :: vect(:)
  real(REALK), intent(in)   :: T_mats(:,:,:)

  real(REALK) :: fmm_contract_multi_Tq(ndim,(LMAX+1)**2)
  integer(INTK) :: L, p, q, r, s, u, qmin, qmax
  real(REALK) :: fac

  ! first do L=0 terms
  p = (LMAX+1)*(LMAX+1)

  fac = half*vect(1)
  do r=1,ndim
    fmm_contract_multi_Tq(r,1) = fac*T_mats(r,1,1)
  end do

  do s=2,p
    fac = half*vect(s)
    do r=1,ndim
      fmm_contract_multi_Tq(r,1) = fmm_contract_multi_Tq(r,1)+fac*T_mats(r,s,1)
    end do
  end do

  fac = vect(1)
  do s=2,p
    do r=1,ndim
      fmm_contract_multi_Tq(r,s) = fac*T_mats(r,s,1)
    end do
  end do

  contract: do L=1,LMAX

    u = L*(L+1)+1
    p = (LMAX-L+1)*(LMAX-L+1)

    qmin = u-L
    qmax = min(u+L,p)

    qloop: do q=qmin,qmax

      do s=q,p
        call DAXPY_(ndim,vect(s),T_mats(1,s,q:),1,fmm_contract_multi_Tq(1,q),1)
      end do

      fac = vect(q)
      do s=q+1,p
        do r=1,ndim
          fmm_contract_multi_Tq(r,s) = fmm_contract_multi_Tq(r,s)+fac*T_mats(r,s,q)
        end do
      end do

    end do qloop

    fmm_contract_multi_Tq(:,u) = half*fmm_contract_multi_Tq(:,u)    ! m=0

  end do contract

end function fmm_contract_multi_Tq

!------------------------------------------------------------------------------

end module fmm_multiple_T_worker
