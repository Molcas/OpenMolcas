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

module fmm_shell_pairs

use fmm_global_paras, only: INTK, REALK, LUPRI, fmm_sh_pairs, fmm_basis, Zero, One, Half
use fmm_utils, only: fmm_quit
implicit none
private
! Public procedures
public :: fmm_get_shell_pairs, fmm_free_shell_pairs

type(fmm_sh_pairs), allocatable, target :: sh_pairs(:)

contains

!-------------------------------------------------------------------------------

subroutine fmm_get_shell_pairs(basis,sh_pairs_ptr)

  type(fmm_basis), intent(in) :: basis
  type(fmm_sh_pairs), pointer :: sh_pairs_ptr(:)
  integer(INTK) :: n_pairs

  if (allocated(sh_pairs)) then
    sh_pairs_ptr => sh_pairs
  else
    ! make list of non-vanishing shell pairs
    call fmm_make_shell_pairs(basis,n_pairs)  ! first get n_pairs
    allocate(sh_pairs(n_pairs))
    call fmm_make_shell_pairs(basis,n_pairs)  ! now store list
    sh_pairs_ptr => sh_pairs
    !FIXME add to fmm_stats
    write(LUPRI,*) 'Number of shell pairs =',size(sh_pairs)
  end if

end subroutine fmm_get_shell_pairs

!-------------------------------------------------------------------------------

subroutine fmm_free_shell_pairs()

  implicit none

  if (allocated(sh_pairs)) deallocate(sh_pairs)

end subroutine fmm_free_shell_pairs

!-------------------------------------------------------------------------------

function fmm_extent(ExpPI)

  !use fmm_md4_globals, only: fmm_shell_pair_ThrFac, fmm_X0
  !use fmm_md4_globals, only: fmm_grain_inv, fmm_extent_min
  !use fmm_box_utils, only: fmm_branch

  implicit none

  real(REALK), intent(in) :: ExpPI

  real(REALK), parameter :: root3 = sqrt(3.0_REALK)
  !real(REALK)   :: tmp1, tmp2
  real(REALK) :: fmm_extent
  !integer(INTK) :: bra_min

  !real(REALK), parameter :: erfc_inv = 3.4589_REALK  ! 1e-6
  !real(REALK), parameter :: erfc_inv = 4.0522_REALK  ! 1e-8
  !real(REALK), parameter :: erfc_inv = 4.5728_REALK  ! 1e-10
  !real(REALK), parameter :: erfc_inv = 5.6739_REALK  ! 1e-15

  !fmm_extent = fmm_X0*ExpPI

  ! Extent based on Classical overlap
  !tmp1 = sqrt(ExpPI)*erfc_inv

  ! Artificial extent to keep all near-field exact integrals
  !bra_min = fmm_branch(fmm_extent_min,fmm_grain_inv)
  !tmp2 = half*root3*(bra_min+1)/fmm_grain_inv

  !tmp2 = zero
  !fmm_extent = max(tmp1,tmp2)

  fmm_extent = sqrt(ExpPI*log(1.0e10_REALK))

end function fmm_extent

!-------------------------------------------------------------------------------

subroutine fmm_make_shell_pairs(basis,n_pairs)

  !use fmm_md4_globals, only: fmm_shell_pair_ThrFac

  implicit none

  type(fmm_basis), intent(in) :: basis
  integer(INTK), intent(out)  :: n_pairs

  !real(REALK), parameter :: ThrFac = fmm_shell_pair_ThrFac
  real(REALK), parameter :: ThrFac = 1.0e-12_REALK

  integer(INTK) :: II, JJ, I, J
  integer(INTK) :: IPrim1, JPrim1
  integer(INTK) :: IPrim2, JPrim2, IJPrim, JPTemp2
  real(REALK) :: Acentr(3), Bcentr(3), ABcentr(3), P(3), PM(3), RAB(3)
  real(REALK) :: ExpA, ExpB, ExpP, ExpPI, ExpAR2
  real(REALK) :: R2AB, ExpKAB, tmp_ext, extent

  n_pairs = 0

  Ishel: do II=1,basis%nshells

    Acentr(:) = basis%Centr(:,basis%KAtom(II))
    IPrim1 = basis%KStart(II)
    IPrim2 = IPrim1+basis%KontG(II)-1

    Jshel: do JJ=1,II

      Bcentr(:) = basis%Centr(:,basis%KAtom(JJ))
      JPrim1 = basis%KStart(JJ)
      JPrim2 = JPrim1+basis%KontG(JJ)-1

      ABcentr(:) = Half*(Acentr(:)+Bcentr(:))
      RAB(:) = Acentr(:)-Bcentr(:)
      R2AB = dot_product(RAB,RAB)

      extent = zero
      IJPrim = 0
      do I=IPrim1,IPrim2
        ExpA = basis%Expnt(I)
        ExpAR2 = ExpA*R2AB
        JPTemp2 = JPrim2
        if (II == JJ) JPTemp2 = I
        do J=JPrim1,JPTemp2
          ExpB = basis%Expnt(J)
          ExpP = ExpA+ExpB
          ExpPI = One/ExpP
          ExpKAB = -ExpAR2*ExpB*ExpPI
          if (ExpKAB >= ThrFac) then
            IJPrim = IJPrim+1
            P(:) = (ExpA*Acentr(:)+ExpB*Bcentr(:))*ExpPI
            PM(:) = P(:)-ABcentr(:)
            tmp_ext = fmm_extent(ExpPI)
            tmp_ext = tmp_ext+sqrt(dot_product(PM,PM))
            extent = max(extent,tmp_ext)
          end if
        end do
      end do

      if (IJPrim > 0) then
        n_pairs = n_pairs+1
        if (allocated(sh_pairs)) then
          if (n_pairs > size(sh_pairs)) call fmm_quit('get_sh_pairs')
          sh_pairs(n_pairs)%I = II
          sh_pairs(n_pairs)%J = JJ
          sh_pairs(n_pairs)%extent = extent
          sh_pairs(n_pairs)%centre(:) = ABcentr(:)
        end if
      end if

    end do Jshel

  end do Ishel

end subroutine fmm_make_shell_pairs

!-------------------------------------------------------------------------------

end module fmm_shell_pairs
