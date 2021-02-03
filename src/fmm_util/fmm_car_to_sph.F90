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

module fmm_car_to_sph

use fmm_global_paras, only: INTK, REALK, Zero, One
use fmm_utils, only: fmm_quit
implicit none
private
! public procedures
public :: fmm_init_car_to_sph, &
          fmm_free_car_to_sph, &
          fmm_transform_car_to_sph

real(REALK), allocatable, save :: SphCoef(:,:,:)

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_car_to_sph(LMAX)

  implicit none
  integer(INTK), intent(in) :: LMAX

  if (allocated(SphCoef)) then
    call fmm_quit('mm_car_to_sph not freed correctly!')
  end if

  allocate (SphCoef(0:LMax*2+1,(LMax+1)*(LMax+2)/2,0:LMax))
  call fmm_build_sphcoef(LMAX)

end subroutine fmm_init_car_to_sph

!-------------------------------------------------------------------------------

subroutine fmm_free_car_to_sph()

  implicit none

  deallocate(SphCoef)

end subroutine fmm_free_car_to_sph

!-------------------------------------------------------------------------------

subroutine fmm_build_sphcoef(LMAX)

  implicit none
  integer(INTK), intent(in) :: LMAX

  integer(INTK) :: L, LL, I, M
  integer(INTK) :: IX, IY
  integer(INTK) :: IPX1, IPY1, IPZ1, IPY2, IPZ2, IPX2
  real(REALK)   :: Denom, PreFac, Arg

  ! --- Initialization ---

  SphCoef(:,:,:) = Zero

  ! --- s coefficients ---

  SphCoef(1,1,0) = One
  if (LMax == 0) return

  ! --- p coefficients ---

  SphCoef(1,2,1) = One
  SphCoef(2,3,1) = One
  SphCoef(3,1,1) = One
  if (LMax /= 1) then

  !SphCoef(1,2,2) =  Sqrt3
  !SphCoef(2,5,2) =  Sqrt3
  !SphCoef(3,1,2) = -Half
  !SphCoef(3,4,2) = -Half
  !SphCoef(3,6,2) =  One
  !SphCoef(4,3,2) =  Sqrt3
  !SphCoef(5,1,2) =  Half * Sqrt3
  !SphCoef(5,4,2) = -Half * Sqrt3

  ! --- Higher than d ---

  do LL=2,LMax
    L = LL-1
    Arg = real(L+L+1,kind=REALK)/real(L+L+2,kind=REALK)
    PreFac = sqrt(Arg)

    I = 0
    do IX=L,0,-1
      do IY=L-IX,0,-1
        I = I+1
        IPX1 = I
        IPY1 = I+(L-IX+1)              ! Increase 1 Y
        IPZ1 = I+(L-IX+2)              ! Increase 1 Z

        ! --- Diagonal recurrence ---

        SphCoef(LL+LL+1,IPX1,LL) = SphCoef(LL+LL+1,IPX1,LL)+PreFac*SphCoef(L+L+1,I,L)
        SphCoef(LL+LL+1,IPY1,LL) = SphCoef(LL+LL+1,IPY1,LL)-PreFac*SphCoef(1,I,L)

        SphCoef(1,IPY1,LL) = SphCoef(1,IPY1,LL)+PreFac*SphCoef(L+L+1,I,L)
        SphCoef(1,IPX1,LL) = SphCoef(1,IPX1,LL)+PreFac*SphCoef(1,I,L)

        ! --- Vertical recurrence (1) ---

        do M=-LL+1,LL-1
          Denom = sqrt(real((L+M+1)*(L-M+1),kind=REALK))
          SphCoef(LL+M+1,IPZ1,LL) = SphCoef(LL+M+1,IPZ1,LL)+(real(L+L+1,kind=REALK)/Denom)*SphCoef(L+M+1,I,L)
        end do

      end do   ! IY
    end do   ! IX

    ! --- Vertical recurrence (2) ---

    I = 0
    do IX=L-1,0,-1
      do IY=L-1-IX,0,-1
        I = I+1
        IPX2 = I
        IPY2 = I+(L-1-IX+1)+(L-1-IX+2)   ! Increase 2 Y
        IPZ2 = I+(L-1-IX+2)+(L-1-IX+3)   ! Increase 2 Z

        do M=-LL+1,LL-1
          Arg = real((L+M)*(L-M),kind=REALK)/real((L+M+1)*(L-M+1),kind=REALK)
          PreFac = sqrt(Arg)
          SphCoef(LL+M+1,IPX2,LL) = SphCoef(LL+M+1,IPX2,LL)-PreFac*SphCoef(L+M,I,L-1)
          SphCoef(LL+M+1,IPY2,LL) = SphCoef(LL+M+1,IPY2,LL)-PreFac*SphCoef(L+M,I,L-1)
          SphCoef(LL+M+1,IPZ2,LL) = SphCoef(LL+M+1,IPZ2,LL)-PreFac*SphCoef(L+M,I,L-1)
        end do   ! M

      end do   ! IY
    end do   ! IX

  end do   ! LL

  ! o p-1, p0, p+1 -> px, py, pz

  end if
  SphCoef(:,:,1) = Zero
  SphCoef(1,1,1) = One
  SphCoef(2,2,1) = One
  SphCoef(3,3,1) = One

end subroutine fmm_build_sphcoef

!-------------------------------------------------------------------------------

subroutine fmm_transform_car_to_sph(CarMpole,SphMpole,ndim,LMAX)

  implicit none
  integer(INTK), intent(in) :: ndim, LMAX
  real(REALK), intent(in)   :: CarMpole(ndim, (LMAX+1)*(LMAX+2)/2, 0:LMAX)
  real(REALK), intent(out)  :: SphMpole(ndim, 2*LMAX+1, 0:LMAX)

  integer(INTK) :: MType, IS, IC, IJAO
  real(REALK)  :: Tmp, TempSphCoef, SphCoef_yzx(3,3)

  SphMpole(:,:,:) = Zero
  SphCoef_yzx(:,:) = Zero
  SphCoef_yzx(1,2) = One
  SphCoef_yzx(2,3) = One
  SphCoef_yzx(3,1) = One

  do MType=0,LMAX
    do IS=1,(MType+MType+1)
      do IC=1,(MType+1)*(MType+2)/2

        if(MType /= 1) then
          TempSphCoef = SphCoef(IS,IC,MType)
        else
          TempSphCoef = SphCoef_yzx(IS,IC)
        endif

        do IJAO=1,ndim
          Tmp = TempSphCoef*CarMpole(IJAO,IC,MType)
          SphMpole(IJAO,IS,MType) = SphMpole(IJAO,IS,MType)+Tmp
        end do

      end do
    end do
  end do

end subroutine fmm_transform_car_to_sph

!-------------------------------------------------------------------------------

end module fmm_car_to_sph
