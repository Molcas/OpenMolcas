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

! In this routine H_0 in RASSI basis is constructed, possibly with external perturbation added on.
subroutine RasH0(nB)

use qmstat_global, only: AddExt, AvRed, BigT, ExtLabel, HmatState, iCompExt, iPrint, MoAveRed, nExtAddOns, nRedMo, nState, ScalExt
use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nB
integer(kind=iwp) :: i, iExt, iopt, irc, iS1, iS2, iSmLbl, kaunter, Lu_One, nBTri, nSize
real(kind=wp) :: Element
real(kind=wp), allocatable :: AOG(:), AOx(:), AUX(:,:), DiagH0(:), MOG(:), MOx(:), SqAO(:,:), SqMO(:,:)
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: Ddot_
#include "warnings.h"

nBTri = nTri_Elem(nB)
if (.not. AddExt) then
  call mma_allocate(DiagH0,nState,label='DiagH0')
  do i=1,nState
    DiagH0(i) = HmatState(nTri_Elem(i))
  end do
  write(u6,*) '     -----RASSI H_0 eigenvalues:'
  write(u6,99) DiagH0(:)
  call mma_deallocate(DiagH0)
else

  ! Collect one-electron perturbations.

  Lu_One = IsFreeUnit(49)
  iopt = 0
  call OpnOne(irc,iopt,'ONEINT',Lu_One)
  call mma_allocate(AOx,nBTri,label='AOExt')
  do iExt=1,nExtAddOns
    irc = -1
    iopt = ibset(ibset(0,sNoOri),sNoNuc)
    iSmLbl = 0
    call RdOne(irc,iopt,ExtLabel(iExt),iCompExt(iExt),AOx,iSmLbl)
    AOx(:) = AOx*ScalExt(iExt)
    if (irc /= 0) then
      write(u6,*)
      write(u6,*) 'ERROR when reading ',ExtLabel(iExt),'.'
      write(u6,*) 'Have Seward computed this integral?'
      call Quit(_RC_IO_ERROR_READ_)
    end if

    ! We need to know in which basis the TDM is and then transform
    ! the one-electron integrals to RASSI-basis.

    if (.not. MoAveRed) then
      call mma_allocate(AOG,nBTri,label='Transition')
      kaunter = 0
      do iS1=1,nState
        do iS2=1,iS1
          kaunter = kaunter+1
          AOG(:) = BigT(:,kaunter)
          Element = Ddot_(nBTri,AOG,1,AOx,1)
          HmatState(kaunter) = HmatState(kaunter)+Element
        end do
      end do
      call mma_deallocate(AOG)
    else
      nSize = nTri_Elem(nRedMO)
      call mma_allocate(MOG,nSize,label='Transition')
      call mma_allocate(AUX,nRedMO,nB,label='AUX')
      call mma_allocate(SqAO,nB,nB,label='SquareAO')
      call mma_allocate(SqMO,nRedMO,nRedMO,label='SquareMO')
      call mma_allocate(MOx,nSize,label='MOExt')
      call Square(AOx,SqAO,1,nB,nB)
      call Dgemm_('T','N',nRedMO,nB,nB,One,AvRed,nB,SqAO,nB,Zero,AUX,nRedMO)
      call Dgemm_('N','N',nRedMO,nRedMO,nB,One,AUX,nRedMO,AvRed,nB,Zero,SqMO,nRedMO)
      call SqToTri_Q(SqMO,MOx,nRedMO)
      kaunter = 0
      do iS1=1,nState
        ! This was 1,nState before... I think that was a bug, because HMatState is triangular
        do iS2=1,iS1
          kaunter = kaunter+1
          MOG(:) = BigT(1:nSize,kaunter)
          Element = Ddot_(nSize,MOG,1,MOx,1)
          HmatState(kaunter) = HmatState(kaunter)+Element
        end do
      end do
      call mma_deallocate(MOG)
      call mma_deallocate(AUX)
      call mma_deallocate(SqAO)
      call mma_deallocate(SqMO)
      call mma_deallocate(MOx)
    end if
  end do
  call mma_deallocate(AOx)
  call ClsOne(irc,Lu_One)

  ! If sufficient print level, print HmatState with perturbation added.

  if (iPrint >= 5) then
    write(u6,*)
    call TriPrt('H_0+External perturbation',' ',HmatState,nState)
  end if
end if

return

99 format('            ',9(F12.7,'  '))

end subroutine RasH0
