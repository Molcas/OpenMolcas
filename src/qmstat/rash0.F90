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

use qmstat_global, only: AddExt, ExtLabel, HmatState, iBigT, iCompExt, ipAvRed, iPrint, MoAveRed, nExtAddOns, nRedMo, nState, &
                         ScalExt
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp) :: nB
#include "maxi.fh"
#include "WrkSpc.fh"
#include "warnings.h"
integer(kind=iwp) :: i, iAUX, iExt, iopt, ipAOG, ipAOx, ipMOG, ipMOx, irc, iS1, iS2, iSmLbl, iSqAO, iSqMO, j, kaunter, Lu_One, &
                     nBTri, nSize
real(kind=wp) :: Element
real(kind=wp), allocatable :: DiagH0(:)
integer(kind=iwp), external :: IsFreeUnit
real(kind=r8), external :: Ddot_

nBTri = nTri_Elem(nB)
if (.not. AddExt) then
  call mma_allocate(DiagH0,nState,label='DiagH0')
  kaunter = 0
  do i=1,nState
    do j=1,i
      kaunter = kaunter+1
    end do
    DiagH0(i) = HmatState(kaunter)
  end do
  write(u6,*) '     -----RASSI H_0 eigenvalues:'
  write(u6,99) DiagH0(:)
  call mma_deallocate(DiagH0)
else

  ! Collect one-electron perturbations.

  Lu_One = 49
  Lu_One = IsFreeUnit(Lu_One)
  call OpnOne(irc,0,'ONEINT',Lu_One)
  call GetMem('AOExt','Allo','Real',ipAOx,nBTri+4)
  do iExt=1,nExtAddOns
    irc = -1
    iopt = 0
    iSmLbl = 0
    call RdOne(irc,iopt,ExtLabel(iExt),iCompExt(iExt),Work(ipAOx),iSmLbl)
    call dscal_(nBTri,ScalExt(iExt),Work(ipAOx),1)
    if (irc /= 0) then
      write(u6,*)
      write(u6,*) 'ERROR when reading ',ExtLabel(iExt),'.'
      write(u6,*) 'Have Seward computed this integral?'
      call Quit(_RC_IO_ERROR_READ_)
    end if

    ! We need to know in which basis the TDM is and then transform
    ! the one-electron integrals to RASSI-basis.

    if (.not. MoAveRed) then
      call GetMem('Transition','Allo','Real',ipAOG,nBTri)
      kaunter = 0
      do iS1=1,nState
        do iS2=1,iS1
          call dcopy_(nBTri,Work(iBigT+nBTri*kaunter),1,Work(ipAOG),1)
          Element = Ddot_(nBTri,Work(ipAOG),1,Work(ipAOx),1)
          kaunter = kaunter+1
          HmatState(kaunter) = HmatState(kaunter)+Element
        end do
      end do
      call GetMem('Transition','Free','Real',ipAOG,nBTri)
    else
      nSize = nTri_Elem(nRedMO)
      call GetMem('Transition','Allo','Real',ipMOG,nSize)
      call GetMem('AUX','Allo','Real',iAUX,nRedMO*nB)
      call GetMem('SquareAO','Allo','Real',iSqAO,nB**2)
      call GetMem('SquareMO','Allo','Real',iSqMO,nRedMO**2)
      call GetMem('MOExt','Allo','Real',ipMOx,nSize)
      call Square(Work(ipAOx),Work(iSqAO),1,nB,nB)
      call Dgemm_('T','N',nRedMO,nB,nB,One,Work(ipAvRed),nB,Work(iSqAO),nB,Zero,Work(iAUX),nRedMO)
      call Dgemm_('N','N',nRedMO,nRedMO,nB,One,Work(iAUX),nRedMO,Work(ipAvRed),nB,Zero,Work(iSqMO),nRedMO)
      call SqToTri_Q(Work(iSqMO),Work(ipMOx),nRedMO)
      kaunter = 0
      do iS1=1,nState
        do iS2=1,nState
          call dcopy_(nSize,Work(iBigT+nSize*kaunter),1,Work(ipMOG),1)
          Element = Ddot_(nSize,Work(ipMOG),1,Work(ipMOx),1)
          kaunter = kaunter+1
          HmatState(kaunter) = HmatState(kaunter)+Element
        end do
      end do
      call GetMem('Transition','Free','Real',ipMOG,nSize)
      call GetMem('AUX','Free','Real',iAUX,nRedMO*nB)
      call GetMem('SquareAO','Free','Real',iSqAO,nB**2)
      call GetMem('SquareMO','Free','Real',iSqMO,nRedMO**2)
      call GetMem('MOExt','Free','Real',ipMOx,nSize)
    end if
  end do
  call GetMem('AOExt','Free','Real',ipAOx,nBTri+4)
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
