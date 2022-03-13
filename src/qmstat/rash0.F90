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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "numbers.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
#include "warnings.h"
dimension DiagH0(MxState)

nBTri = nB*(nB+1)/2
if (.not. AddExt) then
  kaunter = 0
  do i=1,nState
    do j=1,i
      kaunter = kaunter+1
    end do
    DiagH0(i) = HmatState(kaunter)
  end do
  write(6,*) '     -----RASSI H_0 eigenvalues:'
  write(6,99) (DiagH0(k),k=1,nState)
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
    call dscal_(nBTri,ScalExt(iExt),Work(ipAOx),iONE)
    if (irc /= 0) then
      write(6,*)
      write(6,*) 'ERROR when reading ',ExtLabel(iExt),'.'
      write(6,*) 'Have Seward computed this integral?'
      call Quit(_RC_IO_ERROR_READ_)
    end if

    ! We need to know in which basis the TDM is and then transform
    ! the one-electron integrals to RASSI-basis.

    if (.not. MoAveRed) then
      call GetMem('Transition','Allo','Real',ipAOG,nBTri)
      kaunter = 0
      do iS1=1,nState
        do iS2=1,iS1
          call dcopy_(nBTri,Work(iBigT+nBTri*kaunter),iONE,Work(ipAOG),iONE)
          Element = Ddot_(nBTri,Work(ipAOG),iONE,Work(ipAOx),iONE)
          kaunter = kaunter+1
          HmatState(kaunter) = HmatState(kaunter)+Element
        end do
      end do
      call GetMem('Transition','Free','Real',ipAOG,nBTri)
    else
      nSize = nRedMO*(nRedMO+1)/2
      call GetMem('Transition','Allo','Real',ipMOG,nSize)
      call GetMem('AUX','Allo','Real',iAUX,nRedMO*nB)
      call GetMem('SquareAO','Allo','Real',iSqAO,nB**2)
      call GetMem('SquareMO','Allo','Real',iSqMO,nRedMO**2)
      call GetMem('MOExt','Allo','Real',ipMOx,nSize)
      call Square(Work(ipAOx),Work(iSqAO),iONE,nB,nB)
      call Dgemm_('T','N',nRedMO,nB,nB,ONE,Work(ipAvRed),nB,Work(iSqAO),nB,ZERO,Work(iAUX),nRedMO)
      call Dgemm_('N','N',nRedMO,nRedMO,nB,ONE,Work(iAUX),nRedMO,Work(ipAvRed),nB,ZERO,Work(iSqMO),nRedMO)
      call SqToTri_Q(Work(iSqMO),Work(ipMOx),nRedMO)
      kaunter = 0
      do iS1=1,nState
        do iS2=1,nState
          call dcopy_(nSize,Work(iBigT+nSize*kaunter),iONE,Work(ipMOG),iONE)
          Element = Ddot_(nSize,Work(ipMOG),iONE,Work(ipMOx),iONE)
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
    write(6,*)
    call TriPrt('H_0+External perturbation',' ',HmatState,nState)
  end if
end if

99 format('            ',9(F12.7,'  '))

return

end subroutine RasH0
