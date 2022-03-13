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

subroutine MomentMod(ipRe,ipNRe,iCmo,nBRe,nBNRe,LindMOs,iS1,iS2,First,DiffMax)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: ipRe, ipNRe, iCmo, nBRe, nBNRe, iS1, iS2
logical(kind=iwp) :: LindMOs(MxBas), First
real(kind=wp) :: DiffMax
integer(kind=iwp) :: i, icomp, iopt, ipDx, ipDxM, ipDxRe, ipDxsq, ipDy, ipDyM, ipDyRe, ipDysq, ipDz, ipDzM, ipDzRe, ipDzsq, &
                     ipTEMP, irc, iSmLbl, j, kaunt1, kaunt2, nSize1, nSize2
real(kind=wp) :: Diffx, Diffy, Diffz, DipRe(3), DipNRe(3)
real(kind=r8), external :: Ddot_

if (First .and. (iPrint >= 5)) then
  write(u6,*)
  write(u6,*) '     Modifications of dipoles by renormalization and basis reduction.'
  write(u6,*)
  write(u6,*) '     State pair    |  Difference '
  write(u6,*) '     --------------|---------------------'
  First = .false.
end if

nSize1 = nBNRe*(nBNRe+1)/2
nSize2 = nBRe*(nBRe+1)/2
call GetMem('DipX','Allo','Real',ipDx,nSize1+4)
call GetMem('DipY','Allo','Real',ipDy,nSize1+4)
call GetMem('DipZ','Allo','Real',ipDz,nSize1+4)
call GetMem('DipXre','Allo','Real',ipDxRe,nSize2)
call GetMem('DipYre','Allo','Real',ipDyRe,nSize2)
call GetMem('DipZre','Allo','Real',ipDzRe,nSize2)
call GetMem('DipXsq','Allo','Real',ipDxsq,nBNRe**2)
call GetMem('DipYsq','Allo','Real',ipDysq,nBNRe**2)
call GetMem('DipZsq','Allo','Real',ipDzsq,nBNRe**2)
call GetMem('DipXm','Allo','Real',ipDxM,nBNRe**2)
call GetMem('DipYm','Allo','Real',ipDyM,nBNRe**2)
call GetMem('DipZm','Allo','Real',ipDzM,nBNRe**2)
call GetMem('TEMP','Allo','Real',ipTEMP,nBNRe**2)
irc = -1
iopt = 0
iSmLbl = 0
! X
icomp = 1
call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDx),iSmLbl)
call Square(Work(ipDx),Work(ipDxsq),1,nBNRe,nBNRe)
call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,One,Work(iCmo),nBNRe,Work(ipDxsq),nBNRe,Zero,Work(ipTEMP),nBNRe)
call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,One,Work(ipTEMP),nBNRe,Work(iCmo),nBNRe,Zero,Work(ipDxM),nBNRe)
! Y
icomp = 2
call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDy),iSmLbl)
call Square(Work(ipDy),Work(ipDysq),1,nBNRe,nBNRe)
call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,One,Work(iCmo),nBNRe,Work(ipDysq),nBNRe,Zero,Work(ipTEMP),nBNRe)
call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,One,Work(ipTEMP),nBNRe,Work(iCmo),nBNRe,Zero,Work(ipDyM),nBNRe)
! Z
icomp = 3
call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDz),iSmLbl)
call Square(Work(ipDz),Work(ipDzsq),1,nBNRe,nBNRe)
call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,One,Work(iCmo),nBNRe,Work(ipDzsq),nBNRe,Zero,Work(ipTEMP),nBNRe)
call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,One,Work(ipTEMP),nBNRe,Work(iCmo),nBNRe,Zero,Work(ipDzM),nBNRe)
! Triangularize and reduce.
kaunt1 = 0
kaunt2 = 0
do i=1,nBNRe
  do j=1,nBNRe
    if (j <= i) then
      if (LindMOs(i) .and. LindMOs(j)) then
        Work(ipDxRe+kaunt1) = Work(ipDxM+kaunt2)
        Work(ipDyRe+kaunt1) = Work(ipDyM+kaunt2)
        Work(ipDzRe+kaunt1) = Work(ipDzM+kaunt2)
        kaunt1 = kaunt1+1
      end if
    end if
    kaunt2 = kaunt2+1
  end do
end do
! Density
DipNRe(1) = Ddot_(nBNRe**2,Work(ipDxM),1,Work(ipNRe),1)
DipNRe(2) = Ddot_(nBNRe**2,Work(ipDyM),1,Work(ipNRe),1)
DipNRe(3) = Ddot_(nBNRe**2,Work(ipDzM),1,Work(ipNRe),1)
DipRe(1) = Ddot_(nSize2,Work(ipDxRe),1,Work(ipRe),1)
DipRe(2) = Ddot_(nSize2,Work(ipDyRe),1,Work(ipRe),1)
DipRe(3) = Ddot_(nSize2,Work(ipDzRe),1,Work(ipRe),1)
Diffx = abs(DipRe(1)-DipNRe(1))
Diffy = abs(DipRe(2)-DipNRe(2))
Diffz = abs(DipRe(3)-DipNRe(3))
if (iPrint >= 5) then
  write(u6,99) iS1,iS2,'(',Diffx,',',Diffy,',',Diffz,')'
end if
! Return number
DiffMax = Diffy
if (Diffx >= Diffy) DiffMax = Diffx
if ((Diffz >= Diffx) .and. (Diffz >= Diffy)) DiffMax = Diffz
! Deallocate en masse.
call GetMem('DipX','Free','Real',ipDx,nSize1+4)
call GetMem('DipY','Free','Real',ipDy,nSize1+4)
call GetMem('DipZ','Free','Real',ipDz,nSize1+4)
call GetMem('DipXre','Free','Real',ipDxRe,nSize2)
call GetMem('DipYre','Free','Real',ipDyRe,nSize2)
call GetMem('DipZre','Free','Real',ipDzRe,nSize2)
call GetMem('DipXsq','Free','Real',ipDxsq,nBNRe**2)
call GetMem('DipYsq','Free','Real',ipDysq,nBNRe**2)
call GetMem('DipZsq','Free','Real',ipDzsq,nBNRe**2)
call GetMem('DipXm','Free','Real',ipDxM,nBNRe**2)
call GetMem('DipYm','Free','Real',ipDyM,nBNRe**2)
call GetMem('DipZm','Free','Real',ipDzM,nBNRe**2)
call GetMem('TEMP','Free','Real',ipTEMP,nBNRe**2)

return

99 format('     ',2I3,'          ',3(A,F10.7),A)

end subroutine MomentMod
