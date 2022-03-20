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

! PUT NUMBERS ON SAMPFILE.
subroutine Put9(Etot,Ract,iHowMSamp,Gmma,Gam,Esav,iDisk)

use qmstat_global, only: Cordst, iLuSaUt, iPrint, iTcSim, nCent, nPart
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
#include "WrkSpc.fh"
real(kind=wp) :: Etot, Ract, Gmma, Gam, Esav
integer(kind=iwp) :: i, iCT, iDiskHead, iDiskOld, iHowMSamp, iDisk, j
character(len=200) :: Head

iHowMSamp = iHowMSamp+1
iDiskOld = iDisk
call WrRdSim(iLuSaUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gmma,Gam,Esav) !A header
iTcSim(1) = iDisk
do i=1,3
  call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
  do j=1,nCent*nPart
    Work(iCT+j-1) = Cordst(j,i)
  end do
  call dDaFile(iLuSaUt,1,Work(iCT),nPart*nCent,iDisk)
  !The solvent coordinates.
  call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
  iTcSim(i+1) = iDisk
end do
!do i=1,3
!  call dDaFile(iLuSaUt,1,-Work(iDT(i)),nPart*nPol,iDisk)
!end do
!iTcSim(5) = iDisk
iDiskHead = iDiskOld
! Put header again, but now with a meaningful iTcSim vector that contains the table of contents which simplifies reading
call WrRdSim(iLuSaUt,1,iDiskHead,iTcSim,64,Etot,Ract,nPart,Gmma,Gam,Esav)

if (iPrint >= 15) then
  write(Head,*) ' Coordinates put on sampfile.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine PUT9
