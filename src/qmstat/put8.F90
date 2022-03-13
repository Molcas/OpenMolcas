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

! PUT NUMBERS ON STARTFILE.
subroutine Put8(Ract,Etot,Gmma,Gam,Esav)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Ract, Etot, Gmma, Gam, Esav
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iCT, iDisk, j
character(len=200) :: Head

call DaName(iLuStUt,StFilUt) !Here follows a general output to the startfile
iDisk = 0
call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gmma,Gam,Esav)
iTcSim(1) = iDisk
! In this loop the coordinates are put on file.
! The loop is needed due to how Cordst is statically allocated.
do i=1,3
  call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
  do j=1,nPart*nCent
    Work(iCT+j-1) = Cordst(j,i)
  end do
  call dDaFile(iLuStUt,1,Work(iCT),nPart*nCent,iDisk)
  iTcSim(1+i) = iDisk
  call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
end do
iDisk = 0
call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gmma,Gam,Esav)
call DaClos(iLuStUt)
if (iPrint >= 10) then !Print the stored configuration.
  write(Head,*) ' Coordinates put on the startfile solvent configuration.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine PUT8
