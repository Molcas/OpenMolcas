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

use qmstat_global, only: Cordst, iLuStut, iPrint, iTcSim, nCent, nPart, StFilUt
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_IN_) :: Ract, Etot, Gmma, Gam, Esav
integer(kind=iwp) :: i, iDisk
character(len=200) :: Head
real(kind=wp), allocatable :: CT(:)

call DaName(iLuStUt,StFilUt) !Here follows a general output to the startfile
iDisk = 0
call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gmma,Gam,Esav)
iTcSim(1) = iDisk
! In this loop the coordinates are put on file.
! The loop is needed due to how Cordst is statically allocated.
call mma_allocate(CT,nPart*nCent,label='CTemp')
do i=1,3
  CT(:) = Cordst(i,1:nPart*nCent)
  call dDaFile(iLuStUt,1,CT,nPart*nCent,iDisk)
  iTcSim(1+i) = iDisk
end do
call mma_deallocate(CT)
iDisk = 0
call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gmma,Gam,Esav)
call DaClos(iLuStUt)
if (iPrint >= 10) then !Print the stored configuration.
  write(Head,*) ' Coordinates put on the startfile solvent configuration.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine PUT8
