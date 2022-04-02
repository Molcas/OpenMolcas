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
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_IN_) :: Etot, Ract, Gmma, Gam, Esav
integer(kind=iwp), intent(inout) :: iHowMSamp, iDisk
integer(kind=iwp) :: i, iDiskHead, iDiskOld
character(len=200) :: Head
real(kind=wp), allocatable :: CT(:)

iHowMSamp = iHowMSamp+1
iDiskOld = iDisk
call WrRdSim(iLuSaUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gmma,Gam,Esav) !A header
iTcSim(1) = iDisk
call mma_allocate(CT,nPart*nCent)
do i=1,3
  CT(:) = Cordst(i,1:nCent*nPart)
  call dDaFile(iLuSaUt,1,CT,nPart*nCent,iDisk)
  ! The solvent coordinates.
  iTcSim(i+1) = iDisk
end do
call mma_deallocate(CT)
!call dDaFile(iLuSaUt,1,-DT,3*nPart*nPol,iDisk)
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
