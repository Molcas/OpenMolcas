!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************

subroutine mk_MOs(SOValue,mAO,nCoor,MOValue,nMOs,CMOs,nCMO)

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: mAO, nCoor, nMOs, nCMO
real(kind=wp), intent(in) :: SOValue(mAO*nCoor,nMOs), CMOs(nCMO)
real(kind=wp), intent(out) :: MOValue(mAO*nCoor,nMOs)
integer(kind=iwp) :: iCMO, iIrrep, iSO
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iOff
character(len=80) :: Label
#endif

#ifdef _DEBUGPRINT_
write(u6,*) 'mk_MOs: MO-Coefficients'
iOff = 1
do iIrrep=0,nIrrep-1
  if (nBas(iIrrep) > 0) then
    write(u6,*) ' Symmetry Block',iIrrep
    call RecPrt(' ',' ',CMOs(iOff),nBas(iIrrep),nBas(iIrrep))
  end if
  iOff = iOff+nBas(iIrrep)**2
end do
#endif

! Compute some offsets

iSO = 1
iCMO = 1
do iIrrep=0,nIrrep-1
  if (nBas(iIrrep) == 0) cycle
  call DGeMM_('N','N',mAO*nCoor,nBas(iIrrep),nBas(iIrrep),One,SOValue(:,iSO:),mAO*nCoor,CMOs(iCMO:),nBas(iIrrep),Zero, &
              MOValue(:,iSO:),mAO*nCoor)
  iSO = iSO+nBas(iIrrep)
  iCMO = iCMO+nBas(iIrrep)*nBas(iIrrep)
end do

#ifdef _DEBUGPRINT_
write(Label,'(A)') 'mk_MOs: MOValue(mAO*nCoor,nMOs)'
call RecPrt(Label,' ',MOValue(1,1),mAO*nCoor,nMOs)
#endif

return

end subroutine mk_MOs
