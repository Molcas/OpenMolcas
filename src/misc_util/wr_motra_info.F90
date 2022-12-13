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

subroutine WR_MOTRA_Info(Lu,iOpt,iDisk,TCONEMO,nTCONEMO,ECOR,NSYM,NBAS,NORB,NFRO,NDEL,MxSym,BSLBL,nBSLBL)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, nTCONEMO, MxSym, nBSLBL
integer(kind=iwp), intent(inout) :: iDisk, TCONEMO(nTCONEMO), NSYM, nBas(MxSym), nOrb(MxSym), nFro(MxSym), nDel(MxSym)
real(kind=wp), intent(inout) :: ECOR
character, intent(inout) :: BSLBL(nBSLBL)
integer(kind=iwp) :: idum(1)
real(kind=wp) :: dum(1)

call iDafile(Lu,iOpt,TCONEMO,nTCONEMO,iDisk)
dum(1) = ECor
call dDafile(Lu,iOpt,dum,1,iDisk)
ECor = dum(1)
idum(1) = nSym
call iDafile(Lu,iOpt,idum,1,iDisk)
nSym = idum(1)
call iDafile(Lu,iOpt,nBas,MxSym,iDisk)
call iDafile(Lu,iOpt,nOrb,MxSym,iDisk)
call iDafile(Lu,iOpt,nFro,MxSym,iDisk)
call iDafile(Lu,iOpt,nDel,MxSym,iDisk)
call cDafile(Lu,iOpt,BSLBL,nBSLBL,iDisk)

return

end subroutine WR_MOTRA_Info
