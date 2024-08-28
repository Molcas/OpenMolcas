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

!#define _DEBUGPRINT_
subroutine SODist(SOValue,mAO,nCoor,mBas,nCmp,nDeg,MOValue,nMOs,iAO,CMOs,nCMO,DoIt)

use SOAO_Info, only: iAOtSO
use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Constants
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer mAO, nCoor, mBas, nCmp, nDeg, nCMO, nMOs
real*8 SOValue(mAO*nCoor,mBas,nCmp*nDeg), MOValue(mAO*nCoor,nMOs), CMOs(nCMO)
integer DoIt(nMOs)
integer iOff_MO(0:7), iOff_CMO(0:7)
integer iIrrep, itmp1, itmp2, i1, iDeg, iSO, iOff, iMO, iCMO, iAO
#ifdef _DEBUGPRINT_
character(len=80) Label
#endif

#ifdef _DEBUGPRINT_
write(6,*) 'SODist: MO-Coefficients'
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

itmp1 = 1
itmp2 = 0
do iIrrep=0,nIrrep-1
  iOff_MO(iIrrep) = itmp1
  iOff_CMO(iIrrep) = itmp2
  itmp1 = itmp1+nBas(iIrrep)
  itmp2 = itmp2+nBas(iIrrep)*nBas(iIrrep)
end do

do i1=1,nCmp
  iDeg = 0
  do iIrrep=0,nIrrep-1
    iSO = iAOtSO(iAO+i1,iIrrep)
    if (iSO < 0) cycle
    iDeg = iDeg+1
    iOff = (i1-1)*nDeg+iDeg

    ! Distribute contribution to all MO's in this irrep

    iMO = iOff_MO(iIrrep)
    iCMO = iOff_CMO(iIrrep)+iSO
    call MyDGeMM(DoIt(iMO),mAO*nCoor,nBas(iIrrep),mBas,SOValue(1,1,iOff),mAO*nCoor,CMOs(iCMO),nBas(iIrrep),MOValue(1,iMO),mAO*nCoor)
  end do
end do

#ifdef _DEBUGPRINT_
write(Label,'(A)') 'SODist: MOValue(mAO*nCoor,nMOs)'
call RecPrt(Label,' ',MOValue(1,1),mAO*nCoor,nMOs)
#endif

end subroutine SODist
