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

subroutine SymAdO(ArrIn,nZeta,la,lb,nComp,ArrOut,nIC,iDCRT,lOper,iChO,Factor)

use Index_Functions, only: nTri_Elem1
use Symmetry_Info, only: iChTbl, iOper, nIrrep, Prmt
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nComp, nIC, iDCRT, lOper(nComp), iChO(nComp)
real(kind=wp), intent(in) :: ArrIn(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp), Factor
real(kind=wp), intent(inout) :: ArrOut(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nIC)
integer(kind=iwp) :: iComp, iIC, iIrrep
real(kind=wp) :: pO, Xg

!nA = nTri_Elem1(la)
!nB = nTri_Elem1(lb)
!call RecPrt('SymAdO: ArrIn',' ',ArrIn,nZeta*nA*nB, nComp)

! Accumulate contributions

iIC = 0
do iComp=1,nComp
  pO = Prmt(iOper(iDCRT),iChO(iComp))
  do iIrrep=0,nIrrep-1
    if (.not. btest(lOper(iComp),iIrrep)) cycle
    iIC = iIC+1
    Xg = real(iChTbl(iIrrep,iDCRT),kind=wp)
    ArrOut(:,:,:,iIC) = ArrOut(:,:,:,iIC)+Xg*pO*Factor*ArrIn(:,:,:,iComp)
  end do
end do
if (iIC /= nIC) then
  call WarningMessage(2,' Abend in SymAdO: iIC /= nIC')
  write(u6,*) 'iIC,nIC=',iIC,nIC
  call Abend()
end if
!call RecPrt('SymAdO: ArrOut',' ',ArrOut,nZeta*nA*nB,nIC)

end subroutine SymAdO
