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

use Symmetry_Info, only: iChTbl, iOper, nIrrep, Prmt
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nComp, nIC, iDCRT, lOper(nComp), iChO(nComp)
real(kind=wp), intent(in) :: ArrIn(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Factor
real(kind=wp), intent(inout) :: ArrOut(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC)
integer(kind=iwp) :: iComp, iIC, iIrrep
real(kind=wp) :: pO, Xg
integer(kind=iwp), parameter :: iTwoj(0:7) = [1,2,4,8,16,32,64,128]
! Statement function for Cartesian index
integer(kind=iwp) :: ixyz, nElem
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!nA = (la+1)*(la+2)/2
!nB = (lb+1)*(lb+2)/2
!call RecPrt('SymAdO: ArrIn',' ',ArrIn,nZeta*nA*nB, nComp)

! Accumulate contributions

iIC = 0
do iComp=1,nComp
  pO = Prmt(iOper(iDCRT),iChO(iComp))
  do iIrrep=0,nIrrep-1
    if (iand(lOper(iComp),iTwoj(iIrrep)) == 0) cycle
    iIC = iIC+1
    Xg = real(iChTbl(iIrrep,iDCRT),kind=wp)
    call DaXpY_(nZeta*nElem(la)*nElem(lb),Xg*pO*Factor,ArrIn(1,1,1,iComp),1,ArrOut(1,1,1,iIC),1)
  end do
end do
if (iIC /= nIC) then
  call WarningMessage(2,' Abend in SymAdO: iIC /= nIC')
  write(u6,*) 'iIC,nIC=',iIC,nIC
  call Abend()
end if
!call RecPrt('SymAdO: ArrOut',' ',ArrOut,nZeta*nA*nB,nIC)

end subroutine SymAdO
