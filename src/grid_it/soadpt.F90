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

subroutine SOAdpt(AOValue,mAO,nCoor,mBas,nCmp,nOp,SOValue,nDeg,iAO)

use Symmetry_Info, only: nIrrep, iChTbl
use SOAO_Info, only: iAOtSO
use Basis_Info, only: MolWgh
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mAO, nCoor, mBas, nCmp, nOp, nDeg, iAO
real(kind=wp), intent(in) :: AOValue(mAO,nCoor,mBas,nCmp)
real(kind=wp), intent(inout) :: SOValue(mAO,nCoor,mBas,nCmp*nDeg)
#include "print.fh"
integer(kind=iwp) :: i1, iAux, iCmp, iPrint, iRout, iSO, j1
real(kind=wp) :: Aux(8), Fact, xa
character(len=80) :: Label

iRout = 133
iPrint = nPrint(iRout)
if (MolWgh == 0) then
  Fact = One/nDeg
else if (MolWgh == 1) then
  Fact = One
else
  Fact = One/sqrt(real(nDeg,kind=wp))
end if
iSO = 1
do i1=1,nCmp
  iaux = 0
  do j1=0,nIrrep-1
    if (iAOtSO(iAO+i1,j1) < 0) cycle
    iaux = iaux+1
    xa = iChTbl(j1,nOp)
    Aux(iAux) = Fact*xa
  end do
  if (iPrint >= 49) call RecPrt('Aux',' ',Aux,1,iAux)
  call DnaXpY(iAux,mAO*nCoor*mBas,Aux,1,AOValue(1,1,1,i1),1,0,SOValue(1,1,1,iSO),1,mAO*nCoor*mBas)
  iSO = iSO+iAux
end do

if (iPrint >= 49) then
  do iCmp=1,nCmp*nDeg
    write(Label,'(A,I2,A)') 'SOValue(mAO,nCoor,mBas,',iCmp,')'
    call RecPrt(Label,' ',SOValue(1,1,1,iCmp),mAO*nCoor,mBas)
  end do
end if

return

end subroutine SOAdpt
