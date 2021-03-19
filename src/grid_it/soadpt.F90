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

implicit real*8(a-h,o-z)
#include "print.fh"
#include "real.fh"
real*8 AOValue(mAO,nCoor,mBas,nCmp), SOValue(mAO,nCoor,mBas,nCmp*nDeg), Aux(8)
character*80 Label

iRout = 133
iPrint = nPrint(iRout)
!call GetMem('SOAdpt_E','CHEC','REAL',iDum,iDum)

if (MolWgh == 0) then
  Fact = One/dble(nDeg)
else if (MolWgh == 1) then
  Fact = One
else
  Fact = One/sqrt(dble(nDeg))
end if
iSO = 1
do i1=1,nCmp
  iaux = 0
  do j1=0,nIrrep-1
    if (iAOtSO(iAO+i1,j1) < 0) cycle
    iaux = iaux+1
    xa = dble(iChTbl(j1,nOp))
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

!call GetMem('SOAdpt_X','CHEC','REAL',iDum,iDum)

return

end subroutine SOAdpt
