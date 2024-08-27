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

subroutine inter1(Label,iBas_Lab,Coor,ZNUC,N_Cent)

use Basis_Info, only: nCnttp, DBSC
use Center_Info, only: DC
use Symmetry_Info, only: nIrrep

implicit none
#include "Molcas.fh"
character(len=LENIN) Label(*)
integer Ibas_Lab(*)
real*8 Coor(3,*), ZNUC(*)
integer n_Cent
real*8 A(3)
character(len=LENIN) Lbl
logical DSCF
integer nDiff, mdc, ndc, iCnttp, iCnt, iCo, kOp

DSCF = .false.
nDiff = 0
call IniSew(DSCF,nDiff)

mdc = 0
ndc = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%pChrg) then
    mdc = mdc+dbsc(iCnttp)%nCntr
    Go To 99
  end if
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    Lbl = dc(mdc)%LblCnt(1:LENIN)
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
    do iCo=0,nIrrep/dc(mdc)%nStab-1
      ndc = ndc+1
      kop = dc(mdc)%iCoSet(iCo,0)
      call OA(kOp,A,Coor(1:3,ndc))
      Label(ndc) = Lbl(1:LENIN)
      iBas_Lab(ndc) = iCnttp
      ZNUC(ndc) = dble(dbsc(iCnttp)%AtmNr)
    end do
  end do
99 continue
end do
n_cent = ndc

return

end subroutine inter1
