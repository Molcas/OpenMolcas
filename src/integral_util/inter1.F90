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

use Basis_Info, only: DBSC, nCnttp
use Center_Info, only: DC
use Symmetry_Info, only: nIrrep
use Molcas, only: LenIn
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
character(len=LenIn), intent(_OUT_) :: Label(*)
integer(kind=iwp), intent(_OUT_) :: Ibas_Lab(*)
real(kind=wp), intent(_OUT_) :: Coor(3,*), ZNUC(*)
integer(kind=iwp), intent(out) :: n_Cent
integer(kind=iwp) :: iCnt, iCnttp, iCo, kOp, mdc, ndc, nDiff
real(kind=wp) :: A(3)
logical(kind=iwp) :: DSCF
character(len=LenIn) :: Lbl

DSCF = .false.
nDiff = 0
call IniSew(DSCF,nDiff)

mdc = 0
ndc = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag .or. dbsc(iCnttp)%pChrg) then
    mdc = mdc+dbsc(iCnttp)%nCntr
    cycle
  end if
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    Lbl = dc(mdc)%LblCnt(1:LenIn)
    A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
    do iCo=0,nIrrep/dc(mdc)%nStab-1
      ndc = ndc+1
      kop = dc(mdc)%iCoSet(iCo,0)
      call OA(kOp,A,Coor(:,ndc))
      Label(ndc) = Lbl(1:LenIn)
      iBas_Lab(ndc) = iCnttp
      ZNUC(ndc) = real(dbsc(iCnttp)%AtmNr,kind=wp)
    end do
  end do
end do
n_cent = ndc

return

end subroutine inter1
