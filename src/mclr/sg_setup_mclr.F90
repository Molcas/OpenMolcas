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

subroutine SG_Setup_MCLR(pState_Sym)

use molcas, only: MxLev
use sguga, only: SGS, CIS, EXS, MkCOT, MkSGNum, SG_Init_Simple
use input_mclr, only: iSpin, nActEl, nElec3, nHole1, nRS1, nRS2, nRS3, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in):: pState_Sym
integer(kind=iwp) :: iBas, nLev, iSym, ISM(1:MxLev), Level(MxLev), iq

nLev = 0
do iSym=1,nSym
  do iBas=1,nRs1(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do
do iSym=1,nSym
  do iBas=1,nRs2(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do
do iSym=1,nSym
  do iBas=1,nRs3(iSym)
    nLev = nLev+1
    ISM(nLev) = iSym
  end do
end do

Level(1:MxLev)=[(iq,iq=1,MxLev)]

Call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,EXS,nHole1,nElec3,nRs1,nRs2,nRs3, &
                    xLevel=Level, xL2Act=Level,                                 &
                    xNLEV=nLev, xNSM=ISM)

! PURPOSE: FREE THE GUGA TABLES
! FORM VARIOUS OFFSET TABLES:
! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

! CONSTRUCT THE CASE LIST
call MKCOT(SGS,CIS)

! SET UP ENUMERATION TABLES

call MKSGNUM(pState_Sym,SGS,CIS,EXS)

end subroutine SG_Setup_MCLR
