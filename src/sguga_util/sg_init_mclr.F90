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
subroutine SG_Init_MCLR(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,ksym)

use sguga, only: CIStruct, EXStruct, SGStruct, MkCOT, MkCList, MkSGNum, SG_Init_Simple
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, iSpin, nActEl, nHole1, nElec3, nRs1(nSym), nRs2(nSym), nRs3(nSym), ksym
type(SGStruct), intent(out) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS

Call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,EXS,nHole1,nElec3,nRs1,nRs2,nRs3)

! PURPOSE: FREE THE GUGA TABLES
! FORM VARIOUS OFFSET TABLES:
! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

call MKCOT(SGS,CIS)

! CONSTRUCT THE CASE LIST

call MKCLIST(SGS,CIS)

! SET UP ENUMERATION TABLES

call MKSGNUM(kSYM,SGS,CIS,EXS)

end subroutine SG_Init_MCLR
