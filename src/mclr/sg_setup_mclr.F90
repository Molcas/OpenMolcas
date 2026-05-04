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

use sguga, only: SGS, CIS, EXS, MkCOT, MkCList, MkSGNum, SG_Init_Simple
use input_mclr, only: iSpin, nActEl, nConf, nCSF, nElec3, nHole1, nRS1, nRS2, nRS3, nSym, State_Sym
use Definitions, only: iwp

implicit none
integer(kind=iwp) pState_Sym

Call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,EXS,nHole1,nElec3,nRs1,nRs2,nRs3)

! PURPOSE: FREE THE GUGA TABLES
! FORM VARIOUS OFFSET TABLES:
! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

call MKCOT(SGS,CIS)

! CONSTRUCT THE CASE LIST

call MKCLIST(SGS,CIS)

! SET UP ENUMERATION TABLES

call MKSGNUM(pState_Sym,SGS,CIS,EXS)

end subroutine SG_Setup_MCLR
