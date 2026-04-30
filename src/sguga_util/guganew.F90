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
subroutine GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,ksym)

use sguga, only: CIStruct, EXStruct, SGStruct, MkSGUGA, MkCOT
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, iSpin, nActEl, nHole1, nElec3, nRs1(nSym), nRs2(nSym), nRs3(nSym), ksym
type(SGStruct), intent(out) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: iS, nRas1T, nRas2T, nRas3T

nRas1T = sum(nRs1(1:nSym))
nRas2T = sum(nRs2(1:nSym))
nRas3T = sum(nRs3(1:nSym))

SGS%nSym = nSym
SGS%iSpin = iSpin
SGS%nActEl = nActEl

! COMPUTE RAS RESTRICTIONS ON VERTICES:

SGS%LV1RAS = NRAS1T
SGS%LV3RAS = nRas1T+NRAS2T
SGS%LM1RAS = 2*nRas1T-NHOLE1
SGS%LM3RAS = NACTEL-nElec3

! SET IFRAS FLAG
! IFRAS = 0 : THIS IS A CAS CALCULATION
! IFRAS = 1 : THIS IS A RAS CALCULATION

if ((NRAS1T+NRAS3T) /= 0) then
  SGS%IFRAS = 1
else
  SGS%IFRAS = 0
end if
do IS=1,NSYM
  if ((SGS%IFRAS /= 0) .and. (nRs2(IS) /= 0)) SGS%IFRAS = SGS%IFRAS+1
end do

Call MkSGUGA(SGS,CIS)

! PURPOSE: FREE THE GUGA TABLES
! FORM VARIOUS OFFSET TABLES:
! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

call MKCOT(SGS,CIS)

! CONSTRUCT THE CASE LIST

call MKCLIST(SGS,CIS)

! SET UP ENUMERATION TABLES

call MKSGNUM(kSYM,SGS,CIS,EXS)

end subroutine GugaNew
