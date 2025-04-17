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
subroutine GugaNew(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,CIL,imode,ksym,State_Sym)

use Str_Info, only: CFTP, CNSM
use gugx, only: CIStruct, EXStruct, SGStruct
use MkGUGA_mod, only: MKGUGA
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nSym, iSpin, nActEl, nHole1, nElec3, nRs1(nSym), nRs2(nSym), nRs3(nSym), imode, ksym, State_Sym
type(SGStruct), intent(out) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
real(kind=wp), intent(inout) :: CIL(*)
integer(kind=iwp) :: iS, iss, NCONF, nRas1T, nRas2T, nRas3T
#ifdef _DEBUGPRINT_
real(kind=wp), parameter :: PRWTHR = 0.05_wp
#endif

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

call mkGUGA(SGS,CIS)

! PURPOSE: FREE THE GUGA TABLES
! FORM VARIOUS OFFSET TABLES:
! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

call MKCOT(SGS,CIS)

! CONSTRUCT THE CASE LIST

call MKCLIST(SGS,CIS)

! SET UP ENUMERATION TABLES

call MKSGNUM(kSYM,SGS,CIS,EXS)

nConf = CIS%nCSF(kSym)

iss = 1
if (ksym /= state_sym) iss = 2

#ifdef _DEBUGPRINT_
write(u6,101)
write(u6,102) PRWTHR
call SGPRWF(SGS,CIS,ksym,PRWTHR,SGS%iSpin,CIL,nConf,.false.,-99)
write(u6,103)
#endif

call REORD(SGS,CIS,EXS,NCONF,iMode,CNSM(iss)%ICONF,CFTP,kSym,CIL)

#ifdef _DEBUGPRINT_
call SGPRWF(SGS,CIS,ksym,PRWTHR,SGS%iSpin,CIL,nConf,.false.,-99)
101 format(/,6X,100('-'),/,6X,29X,'Wave function printout: Split Graph format',/, &
           6X,8X,'in parenthesis: midvertex, upper-walk symmetry upper- and lower-walk serial numbers',/,6X,100('-'),/)
102 format(6X,'printout of CI-coefficients larger than',F6.2)
103 format(/,6X,100('-'),/)
#endif

end subroutine GugaNew
