!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Per Ake Malmqvist                                      *
!               Markus P. Fuelscher                                    *
!***********************************************************************

subroutine GUGACTL(nSym,iSpin,nActEl,nHole1,nElec3,nRs1,nRs2,nRs3,SGS,CIS,EXS,STSYM,DoblockDMRG)
! PURPOSE: CONTROL ROUTINE TO SET UP GUGA TABLES
! AUTHOR:  P.-AA. MALMQVIST
!
! MODIFIED TO FIT THE DETRAS PROGRAM BY M.P. FUELSCHER

use gugx, only: CIStruct, EXStruct, SGStruct
use MkGUGA_mod, only: MKGUGA
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nSym, iSpin, nActEl, nHole1, nElec3, nRs1(nSym), nRs2(nSym), nRs3(nSym), STSYM
type(SGStruct), intent(out) :: SGS
type(CIStruct), intent(out) :: CIS
type(EXStruct), intent(out) :: EXS
logical(kind=iwp), intent(in) :: DoBlockDMRG
integer(kind=iwp) :: IS, nRas1T, nRas2T, nRas3T

nRas1T = sum(nRs1(1:nSym))
nRas2T = sum(nRs2(1:nSym))
nRas3T = sum(nRs3(1:nSym))

#ifdef _DEBUGPRINT_
write(u6,*) 'nSym,iSpin,nActEl,nHole1,nElec3,nRas1T,nRas2T,nRas3T,STSYM=',nSym,iSpin,nActEl,nHole1,nElec3,nRas1T,nRas2T,nRas3T,STSYM
#endif

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

if (NRAS1T+NRAS3T /= 0) then
  SGS%IFRAS = 1
  do IS=1,NSYM
    if (nRs1(IS)+nRs2(IS)+nRs3(IS) /= 0) SGS%IFRAS = SGS%IFRAS+1
  end do
else
  SGS%IFRAS = 0
end if

! INITIALIZE GUGA TABLES:

call MKGUGA(SGS,CIS)

if (SGS%NVERT0 == 0) then
  CIS%NCSF(STSYM) = 0
  return
end if
if (doBlockDMRG) then
  CIS%NCSF(STSYM) = 1
  return
end if

! PURPOSE: FREE THE GUGA TABLES
! FORM VARIOUS OFFSET TABLES:
! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

call MKCOT(SGS,CIS)

! CONSTRUCT THE CASE LIST

call MKCLIST(SGS,CIS)

! SET UP ENUMERATION TABLES

call MKSGNUM(STSYM,SGS,CIS,EXS)

if (NActEl == 0) CIS%NCSF(STSYM) = 1

! (IFRAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
! IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
! SET UP THE REINDEXING TABLE

call SETSXCI()

end subroutine GUGACTL
