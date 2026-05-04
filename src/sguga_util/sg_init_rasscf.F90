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
!               Roland Lindh                                           *
!***********************************************************************

subroutine SG_Init_RASSCF(nSym,nActEl,iSpin,                         &
                          SGS,CIS,EXS,                               &
                          nHole1,nElec3,nRs1,nRs2,nRs3,              &
                          STSYM,DoblockDMRG)
! PURPOSE: CONTROL ROUTINE TO SET UP GUGA TABLES
! AUTHOR:  P.-AA. MALMQVIST
!
! MODIFIED TO FIT THE DETRAS PROGRAM BY M.P. FUELSCHER

use sguga, only: CIStruct, EXStruct, SGStruct, MkCOT, MkCList, MkSGNum, SG_Init_Simple
use Definitions, only: iwp
use gas_data, only: NGAS, NGSSH
use rasscf_global, only: NSM

implicit none
integer(kind=iwp), intent(in) :: nSym, iSpin, nActEl, nHole1, nElec3, nRs1(nSym), nRs2(nSym), nRs3(nSym), STSYM
type(SGStruct), intent(out) :: SGS
type(CIStruct), intent(out) :: CIS
type(EXStruct), intent(out) :: EXS
logical(kind=iwp), intent(in) :: DoBlockDMRG

integer(kind=iwp) :: IGAS, ISYM, NLEV, NSTA

NLEV = 0
do IGAS=1,NGAS
  do ISYM=1,NSYM
    NSTA = NLEV+1
    NLEV = NLEV+NGSSH(IGAS,ISYM)
    NSM(NSTA:NLEV) = ISYM
  end do
end do

Call SG_Init_Simple(nSym,nActEl,iSpin,                         &
                    SGS,CIS,EXS,                               &
                    nHole1,nElec3,nRs1,nRs2,nRs3,              &
                    xNLEV=NLEV,xNSM=NSM)

if (SGS%NVERT0 == 0) then
  CIS%NCSF(STSYM) = 0
  return
end if

if (doBlockDMRG) then
  CIS%NCSF(STSYM) = 1
  return
end if

! FORM VARIOUS OFFSET TABLES:
! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
!       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

call MKCOT(SGS,CIS)

! CONSTRUCT THE CASE LIST

call MKCLIST(SGS,CIS)

! SET UP ENUMERATION TABLES

call MKSGNUM(STSYM,SGS,CIS,EXS)

if (NActEl == 0) CIS%NCSF(STSYM) = 1

end subroutine SG_Init_RASSCF
