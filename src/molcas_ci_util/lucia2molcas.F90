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

subroutine LUCIA2MOLCAS(NDET_LUCIA,NCSASM_LUCIA,NDTASM_LUCIA,MXPCSM,MXPORB,NCONF_PER_OPEN,NPDTCNF,NPCSCNF, &
                        nCSF_HEXS_LUCIA)
! Transfer arguments to the common blocks used by MOLCAS.

use csfbas, only: maxop_lucia, NAEL, NBEL
use CI_Input, only: ISPIN, NACTEL, STSYM
use rasscf_global, only: NSEL, NDET
use CI_Info, only: CI_St, Space_Inp
use CI_Space_Def, only: Space_Removals
use lucia_data, only: MINOP, NCNFTP, NCSFTP, NDTFTP, NTYP
use Molcas, only: MxSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NDET_LUCIA, MXPCSM, NCSASM_LUCIA(MXPCSM), NDTASM_LUCIA(MXPCSM), MXPORB, &
                                 NCONF_PER_OPEN(MXPORB+1,MXPCSM), NPDTCNF(MXPORB+1), NPCSCNF(MXPORB+1), nCSF_HEXS_LUCIA
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp) :: MS2

! the CI dimensions of the RASSCF state (CI_Info)
CI_St(iState)%nDet(1:MxSym) = NDTASM_LUCIA(1:MxSym)
CI_St(iState)%nCSF(1:MxSym) = NCSASM_LUCIA(1:MxSym)
if (NSEL > CI_St(iState)%nCSF(STSYM)) NSEL = CI_St(iState)%nCSF(STSYM)
! For small calculations - Lasse/MGD
CI_St(iState)%nCSF_Kept = nCSF_HEXS_LUCIA
if (Space_Removals(Space_Inp) .and. (NSEL > CI_St(iState)%nCSF_Kept)) NSEL = CI_St(iState)%nCSF_Kept

MS2 = iSpin-1
! TOTAL NO. OF ALPHA AND BETA EL.
NAEL = (MS2+NACTEL)/2
NBEL = (NACTEL-MS2)/2
! MIN. AND MAX. NO. OF PROTOTYPE CSFS (= MIN. AND. MAX. OPEN ORBS.);
! MINOP = abs(MS2) is Lucia's own (csfdim_gas)

NTYP = maxop_lucia-MINOP+1
NCNFTP(1:NTYP,1:MxSym) = NCONF_PER_OPEN(MINOP+1:maxop_lucia+1,1:MxSym)
NDTFTP(1:NTYP) = NPDTCNF(MINOP+1:maxop_lucia+1)
NCSFTP(1:NTYP) = NPCSCNF(MINOP+1:maxop_lucia+1)

NDET = NDET_LUCIA

end subroutine LUCIA2MOLCAS
