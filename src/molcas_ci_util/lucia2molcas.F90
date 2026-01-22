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

subroutine LUCIA2MOLCAS(KICONF_OCC_LUCIA,KSDREO_I,NDET_LUCIA,NCSASM_LUCIA,NDTASM_LUCIA,NCNASM_LUCIA,MXPCSM,MXPORB,NCONF_PER_OPEN, &
                        NPDTCNF,NPCSCNF,nCSF_HEXS_LUCIA)
! Transfer arguments to the common blocks used by MOLCAS.

use csfbas, only: CONF, CTS, maxop_lucia, NAEL, NBEL
use splitcas_data, only: iDimBlockA
use general_data, only: ISPIN, NACTEL, NELEC3, NHOLE1, NRS1, NRS2, NSEL, NSYM, STSYM
use spinfo, only: I_ELIMINATE_GAS_MOLCAS, MINOP, MS2, NCNASM, NCNFTP, NCSASM, NCSF_HEXS, NCSFTP, NDET, NDTASM, NDTFTP, NTYP
use Molcas, only: MxSym
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: KICONF_OCC_LUCIA(*), KSDREO_I(*), NDET_LUCIA, MXPCSM, NCSASM_LUCIA(MXPCSM), NDTASM_LUCIA(MXPCSM), &
                                 NCNASM_LUCIA(MXPCSM), MXPORB, NCONF_PER_OPEN(MXPORB+1,MXPCSM), NPDTCNF(MXPORB+1), &
                                 NPCSCNF(MXPORB+1), nCSF_HEXS_LUCIA
integer(kind=iwp) :: ICL, IOPEN, ISYM, ITYP, LCONF, LDET, LLCONF, NEL2MN, NEL2MX, NORB1, NORB2

NDTASM(1:MxSym) = NDTASM_LUCIA(1:MxSym)
NCSASM(1:MxSym) = NCSASM_LUCIA(1:MxSym)
NCNASM(1:MxSym) = NCNASM_LUCIA(1:MxSym)
if (NSEL > NCSASM(STSYM)) NSEL = NCSASM(STSYM)
! For small calculations - Lasse/MGD
nCSF_HEXS = nCSF_HEXS_LUCIA
if ((I_ELIMINATE_GAS_MOLCAS > 0) .and. (NSEL > nCSF_HEXS)) NSEL = nCSF_HEXS

if (iDimBlockA > NCSASM(STSYM)) then
  write(u6,*) ''
  write(u6,*) '******************** WARNING *********************'
  write(u6,*) ' AA-Block dimension selected is bigger than the'
  write(u6,*) ' number  of CSFs reachable  within the selected'
  write(u6,*) ' Active Space. The code automatically  reset it'
  write(u6,*) ' to the number  of  CSFs. You  are  allowed  to'
  write(u6,*) ' decrease  this  number  in  your input and run'
  write(u6,*) ' again the calculation.'
  write(u6,'(1X,A,I5)') ' AA-Block dimension selected:',iDimBlockA
  iDimBlockA = NCSASM(STSYM)
  write(u6,'(1X,A,I5)') ' AA-Block dimension reset:',iDimBlockA
  write(u6,*) '**************************************************'
  write(u6,*) ''
end if

! SET INITIAL VALUES FOR LOOP COUNTERS AND ARRAY SIZES

! TOTAL NO. ORBITALS
MS2 = iSpin-1
NORB1 = sum(NRS1(1:NSYM))
NORB2 = sum(NRS2(1:NSYM))
! MIN. AND MAX. NO. OF EL IN RAS2
NEL2MX = NACTEL-2*NORB1+NHOLE1
NEL2MN = NACTEL-2*NORB1-NELEC3
NEL2MX = min(NEL2MX,2*NORB2)
if (NEL2MN < 0) NEL2MN = 0
! TOTAL NO. OF ALPHA AND BETA EL.
NAEL = (MS2+NACTEL)/2
NBEL = (NACTEL-MS2)/2
! MIN. AND MAX. NO. OF PROTOTYPE CSFS (= MIN. AND. MAX. OPEN ORBS.)
MINOP = abs(MS2)

NTYP = maxop_lucia-MINOP+1
NCNFTP(1:NTYP,1:MxSym) = NCONF_PER_OPEN(MINOP+1:maxop_lucia+1,1:MxSym)
NDTFTP(1:NTYP) = NPDTCNF(MINOP+1:maxop_lucia+1)
NCSFTP(1:NTYP) = NPCSCNF(MINOP+1:maxop_lucia+1)

! Memory needed to store ICONF array
LCONF = 0
LDET = 0
do ISYM=1,NSYM
  LLCONF = 0
  do ITYP=1,NTYP
    IOPEN = ITYP+MINOP-1
    ICL = (NACTEL-IOPEN)/2
    LLCONF = LLCONF+NCNFTP(ITYP,ISYM)*(IOPEN+ICL)
  end do
  LCONF = max(LCONF,LLCONF)
  LDET = max(LDET,NDTASM(ISYM))
end do

call mma_allocate(CONF,LCONF,label='CONF')
call mma_allocate(CTS,LDET,label='CTS')

CONF(:) = KICONF_OCC_LUCIA(1:LCONF)

CTS(:) = KSDREO_I(1:LDET)

NDET = NDET_LUCIA

end subroutine LUCIA2MOLCAS
