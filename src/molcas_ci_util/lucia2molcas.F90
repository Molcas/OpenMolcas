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

subroutine LUCIA2MOLCAS(KDFTP_LUCIA,KCFTP_LUCIA,KDTOC_LUCIA,KICONF_OCC_LUCIA,KSDREO_I,NDET_LUCIA,NCSASM_LUCIA,NDTASM_LUCIA, &
                        NCNASM_LUCIA,MXPCSM,MXPORB,NCONF_PER_OPEN,NPDTCNF,NPCSCNF,MULTS_LUCIA,KICTS_POINTER, &
                        nCSF_HEXS_LUCIA)
! Transfer arguments to the common blocks used by MOLCAS.

use csfbas, only: CONF, CTS, KCFTP, KDFTP, KDTOC, maxop_lucia
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MXPCSM, MXPORB, KDFTP_LUCIA, KCFTP_LUCIA, KDTOC_LUCIA, KICONF_OCC_LUCIA(*), KSDREO_I(*), &
                                 NDET_LUCIA, NCSASM_LUCIA(MXPCSM), NDTASM_LUCIA(MXPCSM), NCNASM_LUCIA(MXPCSM), &
                                 NCONF_PER_OPEN(MXPORB+1,MXPCSM), NPDTCNF(MXPORB+1), NPCSCNF(MXPORB+1), MULTS_LUCIA, nCSF_HEXS_LUCIA
integer(kind=iwp), intent(out) :: KICTS_POINTER
integer(kind=iwp) :: I, ICL, IOPEN, IORB2F, IORB2L, ISYM, ITYP, J, LCONF, LDET, LLCONF, LUCIA_TYPE, NEL1MNA, NEL1MNB, NEL2MN, NEL2MX
integer(kind=iwp), external :: ip_of_iWork
#include "rasdim.fh"
#include "ciinfo.fh"
#include "spinfo.fh"
#include "rasscf.fh"
#include "general.fh"
#include "splitcas.fh"
#include "strnum.fh"
#include "lucia_ini.fh"
#include "WrkSpc.fh"

do I=1,MXCISM
  NDTASM(I) = NDTASM_LUCIA(I)
  NCSASM(I) = NCSASM_LUCIA(I)
  NCNASM(I) = NCNASM_LUCIA(I)
end do
if (NSEL > NCSASM(STSYM)) NSEL = NCSASM(STSYM)
! For small calculations - Lasse/MGD
nCSF_HEXS = nCSF_HEXS_LUCIA
if ((I_ELIMINATE_GAS_MOLCAS > 0) .and. (NSEL > nCSF_HEXS)) then
  NSEL = nCSF_HEXS
end if

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
NORB1 = 0
NORB2 = 0
NORB3 = 0
do ISYM=1,NSYM
  NORB1 = NORB1+NRS1(ISYM)
  NORB2 = NORB2+NRS2(ISYM)
  NORB3 = NORB3+NRS3(ISYM)
end do
! FIRST AND LAST ORBITAL OF RAS1, RAS2 AND RAS3 SPACE
IORB1F = 1
IORB1L = IORB1F+NORB1-1
IORB2F = IORB1L+1
IORB2L = IORB2F+NORB2-1
IORB3F = IORB2L+1
IORB3L = IORB3F+NORB3-1
! MIN. AND MAX. NO. OF EL IN RAS1,RAS2 AND RAS3
NEL1MX = 2*NORB1
NEL1MN = NEL1MX-NHOLE1
NEL2MX = NACTEL-2*NORB1+NHOLE1
NEL2MN = NACTEL-2*NORB1-NELEC3
NEL2MX = min(NEL2MX,2*NORB2)
if (NEL2MN < 0) NEL2MN = 0
NEL3MX = NELEC3
NEL3MN = 0
! TOTAL NO. OF ALPHA AND BETA EL.
NAEL = (MS2+NACTEL)/2
NBEL = (NACTEL-MS2)/2
! MIN. AND MAX. NO. OF PROTOTYPE CSFS (= MIN. AND. MAX. OPEN ORBS.)
MINOP = abs(MS2)
maxop = maxop_lucia

NTYP = MAXOP-MINOP+1
do J=1,MXSM
  do ITYP=1,NTYP
    LUCIA_TYPE = ITYP+MINOP
    NCNFTP(ITYP,J) = NCONF_PER_OPEN(LUCIA_TYPE,J)
    NDTFTP(ITYP) = NPDTCNF(LUCIA_TYPE)
    NCSFTP(ITYP) = NPCSCNF(LUCIA_TYPE)
  end do
end do

! MIN. NO. OF EL. IN ALPHA AND BETA STRING
NL1MNA = max(0,NEL1MN-min(NBEL,NORB1))
NL1MNB = max(0,NEL1MN-min(NAEL,NORB1))
! NO. OF SINGLE EXITATIONS FOR ALPHA AND BETA STRING
NAEXCI = NAEL*(NAC-NAEL+1)
NBEXCI = NBEL*(NAC-NBEL+1)

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
KICTS_POINTER = ip_of_iWork(CTS(1))

CONF(:) = KICONF_OCC_LUCIA(1:LCONF)

CTS(:) = KSDREO_I(1:LDET)

NDET = NDET_LUCIA
MULTS = MULTS_LUCIA
MAXSYM = NSYM

NEL1MNA = max(0,NEL1MN-min(NAEL,NORB1))
NEL1MNB = max(0,NEL1MN-min(NBEL,NORB1))
NOCTPA = (NORB1-NEL1MNA+1)*(NEL3MX+1)
NOCTPB = (NORB1-NEL1MNB+1)*(NEL3MX+1)

KDFTP = KDFTP_LUCIA
KCFTP = KCFTP_LUCIA
KDTOC = KDTOC_LUCIA

return

end subroutine LUCIA2MOLCAS
