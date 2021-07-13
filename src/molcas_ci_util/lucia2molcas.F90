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
                        NCNASM_LUCIA,MXPCSM,MXPORB,NCONF_PER_OPEN,NPDTCNF,NPCSCNF,MULTS_LUCIA,NSSOA,NSSOB,KICTS_POINTER, &
                        nCSF_HEXS_LUCIA)
! Transfer arguments to the common blocks used by MOLCAS.

implicit real*8(A-H,O-Z)
#include "rasdim.fh"
#include "csfbas.fh"
#include "ciinfo.fh"
#include "spinfo.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "splitcas.fh"
#include "strnum.fh"
#include "detbas.fh"
#include "lucia_ini.fh"
#include "WrkSpc.fh"
dimension NCSASM_LUCIA(MXPCSM), KICONF_OCC_LUCIA(*)
dimension KSDREO_I(*)
dimension NDTASM_LUCIA(MXPCSM), NCNASM_LUCIA(MXPCSM)
dimension NCONF_PER_OPEN(MXPORB+1,MXPCSM)
dimension NPDTCNF(MXPORB+1), NPCSCNF(MXPORB+1)
dimension NSSOA(*), NSSOB(*)

do I=1,MXCISM
  NDTASM(I) = NDTASM_LUCIA(I)
  NCSASM(I) = NCSASM_LUCIA(I)
  NCNASM(I) = NCNASM_LUCIA(I)
end do
if (NSEL > NCSASM(STSYM)) NSEL = NCSASM(STSYM)
! For small calculations - Lasse/MGD
nCSF_HEXS = nCSF_HEXS_LUCIA
if ((N_ELIMINATED_GAS_MOLCAS > 0) .and. (NSEL > nCSF_HEXS)) then
  NSEL = nCSF_HEXS
end if

if (iDimBlockA > NCSASM(STSYM)) then
  write(6,*) ''
  write(6,*) '******************** WARNING *********************'
  write(6,*) ' AA-Block dimension selected is bigger than the '
  write(6,*) ' number  of CSFs reachable  within the selected '
  write(6,*) ' Active Space. The code automatically  reset it '
  write(6,*) ' to the number  of  CSFs. You  are  allowed  to '
  write(6,*) ' decrease  this  number  in  your input and run '
  write(6,*) ' again the calculation.'
  write(6,'(1X,A,I5)') ' AA-Block dimension selected:',iDimBlockA
  iDimBlockA = NCSASM(STSYM)
  write(6,'(1X,A,I5)') ' AA-Block dimension reset:',iDimBlockA
  write(6,*) '**************************************************'
  write(6,*) ''
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
  LDET = max(LDET,NDTASM(ISYM))
  do ITYP=1,NTYP
    IOPEN = ITYP+MINOP-1
    ICL = (NACTEL-IOPEN)/2
    LLCONF = LLCONF+NCNFTP(ITYP,ISYM)*(IOPEN+ICL)
  end do
  LCONF = max(LCONF,LLCONF)
end do

call GetMem('KICONF','Allo','Integer',KICONF(1),LCONF)
call GetMem('KICTS','Allo','Integer',KICTS(1),LDET)
KICTS_POINTER = KICTS(1)

do I=1,LCONF
  IWORK(KICONF(1)+I-1) = KICONF_OCC_LUCIA(I)
end do

do I=1,LDET
  IWORK(KICTS(1)+I-1) = KSDREO_I(I)
end do

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
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(NSSOA)
  call Unused_integer_array(NSSOB)
end if

end subroutine LUCIA2MOLCAS

subroutine LUCIA2MOLCAS_FREE

implicit real*8(A-H,O-Z)
#include "rasdim.fh"
#include "csfbas.fh"
#include "ciinfo.fh"
#include "spinfo.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "splitcas.fh"
#include "strnum.fh"
#include "detbas.fh"
#include "lucia_ini.fh"
#include "WrkSpc.fh"

! Memory needed to store ICONF array
LCONF = 0
LDET = 0
do ISYM=1,NSYM
  LLCONF = 0
  LDET = max(LDET,NDTASM(ISYM))
  do ITYP=1,NTYP
    IOPEN = ITYP+MINOP-1
    ICL = (NACTEL-IOPEN)/2
    LLCONF = LLCONF+NCNFTP(ITYP,ISYM)*(IOPEN+ICL)
  end do
  LCONF = max(LCONF,LLCONF)
end do

call GetMem('KICONF','Free','Integer',KICONF(1),LCONF)
call GetMem('KICTS','Free','Integer',KICTS(1),LDET)

end subroutine LUCIA2MOLCAS_FREE
