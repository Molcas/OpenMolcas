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

subroutine ABCD_MRCI(INTSYM,indx,ISAB,C,S,ACBDS,ACBDT,BUFIN)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension INTSYM(*), indx(*), ISAB(NVIRT,NVIRT), C(*), S(*), ACBDS(*), ACBDT(*), BUFIN(*)

!vv hand-made loop unrolling to fix a bug in GCC 3.x
IAD16 = 0
INSIN = KBUFF1
call CSCALE(indx,INTSYM,C,SQ2)
call CSCALE(indx,INTSYM,S,SQ2INV)
NVT = IROW(NVIRT+1)
NOV = (NVT-1)/IPASS+1
IACMAX = 0
!do ISTEP=1,IPASS
ISTEP = 1
if (IPASS < 1) goto 670

770 IACMIN = IACMAX+1
IACMAX = IACMAX+NOV
if (IACMAX > NVT) IACMAX = NVT
if (IACMIN > IACMAX) GO TO 70
!do ISYM=1,NSYM
ISYM = 1
if (NSYM < 1) goto 640
740 IST1 = IRC(3)+JJS(ISYM+9)+1
IFIN1 = IRC(3)+JJS(ISYM+10)
INPS = IFIN1-IST1+1
IST2 = IRC(2)+JJS(ISYM)+1
IFIN2 = IRC(2)+JJS(ISYM+1)
INPT = IFIN2-IST2+1
ITAIL = INPS+INPT
if (ITAIL == 0) GO TO 40
IN1 = -NVIRT
!do NA=1,NVIRT
NA = 1
if (NVIRT < 1) goto 650
750 IN1 = IN1+NVIRT
!do NC=1,NA
NC = 1
if (NA < 1) goto 660

760 IAC = IROW(NA)+NC
if (IAC < IACMIN) GO TO 60
if (IAC > IACMAX) GO TO 60
if (NA == 1) GO TO 60
NSAC = MUL(NSM(LN+NA),NSM(LN+NC))
NSACL = MUL(NSAC,LSYM)
if (NSACL /= ISYM) GO TO 60
ISAC = ISAB(NA,NC)
NSC = NSM(LN+NC)
NDMAX = NVIRP(NSC)+NVIR(NSC)
if (NDMAX > NA) NDMAX = NA
INS = ISAB(NA,NDMAX)
! MOVE INS ITEMS FROM FILE, UNIT 16, VIA BUFFER, INTO ACBDS,
! AND THEN INTO ACBDT:
ILOOP = 0
72 INSB = INS
73 if (INSIN < KBUFF1) GO TO 75
! INSB ITEMS REMAIN TO MOVE.
! INSIN ITEMS HAVE ALREADY BEEN MOVED FROM THE BUFFER.
call dDAFILE(Lu_80,2,BUFIN,KBUFF1,IAD16)
INSIN = 0
75 INB = KBUFF1-INSIN
! INB FRESH ITEMS ARE STILL REMAINING IN BUFFER.
INUMB = min(INSB,INB)
! MOVE INUMB ITEMS.
IST = INS-INSB+1
if (ILOOP == 0) call DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDS(IST),1)
if (ILOOP == 1) call DCOPY_(INUMB,BUFIN(INSIN+1),1,ACBDT(IST),1)
INSIN = INSIN+INUMB
INSB = INSB-INUMB
if (INSB > 0) GO TO 73
ILOOP = ILOOP+1
if (ILOOP == 1) GO TO 72
! INS ITEMS HAVE BEEN TRANSFERRED TO ACBDS AND TO ACBDT.
if (INPS == 0) GO TO 11
INDA = IST1
if (IFIN1 < IST1) goto 610
!do INDA=IST1,IFIN1
710 TERM = DDOT_(INS,C(indx(INDA)+1),1,ACBDS,1)
S(indx(INDA)+ISAC) = S(indx(INDA)+ISAC)+TERM
call DAXPY_(INS,C(indx(INDA)+ISAC),ACBDS,1,S(indx(INDA)+1),1)
!end do
INDA = INDA+1
if (INDA <= IFIN1) goto 710
610 continue
11 if ((INPT == 0) .or. (NA == NC)) GO TO 60
INDA = IST2
if (IFIN2 < IST2) goto 630
!do INDA=IST2,IFIN2
730 TERM = DDOT_(INS,C(indx(INDA)+1),1,ACBDT,1)
S(indx(INDA)+ISAC) = S(indx(INDA)+ISAC)+TERM
call DAXPY_(INS,C(indx(INDA)+ISAC),ACBDT,1,S(indx(INDA)+1),1)
!end do
INDA = INDA+1
if (INDA <= IFIN2) goto 730
630 continue
60 continue
!end do
NC = NC+1
if (NC <= NA) goto 760
660 continue
!vv end of unrolling loop
!NC = NC+1
!if (NC /= NA) goto 61
!end do
NA = NA+1
if (NA <= NVIRT) goto 750
650 continue
40 continue
!end do
ISYM = ISYM+1
if (ISYM <= NSYM) goto 740
640 continue

70 continue
!end do
ISTEP = ISTEP+1
if (ISTEP <= IPASS) goto 770
670 continue
call CSCALE(indx,INTSYM,C,SQ2INV)
call CSCALE(indx,INTSYM,S,SQ2)

return

end subroutine ABCD_MRCI
