************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      INTEGER FUNCTION NEWFSBTAB(NACTEL,MSPIN2,LSYM,REST,SSTAB)
      use stdalloc, only: mma_deallocate
      use rassi_global_arrays, only: FSBARR
      IMPLICIT NONE
      INTEGER NACTEL,MSPIN2,LSYM
      INTEGER SSTAB(*), REST(*)

      INTEGER LFSBTAB,LSSTARR,NSIZE,ITYPE
      INTEGER NASPRT
      INTEGER NSSTARR,NPART,NSYM
      INTEGER NRDETS,NRDETS0,NFSB,NFSB0
      INTEGER KORB,KREST,IFSB,IERR
      INTEGER NHEAD,NHSHMAP,KHSHMAP,LHSHMAP,NULL,JFSB
#include "WrkSpc.fh"
C Purpose: Construct an FSB table and return its address in the
C IWORK array.
C ITYPE=73 is the check code for this table.
      ITYPE=73
      IF(SSTAB(2).NE.19) THEN
        WRITE(6,*)' NEWFSBTAB error: Not a Substring Table.'
        CALL ABEND()
      END IF
      IF(REST(2).NE.91) THEN
        WRITE(6,*)' NEWFSBTAB error: Not a GAS Restriction Table.'
        CALL ABEND()
      END IF
      NSYM  =SSTAB(4)
      NASPRT=SSTAB(5)
      NPART=REST(3)
      KORB=5
      KREST=KORB+(NSYM+1)*(NPART+1)
C Table consists of a 6-word header with data (see below), then
C an array dimensioned (NASPRT+2)*NFSB, and finally a hash map
C with suitable capacity e.g. 2*NFSB (50% usage). Each item in a
C hash map takes up 2 integers. Capacity must be at least NFSB+997.
      CALL mkVERTAB(NACTEL,MSPIN2,LSYM,NPART,REST(KORB),REST(KREST),
     &              SSTAB,NFSB0,NRDETS0,NFSB,NRDETS)
      NHEAD=7
      NSSTARR=(NASPRT+2)*NFSB
      NHSHMAP=997+2*NFSB
      NSIZE=NHEAD+NSSTARR+2*NHSHMAP
      CALL GETMEM('FSBTab','Allo','Inte',LFSBTAB,NSIZE)
      CALL ICOPY((NASPRT+2)*NFSB,FSBARR,1,IWORK(LFSBTAB+NHEAD),1)
      Call mma_deallocate(FSBARR)

      LSSTARR=LFSBTAB+NHEAD

C Position of hash table ('Map')
      KHSHMAP=1+NHEAD+NSSTARR
      LHSHMAP=LFSBTAB-1+KHSHMAP
      IWORK(LFSBTAB)=NSIZE
      IWORK(LFSBTAB+1)=ITYPE
      IWORK(LFSBTAB+2)=NFSB
      IWORK(LFSBTAB+3)=NASPRT
      IWORK(LFSBTAB+4)=NRDETS
      IWORK(LFSBTAB+5)=NHSHMAP
      IWORK(LFSBTAB+6)=KHSHMAP
C Make the hash map: NULL is a null marker. Suggested value=-1.
      NULL=-1
      CALL HSHINI(NHSHMAP,IWORK(LHSHMAP),NULL)
C Store values in the map:
      DO IFSB=1,NFSB
        CALL HSHPUT(NASPRT,NASPRT+2,IWORK(LSSTARR),
     &                              NHSHMAP,IWORK(LHSHMAP),IFSB)
      END DO
C Check that they can be obtained back:
      IERR=0
      DO IFSB=1,NFSB
        CALL HSHGET(IWORK(LSSTARR+(NASPRT+2)*(IFSB-1)),NASPRT,
     &        NASPRT+2,IWORK(LSSTARR),NHSHMAP,IWORK(LHSHMAP),JFSB)
        IF(IFSB.NE.JFSB) IERR=IERR+1
      END DO
      IF(IERR.GT.0) THEN
        WRITE(6,*)'NEWFSBTAB Hash index errors. IERR=',IERR
        CALL ABEND()
      END IF
      NEWFSBTAB=LFSBTAB
      RETURN
      END FUNCTION NEWFSBTAB
