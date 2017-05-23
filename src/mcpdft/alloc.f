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
      SUBROUTINE ALLOC_m
C
C     RASSCF: allocation of core memory
C
C     Called from inpctl
C
C     No subroutine calls
C
C     ********** IBM-3090 Release 88 10 11 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='ALLOC   ')
      IPRLEV=IPRLOC(1)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
C
C     Compute space needed for transformed two-electron integrals
C
!AMS
      ISTORD(1)=0
      ISTORP(1)=0
      IORD=0
      IORP=0
      DO NSP=1,NSYM
       NOP=NORB(NSP)
       NAP=NASH(NSP)
       DO NSQ=1,NSYM
        NAQ=NASH(NSQ)
        NSPQ=IEOR(NSP-1,NSQ-1)
        DO NSR=1,NSYM
         NSPQR=IEOR(NSPQ,NSR-1)+1
         NAR=NASH(NSR)
         DO NSS=1,NSR
          IF(NSPQR.NE.NSS) GO TO 11
          NAS=NASH(NSS)
          NRS=NAR*NAS
          IF(NSS.EQ.NSR) NRS=(NAR+NAR**2)/2
          IORD=IORD+NOP*NAQ*NRS
          IORP=IORP+NAP*NAQ*NRS
11       CONTINUE
         END DO
        END DO
       END DO
       ISTORD(NSP+1)=IORD
       ISTORP(NSP+1)=IORP
      END DO
      NFINT=ISTORD(NSYM+1)
C
      IF(IPRLEV.GE.DEBUG) THEN
       Write(LF,'(1X,A,5X,9I5)')'ISTORD-vector:',(ISTORD(I),I=1,NSYM+1)
      END IF
      RETURN
      END
