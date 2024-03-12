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
      SUBROUTINE ALLOC()
!
!     RASSCF: allocation of core memory
!
!     Called from inpctl
!
!     No subroutine calls
!
!     ********** IBM-3090 Release 88 10 11 **********
!
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      IMPLICIT None
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
      Character(LEN=16), Parameter :: ROUTINE='ALLOC   '
      Integer IORD, IORP, NSP, NOP, NAP, NSQ, NSPQ, NSR, NAQ, NSPQR, NAR, NSS, NRS, NAS
#ifdef _DEBUGPRINT_
      WRITE(LF,*)' Entering ',ROUTINE
#endif
!
!     Compute space needed for transformed two-electron integrals
!
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
          IF(NSPQR.NE.NSS) Cycle
          NAS=NASH(NSS)
          NRS=NAR*NAS
          IF(NSS.EQ.NSR) NRS=(NAR+NAR**2)/2
          IORD=IORD+NOP*NAQ*NRS
          IORP=IORP+NAP*NAQ*NRS
         END DO
        END DO
       END DO
       ISTORD(NSP+1)=IORD
       ISTORP(NSP+1)=IORP
      END DO
      NFINT=ISTORD(NSYM+1)
!
#ifdef _DEBUGPRINT_
      Write(LF,'(1X,A,5X,9I5)')'ISTORD-vector:',(ISTORD(I),I=1,NSYM+1)
#endif
      END SUBROUTINE ALLOC
