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
!#define _DEBUGPRINT_
      SUBROUTINE PMAT_RASSCF(P,X)
C
C     RASSCF version IBM-3090: SX section
C
c     Purpose: To construct from a canonically ordered list of
c              2-matrix elements a list ordered as the transformed
c              two-electron integrals. P is the input matrix and X is
c              the output matrix. the matrix is multiplied by two and
c              used to construct the Q-matrix in fock.
C
C          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
C
      use general_data
#ifdef _DEBUGPRINT_
      use Definitions, only: LF => u6
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasscf.fh"
      Real*8 X(*),P(*)
#ifdef _DEBUGPRINT_
      Character(LEN=16), Parameter :: ROUTINE='PMAT            '

C Local print level (if any)
      WRITE(LF,*)' Entering ',ROUTINE
#endif
C
c     Loop over all reordered 2-matrix elements.
C
      LPMAT=ISTORP(NSYM+1)
      X(1:LPMAT)=0.0D0
C
      IAT=0
      DO NSP=1,NSYM
       NAP=NASH(NSP)
       IF(NAP.EQ.0) Cycle
       INDF=ISTORP(NSP)
       NUVX=(ISTORP(NSP+1)-INDF)/NAP
       LROW=0
       IAU=0
       DO NSQ=1,NSYM
        NAQ=NASH(NSQ)
        IF(NAQ.EQ.0) Cycle
        NSPQ=IEOR(NSP-1,NSQ-1)
        IAV=0
        DO NSR=1,NSYM
         NAR=NASH(NSR)
         IF(NAR.EQ.0) Cycle
         NSS=IEOR(NSPQ,NSR-1)+1
         IF(NSS.GT.NSR) Then
           IAV=IAV+NAR
           Cycle
         ENDIF
         NAS=NASH(NSS)
         IF(NAS.EQ.0) Then
           IAV=IAV+NAR
           Cycle
         ENDIF
         IAX=0
         IF(NSS.NE.1) THEN
          NSSM=NSS-1
          DO NSS1=1,NSSM
           IAX=IAX+NASH(NSS1)
          END DO
         ENDIF
         DO NAV=1,NAR
          LAV=NAV+IAV
          NAXE=NAS
          IF(NSR.EQ.NSS) NAXE=NAV
          DO NAX=1,NAXE
           LAX=NAX+IAX
           DO NAU=1,NAQ
            LAU=NAU+IAU
            LROW=LROW+1
            INDX=INDF+LROW-NUVX
            DO NAT=1,NAP
             INDX=INDX+NUVX
             LAT=NAT+IAT
C
c            Compute canonical index ntuvx and find prefactor
C
             LAT1=MAX(LAT,LAU)
             LAU1=MIN(LAT,LAU)
             NTU=ITRI(LAT1)+LAU1
             NVX=ITRI(LAV)+LAX
             NTUVX=ITRI(MAX(NTU,NVX))+MIN(NTU,NVX)
             FAC=2.0D0
             IF(NTU.LT.NVX) THEN
              IF(LAT1.EQ.LAU1.AND.LAV.NE.LAX) FAC=4.0D0
              IF(LAT1.NE.LAU1.AND.LAV.EQ.LAX) FAC=1.0D0
             ENDIF
             X(INDX)=FAC*P(NTUVX)
            END DO
           END DO
          END DO
         END DO
         IAV=IAV+NAR
        END DO
        IAU=IAU+NAQ
       END DO
       IAT=IAT+NAP
       END DO
C
#ifdef _DEBUGPRINT_
      Write(LF,*)' Reordered 2-matrix:'
      Write(LF,'(1X,10F10.6)') (X(I),I=1,LPMAT)
#endif
      END SUBROUTINE PMAT_RASSCF
