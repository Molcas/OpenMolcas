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
      SUBROUTINE CHO_PMAT(P,X)
C
C     RASSCF version IBM-3090: SX section
C
c     Purpose: To construct from a canonically ordered list of
c              2-matrix elements a list ordered as the transformed
c              two-electron integrals. P is the input matrix and X is
c              the output matrix. the matrix is multiplied by two and
c              used to construct the Q-matrix in fock.
c
c      Used in order to obtain the right factor for the 2-body
c           density to be used in Cholesky (algo=2)
C
C          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
      DIMENSION X(*),P(*)
      COMMON /CHOPMAT / ipPL

C     CALL QENTER('PMAT')
C
c     Loop over all reordered 2-matrix elements.
C
      LPMAT=ISTORP(NSYM+1)
      CALL FZERO(X,LPMAT)
C
      IAT=0
      DO 14 NSP=1,NSYM
       NAP=NASH(NSP)
       IF(NAP.EQ.0) GO TO 14
       INDF=ISTORP(NSP)
       NUVX=(ISTORP(NSP+1)-INDF)/NAP
       LROW=0
       IAU=0
       DO 13 NSQ=1,NSYM
        NAQ=NASH(NSQ)
        IF(NAQ.EQ.0) GO TO 13
        NSPQ=IEOR(NSP-1,NSQ-1)
        IAV=0
        DO 12 NSR=1,NSYM
         NAR=NASH(NSR)
         IF(NAR.EQ.0) GO TO 12
         NSS=IEOR(NSPQ,NSR-1)+1
         IF(NSS.GT.NSR) GO TO 121
         NAS=NASH(NSS)
         IF(NAS.EQ.0) GO TO 121
         IAX=0
         IF(NSS.NE.1) THEN
          NSSM=NSS-1
          DO 7 NSS1=1,NSSM
          IAX=IAX+NASH(NSS1)
7         CONTINUE
         ENDIF
         DO 11 NAV=1,NAR
          LAV=NAV+IAV
          NAXE=NAS
          IF(NSR.EQ.NSS) NAXE=NAV
          DO 10 NAX=1,NAXE
           LAX=NAX+IAX
           DO 9 NAU=1,NAQ
            LAU=NAU+IAU
            LROW=LROW+1
            INDX=INDF+LROW-NUVX
            DO 8 NAT=1,NAP
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
             Work(ipPL+ntuvx-1) = 0.5d0*X(INDX)
8           CONTINUE
9          CONTINUE
10        CONTINUE
11       CONTINUE
121      IAV=IAV+NAR
12      CONTINUE
        IAU=IAU+NAQ
13     CONTINUE
       IAT=IAT+NAP
14    CONTINUE
C
      IF(IPR.GE.20) Write(IW,1000) (X(I),I=1,LPMAT)
1000  FORMAT(/1X,'Reordered 2-matrix'/(10F10.6))
C     CALL QEXIT('PMAT')
      RETURN
      END
