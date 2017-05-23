************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE CMBN2DC(RNXYZ,NZETA,LA,LB,ZETA,RKAPPA,
     &                  FINAL,ALPHA,BETA,IFGRAD)
************************************************************************
*
* OBJECT: COMPUTE THE SECOND DERIVATIVE NON-ADIABATIC COUPLING
* MATRIX ELEMENTS, OF TYPE < D/DX CHI_1 | D/DX CHI_2 >
* WITH DIFFERENTIATION WRT NUCLEAR COORDINATES
*
* CALLED FROM: NONA2
*
* CALLING    : QENTER
*              DDOT_  (ESSL)
*              QEXIT
*
*     AUTHOR: PER AKE MALMQVIST NOV 2000
*             DEPT. OF THEORETICAL CHEMISTRY,
*             UNIVERSITY OF LUND, SWEDEN
*             FOLLOWING THE PATTERN OF R. LINDH,
*             SAME PLACE.
*
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FINAL(NZETA,(LA+1)*(LA+2)/2,(LB+1)*(LB+2)/2,1),
     &       ZETA(NZETA), RKAPPA(NZETA), BETA(NZETA),
     &       RNXYZ(NZETA,3,0:LA+1,0:LB+1), ALPHA(NZETA)
      LOGICAL IFGRAD(3)
* STATEMENT FUNCTION FOR CARTESIAN INDEX
      IND(IXYZ,IX,IZ) = (IXYZ-IX)*(IXYZ-IX+1)/2 + IZ + 1
*
C PREFACTOR FOR THE PRIMITIVE OVERLAP MATRIX
      DO IZETA=1,NZETA
        RKAPPA(IZETA)=RKAPPA(IZETA)*(ZETA(IZETA)**(-1.5D0))
      END DO

C LOOP STRUCTURE FOR THE CARTESIAN ANGULAR PARTS
      DO 10 IXA = 0, LA
         IYAMAX=LA-IXA
      DO 10 IXB = 0, LB
         IYBMAX=LB-IXB
         DO 20 IYA = 0, IYAMAX
            IZA = LA-IXA-IYA
            IPA= IND(LA,IXA,IZA)
         DO 20 IYB = 0, IYBMAX
            IZB = LB-IXB-IYB
            IPB= IND(LB,IXB,IZB)

* COMBINE 1-DIM PRIMITIVE OVERLAP INTEGRALS
      IF (IFGRAD(1)) THEN
C COMPUTE INTEGRALS TYPE <D/DX,D/DX>
         DO IZETA=1,NZETA
           DIFFX=4D0*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,1,IXA+1,IXB+1)
           IF(IXB.GT.0) THEN
             DIFFX=DIFFX-2D0*ALPHA(IZETA)*DBLE(IXB)*
     >             RNXYZ(IZETA,1,IXA+1,IXB-1)
             IF(IXA.GT.0) THEN
               DIFFX=DIFFX+DBLE(IXA*IXB)*RNXYZ(IZETA,1,IXA-1,IXB-1)
             END IF
           END IF
           IF(IXA.GT.0) THEN
             DIFFX=DIFFX-DBLE(2*IXA)*BETA(IZETA)*
     >             RNXYZ(IZETA,1,IXA-1,IXB+1)
           END IF
           OVLY=RNXYZ(IZETA,2,IYA,IYB)
           OVLZ=RNXYZ(IZETA,3,IZA,IZB)
           FINAL(IZETA,IPA,IPB,1)=RKAPPA(IZETA)*DIFFX*OVLY*OVLZ
         END DO
      END IF
      IF (IFGRAD(2)) THEN
C COMPUTE INTEGRALS TYPE <D/DY,D/DY>
         DO IZETA=1,NZETA
           DIFFY=4D0*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,2,IYA+1,IYB+1)
           IF(IYB.GT.0) THEN
             DIFFY=DIFFY-2D0*ALPHA(IZETA)*DBLE(IYB)*
     >             RNXYZ(IZETA,2,IYA+1,IYB-1)
             IF(IYA.GT.0) THEN
               DIFFY=DIFFY+DBLE(IYA*IYB)*RNXYZ(IZETA,2,IYA-1,IYB-1)
             END IF
           END IF
           IF(IYA.GT.0) THEN
             DIFFY=DIFFY-DBLE(2*IYA)*BETA(IZETA)*
     >             RNXYZ(IZETA,1,IYA-1,IYB+1)
           END IF
           OVLX=RNXYZ(IZETA,1,IXA,IXB)
           OVLZ=RNXYZ(IZETA,3,IZA,IZB)
           FINAL(IZETA,IPA,IPB,1)=RKAPPA(IZETA)*OVLX*DIFFY*OVLZ
         END DO
      END IF
      IF (IFGRAD(1)) THEN
C COMPUTE INTEGRALS TYPE <D/DZ,D/DZ>
         DO IZETA=1,NZETA
           DIFFZ=4D0*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,1,IZA+1,IZB+1)
           IF(IZB.GT.0) THEN
             DIFFZ=DIFFZ-2D0*ALPHA(IZETA)*DBLE(IZB)*
     >             RNXYZ(IZETA,1,IZA+1,IZB-1)
             IF(IZA.GT.0) THEN
               DIFFZ=DIFFZ+DBLE(IZA*IZB)*RNXYZ(IZETA,1,IZA-1,IZB-1)
             END IF
           END IF
           IF(IZA.GT.0) THEN
             DIFFZ=DIFFZ-DBLE(2*IZA)*BETA(IZETA)*
     >             RNXYZ(IZETA,1,IZA-1,IZB+1)
           END IF
           OVLX=RNXYZ(IZETA,1,IXA,IXB)
           OVLY=RNXYZ(IZETA,2,IYA,IYB)
           FINAL(IZETA,IPA,IPB,1)=RKAPPA(IZETA)*OVLX*OVLY*DIFFZ
         END DO
      END IF

C END OF LOOP NEST OVER CARTESIAN ANGULAR COMPONENT
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
