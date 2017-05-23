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
      SUBROUTINE SDCMRF(     CSD,     CCM,    IWAY,    IATP,    IBTP,
     &                      IASM,    IBSM,      NA,      NB,     IDC,
     &                        PS,      PL,  ISGVST,    LDET,   LCOMB,
     &                    ISCALE,  SCLFAC)
*
* Change a block of coefficients bwtween combination format and
* Slater determinant format
*
*     IWAY = 1 : SD => Combinations
*     IWAY = 2 : Combinations => SD
*
* Input
* =====
* CSD : Block in determinant form
* CCM : Block in combination  form
* IWAY : as above
* IATP,IBTP : type of alpha- and beta- string
* NA,NB : Number of alpha- and beta- strings
* IDC  : Combination type
* PS   : Spin combination sign
* PL   : Ml   combination sign
* ISGVST : Ml reflection of strings
*
*
* If ISCALE .EQ. 0, no overall scaling is performed,
*                   the overall scale factor is returned
*                   as SCLFAC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CSD(*),CCM(*),ISGVST(*)
*
      NTEST = 00
*
      SQRT2  = SQRT(2.0D0)
      SQRT2I = 1.0D0/SQRT2
*
*. Is combination array packed ?
*
      SCLFAC = 1.0D0
      IPACK = 0
      FACTOR = 1.0D0
*
      IF(IDC.EQ.2.OR.IDC.EQ.4) THEN
         SIGN = PS
         FACTOR = SQRT2
         IF(IASM.EQ.IBSM.AND.IATP.EQ.IBTP) IPACK = 1
      ELSE IF( IDC.EQ.4.AND.IASM.EQ.ISGVST(IBSM)) THEN
        IF(IATP.EQ.IBTP) IPACK = 1
        SIGN = PS*PL
        FACTOR = 2.0D0
      END IF
*
      LDET = NA * NB
      IF(NTEST.GE.100) WRITE(6,*) ' SDCMRF : NA, NB =', NA,NB
      IF( IPACK.EQ.0) THEN
        LCOMB = LDET
      ELSE
        LCOMB = NA*(NA+1)/2
      END IF
      IF(IDC.EQ.4.AND.IPACK.EQ.0) FACTOR = SQRT2
      IF(IWAY.EQ.2) FACTOR = 1.0D0/FACTOR
*
*. SD => combination transformation
*
      IF(IWAY .EQ. 1 ) THEN
        IF(IPACK.EQ.1) THEN
*. Pack to triangular form
          CALL TRIPK3(      CSD,      CCM,        1,       NA,       NA,
     &                     SIGN)
C              TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
        ELSE
          CALL COPVEC(CSD,CCM,NA*NB)
        END IF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          IF(ISCALE.EQ.1) THEN
            SCLFAC = 1.0D0
            CALL SCALVE(CCM,FACTOR,LCOMB)
          ELSE
            SCLFAC = FACTOR
          END IF
          IF(IPACK.EQ.1 ) THEN
            CALL SCLDIA(CCM,SQRT2I,NA,1)
          END IF
        END IF
      END IF
*
*. Combination => SD transformation
*
      IF(IWAY.EQ.2) THEN
        IF(IPACK.EQ.1) THEN
*. Unpack from triangular form
          CALL TRIPK3(      CSD,      CCM,        2,       NA,       NA,
     &                     SIGN)
        ELSE
           CALL COPVEC(CCM,CSD,NA*NB)
        END IF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          IF(ISCALE.EQ.1) THEN
            SCLFAC = 1.0D0
            CALL SCALVE(CSD,FACTOR,LDET)
          ELSE
            SCLFAC = FACTOR
          END IF
          IF(IPACK.EQ.1) THEN
             CALL SCLDIA(CSD,SQRT2,NA,0)
          END IF
        END IF
      END IF
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
C     IF(NTEST.NE.0.AND.IWAY.EQ.1) THEN
        WRITE(6,*) ' Information from SDCMRF '

        WRITE(6,'(A,6I4)') ' IWAY IATP IBTP IASM IBSM IDC ',
     &                   IWAY,IATP,IBTP,IASM,IBSM,IDC
        WRITE(6,'(A,I4,3X,2E15.8)') ' IPACK FACTOR SIGN',
     &  IPACK,FACTOR,SIGN
        IF(NTEST.GE. 100 ) THEN
          WRITE(6,*) ' Slater determinant block '
          CALL WRTMAT(CSD,NA,NB,NA,NB)
          WRITE(6,*)
          WRITE(6,*) ' Combination block '
          IF(IPACK.EQ.1) THEN
            CALL PRSM2(CCM,NA)
          ELSE
            CALL WRTMAT(CCM,NA,NB,NA,NB)
          END IF
        END IF
      END IF
*
      RETURN
      END
