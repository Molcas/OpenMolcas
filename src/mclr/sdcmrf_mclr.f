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
      SUBROUTINE SDCMRF_MCLR(CSD,CCM,IWAY,IATP,IBTP,IASM,IBSM,NA,NB,
     &                  IDC,PS,PL,ISGVST,LDET,LCOMB)
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
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CSD(*),CCM(*),ISGVST(*)
*
      SQRT2  = SQRT(2.0D0)
      SQRT2I = 1.0D0/SQRT2
*
*. Is combination array packed ?
*
      IPACK = 0
      FACTOR = 1.0D0
*
      IF((IDC.EQ.2.OR.IDC.EQ.4).AND.IASM.EQ.IBSM) THEN
         SIGN = PS
         FACTOR = SQRT2
         IF(IATP.EQ.IBTP) IPACK = 1
      ELSE IF( IDC.EQ.4.AND.IASM.EQ.ISGVST(IBSM)) THEN
        IF(IATP.EQ.IBTP) IPACK = 1
        SIGN = PS*PL
        FACTOR = 2.0D0
      END IF
*
      LDET = NA * NB
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
          CALL TRIPK2(CSD,CCM,1,NA,NA,SIGN)
C              TRIPK2(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
        ELSE
          CALL DCOPY_(NA*NB,CSD,1,CCM,1)
        END IF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          CALL DSCAL_(LCOMB,FACTOR,CCM,1)
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
          CALL TRIPK2(CSD,CCM,2,NA,NA,SIGN)
        ELSE
           CALL DCOPY_(NA*NB,CCM,1,CSD,1)
        END IF
*. Scale
        IF(FACTOR.NE.1.0D0) THEN
          CALL DSCAL_(LDET,FACTOR,CSD,1)
          IF(IPACK.EQ.1) THEN
             CALL SCLDIA(CSD,SQRT2,NA,0)
          END IF
        END IF
      END IF
      RETURN
      END
