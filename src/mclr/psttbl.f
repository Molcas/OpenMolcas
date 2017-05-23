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
      SUBROUTINE PSTTBL_MCLR(C,CTT,IATP,IASM,IBTP,IBSM,IOCOC,
     &                  NOCTPA,NOCTPB,NSASO,NSBSO,PSIGN,
     &                  ICOOSC,IAC,IDC,LUHC,SCR)
*
* add(IAC = 1) or copy (IAC =2) determinant block (iatp iasm, ibtp ibsm
* to vector packed in combination format
* iatp,iasm , ibtp,ibsm is assumed to be allowed combination block
*
* Combination type is defined by IDC
* IAC = 2  does not work for LUHC.NE.0 !
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),CTT(*),NSASO(NOCTPA,*),NSBSO(NOCTPB,*)
      DIMENSION IOCOC(NOCTPA,NOCTPB),ICOOSC(NOCTPA,NOCTPB,*)
*
      DIMENSION SCR(*)
*
* ======================
* Write directly to disc
* ======================
*
*. Assumes complete block in,
*. copies to lower half, scales  and write out.
      IF(LUHC.NE.0) THEN
         NAST = NSASO(IATP,IASM)
         NBST = NSBSO(IBTP,IBSM)
         CALL SDCMRF_MCLR(CTT,SCR,1,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &               IDC,PSIGN,PLSIGN,ISGVST,LDET,LCOMB)
*. Note : PLSIGN and ISGVST missing in order to make it work for IDC=3,4
         CALL ITODS(LCOMB,1,-1,LUHC)
         CALL TODSC_MCLR(SCR,LCOMB,-1,LUHC)
      ELSE
* ==================
* Add to packed list
* ===================
      IF(IASM.GT.IBSM.OR.IDC.EQ.1
     &   .OR.IDC.EQ.3)THEN
**************
** IASM > IBSM
**************

        IF( IDC .LT. 4 ) THEN
*.. simple copying
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
          IF(IAC .EQ. 1 ) THEN
*           CALL VECSUM(C(IBASE),C(IBASE),CTT,1.0D0,1.0D0,NELMNT)
            Call DaXpY_(NELMNT,1.0d0,CTT,1,C(IBASE),1)
          ELSE IF(IAC .EQ.2 ) THEN
            CALL DCOPY_(NELMNT,CTT,1,C(IBASE),1)
          END IF
        ELSE IF ( IDC .EQ. 4 ) THEN
          IF(IATP.GT.IBTP) THEN
            IBASE = ICOOSC(IATP,IBTP,IASM)
            NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
            IF(IAC .EQ. 1 ) THEN
*             CALL VECSUM(C(IBASE),C(IBASE),CTT,1.0D0,1.0D0,NELMNT)
              Call DaXpY_(NELMNT,1.0d0,CTT,1,C(IBASE),1)
            ELSE IF(IAC .EQ.2 ) THEN
              CALL DCOPY_(NELMNT,CTT,1,C(IBASE),1)
            END IF
          ELSE IF( IATP .EQ. IBTP ) THEN
            IBASE = ICOOSC(IATP,IBTP,IASM)
            NAST = NSASO(IATP,IASM)
            IF( IAC .EQ. 1 ) THEN
              CALL PMPLFM(C(IBASE),CTT,NDIM)
            ELSE
              CALL TRIPK2(CTT,C(IBASE),1,NAST,NAST,PSIGN)
            END IF
          END IF
        END IF
      ELSE IF( IASM.EQ.IBSM) THEN
**************
** IASM = IBSM
**************
        IF(IATP.GT.IBTP) THEN
*.. simple copying
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
          IF(IAC .EQ. 1 ) THEN
*           CALL VECSUM(C(IBASE),C(IBASE),CTT,1.0D0,1.0D0,NELMNT)
            Call DaXpY_(NELMNT,1.0d0,CTT,1,C(IBASE),1)
          ELSE IF(IAC .EQ.2 ) THEN
            CALL DCOPY_(NELMNT,CTT,1,C(IBASE),1)
          END IF
        ELSE IF( IATP.EQ.IBTP) THEN
*.. reform to triangular packed matrix
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NAST = NSASO(IATP,IASM)
          IF( IAC .EQ. 1 ) THEN
            CALL PMPLFM(C(IBASE),CTT,NAST)
          ELSE
            CALL TRIPK2(CTT,C(IBASE),1,NAST,NAST,PSIGN)
          END IF
        END IF
      END IF
      END IF
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(IOCOC)
      END
