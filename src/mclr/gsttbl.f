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
      SUBROUTINE GSTTBL_MCLR(C,CTT,IATP,IASM,IBTP,IBSM,IOCOC,
     &                  NOCTPA,NOCTPB,NSASO,NSBSO,PSSIGN,ICOOSC,IDC,
     &                  PLSIGN,LUC,SCR)
*
* obtain  determinant block (iatp iasm, ibtp ibsm )
* from vector packed in combination format according to IDC
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(*),CTT(*),NSASO(NOCTPA,*),NSBSO(NOCTPB,*)
      DIMENSION IOCOC(NOCTPA,NOCTPB),ICOOSC(NOCTPA,NOCTPB,*)
      DIMENSION SCR(*)
      DIMENSION ISGVST(IBSM)
      DIMENSION IDUM(1)
*
      PSIGN=0.0D0 ! dummy initialize
*
* =================
* Read in from disc
* =================
      IF(LUC.NE.0) THEN
        CALL IFRMDS(IDUM,1,-1,LUC)
        LBL=IDUM(1)
        CALL FRMDSC_MCLR(SCR,LBL,-1,LUC,IMZERO)
        NAST = NSASO(IATP,IASM)
        NBST = NSBSO(IBTP,IBSM)
        IF(LBL.NE.0)
     &  CALL SDCMRF_MCLR(CTT,SCR,2,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &              IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB)
*. ISGVST and PLSIGN missing to make it work for IDC = 3,4
      ELSE
* =================
* Pack out from C
* =================
* Permutation sign
      IF(IDC.EQ.2) THEN
        PSIGN = PSSIGN
      ELSE IF(IDC .EQ. 3 ) THEN
        PSIGN = PLSIGN
      ELSE
         PSIGN=0.0D0
      END IF
      PLSSGN = PLSIGN * PSSIGN
* check for different packing possibilities and unpack
      IF(IASM.GT.IBSM.OR.IDC.EQ.1
     &   .OR.(IDC.EQ.3.AND.IASM.GE.IBSM))THEN
**************
** IASM > IBSM
**************
        IF ( IDC.LT.4 ) THEN
*. Simple copy
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
          CALL DCOPY_(NELMNT,C(IBASE),1,CTT,1)
        ELSE IF( IDC.EQ.4 ) THEN
*. MLMS packed
          IF(IATP.GT.IBTP) THEN
            IBASE = ICOOSC(IATP,IBTP,IASM)
            NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
            CALL DCOPY_(NELMNT,C(IBASE),1,CTT,1)
          ELSE IF(IATP.EQ.IBTP) THEN
            IBASE = ICOOSC(IATP,IATP,IASM)
            NAST = NSASO(IATP,IASM)
            CALL TRIPK2(CTT,C(IBASE),2,NAST,NAST,PLSIGN*PSSIGN)
          ELSE IF( IATP.LT.IBTP) THEN
            IBASE = ICOOSC(IBTP,IATP,IASM)
            NROW  = NSASO(IBTP,IASM)
            NCOL  = NSBSO(IATP,IBSM)
            CALL TRPMAT(C(IBASE),NROW,NCOL,CTT)
            NELMNT = NROW*NCOL
            CALL DSCAL_(NELMNT,PLSIGN*PSSIGN,CTT,1)
          END IF
        END IF
      ELSE IF( IASM.EQ.IBSM) THEN
**************
** IASM = IBSM
**************
        IF(IATP.GT.IBTP.OR.IDC.EQ.3) THEN
*.. simple copying
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NELMNT = NSASO(IATP,IASM)*NSBSO(IBTP,IBSM)
          CALL DCOPY_(NELMNT,C(IBASE),1,CTT,1)
        ELSE IF( IATP.EQ.IBTP) THEN
*.. expand triangular packed matrix
          IBASE = ICOOSC(IATP,IBTP,IASM)
          NAST = NSASO(IATP,IASM)
          CALL TRIPK2(CTT,C(IBASE),2,NAST,NAST,PSSIGN)
        ELSE IF( IATP .LT. IBTP) THEN
*.. transpose ibtp iasm iatp ibsm block
          IBASE = ICOOSC(IBTP,IATP,IASM)
          NRI = NSASO(IBTP,IASM)
          NCI = NSBSO(IATP,IASM)
          CALL TRPMAT(C(IBASE),NRI,NCI,CTT)
          IF(PSSIGN.EQ.-1.0D0) CALL DSCAL_(NRI*NCI,-1.0d0,CTT,1)
        END IF
      ELSE IF( IASM .LT. IBSM ) THEN
**************
** IASM < IBSM
**************
*.. transpose ibtp ibsm iatp iasm block
        IF(IDC.LT.4) THEN
          IBASE = ICOOSC(IBTP,IATP,IBSM)
          NRI = NSASO(IBTP,IBSM)
          NCI = NSBSO(IATP,IASM)
          IF( IDC.EQ.2) THEN
            CALL TRPMAT(C(IBASE),NRI,NCI,CTT)
          ELSE IF( IDC.EQ.3) THEN
            CALL DCOPY_(NRI*NCI,C(IBASE),1,CTT,1)
          END IF
          IF(PSIGN.EQ.-1.0D0) CALL DSCAL_(NRI*NCI,-1.0D0,CTT,1)
        ELSE IF ( IDC .EQ. 4 ) THEN
          IF(IBTP.GT.IATP) THEN
            IBASE = ICOOSC(IBTP,IATP,IBSM)
            NRI = NSASO(IBTP,IBSM)
            NCI = NSBSO(IATP,IASM)
            CALL TRPMAT(C(IBASE),NRI,NCI,CTT)
            IF(PSSIGN.EQ.-1.0D0) CALL DSCAL_(NRI*NCI,-1.0D0,CTT,1)
          ELSE IF (IBTP.EQ.IATP) THEN
            IBASE = ICOOSC(IBTP,IATP,IBSM)
            NRI   = NSASO(IATP,IBSM)
            NCI   = NSBSO(IATP,IASM)
            CALL TRIPK2(CTT,C(IBASE),2,NRI,NCI,PLSSGN)
            IF(PLSIGN.EQ.-1.0D0) CALL DSCAL_(NRI*NCI,-1.0D0,CTT,1)
          ELSE IF( IBTP.LT.IATP) THEN
            IBASE = ICOOSC(IATP,IBTP,IBSM)
            NELMNT = NSASO(IATP,IBSM)*NSBSO(IBTP,IASM)
            CALL DCOPY_(NELMNT,C(IBASE),1,CTT,1)
            IF(PLSIGN.EQ.-1.0D0) CALL DSCAL_(NELMNT,-1.0D0,CTT,1)
          END IF
        END IF
      END IF
      END IF
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(IOCOC)
      END
