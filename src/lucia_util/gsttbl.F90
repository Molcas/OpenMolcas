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
      SUBROUTINE GSTTBL(       C,     CTT,    IATP,    IASM,    IBTP,   &
     &                      IBSM,   IOCOC,  NOCTPA,  NOCTPB,   NSASO,   &
     &                     NSBSO,  PSSIGN,  ICOOSC,     IDC,  PLSIGN,   &
     &                       LUC,     SCR,   NSMST,  ISCALE,  SCLFAC)
!****************************************************************************
! Variables status:
! C     = input CI vector
! CTT   = output CI vector in SD format
! obtain  determinant block (iatp iasm, ibtp ibsm )
! from vector packed in combination format according to IDC
!
!. If ISCALE = 1, the routine scales and returns the block
!  in determinant notmalization, and SCLFAC = 1.0D0
!
! If ISCALE = 0, the routine does not perform any overall
! scaling, and a scale factor is returned in SCLFAC
!
! IF ISCALE = 0, zero blocks are not set explicitly to zero,
! instead  zero is returned in SCLFAC
!
! ISCALE, SCLFAC added May 97
!
      use lucia_data, only: IDISK
      IMPLICIT NONE
      INTEGER IATP,IASM,IBTP,IBSM,NOCTPA,NOCTPB,IDC,LUC,NSMST,ISCALE
      REAL*8 PSSIGN,PLSIGN,SCLFAC
      REAL*8 C(*),CTT(*)
      INTEGER NSASO(NSMST, *),NSBSO(NSMST, *)
      INTEGER IOCOC(NOCTPA,NOCTPB),ICOOSC(NOCTPA,NOCTPB,*)
      REAL*8 SCR(*)

      INTEGER ISGVST(1)
      INTEGER IDUMMY(1)
      INTEGER LBL,NO_ZEROING,IMZERO,NAST,NBST,IBASE,NELMNT,NRI,NCI,     &
     &        IAMPACK,LCOMB,LDET
      REAL*8 PSIGN
!
!?    write(6,*) ' GSTTBL  ,IATP,IASM,IBTP,IBSM,ISCALE'
!?    write(6,*)            IATP,IASM,IBTP,IBSM,ISCALE
! =================
! Read in from disc
! =================
      IF(LUC.NE.0) THEN
        CALL IDAFILE(LUC,2,IDUMMY,1,IDISK(LUC))
        LBL=IDUMMY(1)
        CALL IDAFILE(LUC,2,IDUMMY,1,IDISK(LUC))
!?      write(6,*) ' LBL = ', LBL
         IF(ISCALE.EQ.1) THEN
            CALL FRMDSC(     SCR,     LBL,      -1,     LUC,  IMZERO,   &
     &                   IAMPACK)
         ELSE
            NO_ZEROING = 1
            CALL FRMDSC2(     SCR,     LBL,      -1,     LUC,  IMZERO,  &
     &                    IAMPACK,NO_ZEROING)
         END IF
!
         IF(IMZERO.EQ.1.AND.ISCALE.EQ.0) THEN
            SCLFAC = 0.0D0
         ELSE
            NAST = NSASO(IASM,IATP)
            NBST = NSBSO(IBSM,IBTP)
            IF(LBL.NE.0) THEN
               CALL SDCMRF(     CTT,     SCR,       2,    IATP,    IBTP,&
     &                         IASM,    IBSM,    NAST,    NBST,     IDC,&
     &                       PSSIGN,  PLSIGN,  ISGVST,    LDET,   LCOMB,&
     &                       ISCALE,  SCLFAC)
            ELSE
               SCLFAC = 0.0D0
            END IF
         END IF
!
!?      WRITE(6,*) ' ISCALE and SCLFAC on return in GSTTBL',
!?   &  ISCALE,SCLFAC

!. ISGVST and PLSIGN missing to make it work for IDC = 3,4
      ELSE
! =================
! Pack out from C
! =================
         IF(ISCALE.EQ.0) THEN
            WRITE(6,*) ' GSTTBL : LUC = 0 and ISCALE = 0'
            WRITE(6,*) ' I will scale as normal '
            SCLFAC = 1.0D0
         END IF
! Permutation sign
!. To get rid of annoying compiler warning
         PSIGN = 0.0D0
         IF(IDC.EQ.2) THEN
            PSIGN = PSSIGN
         ELSE IF(IDC .EQ. 3 ) THEN
            PSIGN = PLSIGN
         END IF
!        PLSSGN = PLSIGN * PSSIGN
! check for different packing possibilities and unpack
         IF(IASM.GT.IBSM.OR.IDC.EQ.1                                    &
     &        .OR.(IDC.EQ.3.AND.IASM.GE.IBSM))THEN
!*************
!* IASM > IBSM
!*************
            IF ( IDC.LT.4 ) THEN
!. Simple copy
               IBASE = ICOOSC(IATP,IBTP,IASM)
               NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
               CALL COPVEC(C(IBASE),CTT,NELMNT)
!?       write(6,*) ' simple copy IBASE NELMNT ',IBASE,NELMNT
!?       CALL WRTMAT(CTT,NSASO(IASM,IATP),NSBSO(IBSM,IBTP),
!?   &   NSASO(IASM,IATP),NSBSO(IBSM,IBTP))
!idc            ELSE IF( IDC.EQ.4 ) THEN
!. MLMS packed
!               IF(IATP.GT.IBTP) THEN
!                  IBASE = ICOOSC(IATP,IBTP,IASM)
!                  NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
!                  CALL COPVEC(C(IBASE),CTT,NELMNT)
!               ELSE IF(IATP.EQ.IBTP) THEN
!                  IBASE = ICOOSC(IATP,IATP,IASM)
!                  NAST = NSASO(IASM,IATP)
!                  CALL TRIPK3(CTT,C(IBASE),2,NAST,NAST,PLSIGN*PSSIGN)
!               ELSE IF( IATP.LT.IBTP) THEN
!                  IBASE = ICOOSC(IBTP,IATP,IASM)
!                  NROW  = NSASO(IASM,IBTP)
!                  NCOL  = NSBSO(IBSM,IATP)
!                  CALL TRPMT3(C(IBASE),NROW,NCOL,CTT)
!                  NELMNT = NROW*NCOL
!                  CALL SCALVE(CTT,PLSIGN*PSSIGN,NELMNT)
!idc               END IF
            END IF
         ELSE IF( IASM.EQ.IBSM) THEN
!*************
!* IASM = IBSM
!*************
            IF(IATP.GT.IBTP.OR.IDC.EQ.3) THEN
!.. simple copying
               IBASE = ICOOSC(IATP,IBTP,IASM)
               NELMNT = NSASO(IASM,IATP)*NSBSO(IBSM,IBTP)
               CALL COPVEC(C(IBASE),CTT,NELMNT)
            ELSE IF( IATP.EQ.IBTP) THEN
!.. expand triangular packed matrix
               IBASE = ICOOSC(IATP,IBTP,IASM)
               NAST = NSASO(IASM,IATP)
               CALL TRIPK3(     CTT,C(IBASE),       2,    NAST,    NAST,&
     &                       PSSIGN)
            ELSE IF( IATP .LT. IBTP) THEN
!.. transpose ibtp iasm iatp ibsm block
               IBASE = ICOOSC(IBTP,IATP,IASM)
               NRI = NSASO(IASM,IBTP)
               NCI = NSBSO(IASM,IATP)
               CALL TRPMT3(C(IBASE),NRI,NCI,CTT)
               IF(PSSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
            END IF
         ELSE IF( IASM .LT. IBSM ) THEN
!*************
!* IASM < IBSM
!*************
!.. transpose ibtp ibsm iatp iasm block
            IF(IDC.LT.4) THEN
               IBASE = ICOOSC(IBTP,IATP,IBSM)
               NRI = NSASO(IBSM,IBTP)
               NCI = NSBSO(IASM,IATP)
               IF( IDC.EQ.2) THEN
                  CALL TRPMT3(C(IBASE),NRI,NCI,CTT)
!               ELSE IF( IDC.EQ.3) THEN
!                  CALL COPVEC(C(IBASE),CTT,NRI*NCI)
               END IF
               IF(PSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
!idc            ELSE IF ( IDC .EQ. 4 ) THEN
!               IF(IBTP.GT.IATP) THEN
!                  IBASE = ICOOSC(IBTP,IATP,IBSM)
!                  NRI = NSASO(IBSM,IBTP)
!                  NCI = NSBSO(IASM,IATP)
!                  CALL TRPMT3(C(IBASE),NRI,NCI,CTT)
!                  IF(PSSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
!               ELSE IF (IBTP.EQ.IATP) THEN
!                  IBASE = ICOOSC(IBTP,IATP,IBSM)
!                  NRI   = NSASO(IBSM,IATP)
!                  NCI   = NSBSO(IASM,IATP)
!                  CALL TRIPK3(CTT,C(IBASE),2,NRI,NCI,PLSSGN)
!                  IF(PLSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NRI*NCI)
!               ELSE IF( IBTP.LT.IATP) THEN
!                  IBASE = ICOOSC(IATP,IBTP,IBSM)
!                  NELMNT = NSASO(IBSM,IATP)*NSBSO(IASM,IBTP)
!                  CALL COPVEC(C(IBASE),CTT,NELMNT)
!idc                  IF(PLSIGN.EQ.-1.0D0) CALL SCALVE(CTT,-1.0D0,NELMNT)
!               END IF
            END IF
         END IF
      END IF
!
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer_array(IOCOC)
      END SUBROUTINE GSTTBL
