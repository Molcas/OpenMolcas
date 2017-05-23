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
* Copyright (C) 1996, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE H0INTSPC(  IH0SPC,  NPTSPC,IOCPTSPC,  NOCTPA,  NOCTPB,
     &                        IOCA,    IOCB,    NGAS, MXPNGAS,INTH0SPC,
     &                      NELFTP)
*
* Set up INTH0SPC : Division of CI space, so only determinants
*                   belonging to the same space  have nonvanishing
*                   matrix elements of H0
*
* =====
* Input
* =====
*
* IH0SPC : ne. 0 : Interacting subspaces have been defined
*          .eq.0 : Interacting subspaces not defined, let
*                  evrything interact
* NPTSPC : Number of subspaces defined
* IOCPTSPC : Allowed occumulated occupation of each subspace
* NOCTPA :  Number of alpha occupation types
* NOCTPB : Number of beta occupation types
* IOCA : Occupation  of alpha string
* IOCB : Occupation  of beta string
*
* Jeppe Olsen, January 1996
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      DIMENSION IOCPTSPC(2,MXPNGAS,*)
      DIMENSION IOCA(MXPNGAS,*),IOCB(MXPNGAS,*)
      DIMENSION NELFTP(*)
*. Output
      DIMENSION INTH0SPC(NOCTPA,NOCTPB)
*
      IF(IH0SPC.EQ.0) THEN
*. All interactions allowed
        IONE = 1
        CALL ISETVC(INTH0SPC,IONE,NOCTPA*NOCTPB)
      ELSE
*. Explicit construction of matrix giving partitionning of
*  subspaces
        IZERO = 0
        CALL ISETVC(INTH0SPC,IZERO,NOCTPA*NOCTPB)
*
        DO ISPC = 1, NPTSPC
          DO IATP = 1, NOCTPA
            DO IBTP = 1, NOCTPB
              IAMOKAY = 1
              IEL = 0
C?            WRITE(6,*) ' ISPC IATP IBTP ', ISPC,IATP,IBTP
              DO IGAS = 1, NGAS
               IEL = IEL
     &             + NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
C?             WRITE(6,*) ' IGAS IEL ', IGAS,IEL
C?             WRITE(6,*)
C?   &          ' Limits :',IOCPTSPC(1,IGAS,ISPC),IOCPTSPC(2,IGAS,ISPC)
               IF(IEL.LT.IOCPTSPC(1,IGAS,ISPC).OR.
     &            IEL.GT.IOCPTSPC(2,IGAS,ISPC)    ) IAMOKAY = 0
              END DO
C?            WRITE(6,*) ' IAMOKAY = ', IAMOKAY
*. Allowed
              IF(IAMOKAY.EQ.1.AND.INTH0SPC(IATP,IBTP).EQ.0) THEN
                INTH0SPC(IATP,IBTP) = ISPC
              END IF
            END DO
          END DO
        END DO
      END IF
*
      NTEST = 0
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Output from  H0INTSPC '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        CALL IWRTMA(INTH0SPC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
