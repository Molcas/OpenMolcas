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
* Copyright (C) 1998, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE T_TO_NK_VECS(      T,   KORB,      C,  LUCIN, LUCOUT,
     &                          NSSOA,  NSSOB, NBLOCK, IBLOCK,   NAEL,
     &                           NBEL,  IASTR,  IBSTR,  IBLTP,  NSMST,
     &                         ICISTR,   NORB, IKAOCC, IKBOCC)
*
* Multiply Vector in LUCIN with t **NK_op to yield vector on LUCOUT
*
* Both files are initially rewinded
*
*
* Jeppe Olsen, Feb. 1998
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "io_util.fh"
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
*. Scratch
      DIMENSION C(*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IKAOCC(*),IKBOCC(*)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
      DIMENSION IBLTP(*)
      DIMENSION IDUM(1)
*
      IDISK(LUCIN)=0
      IDISK(LUCOUT)=0
*
      T2 = T**2
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
C?      WRITE(6,*) ' IATP IBTP IASM IBSM ', IATP,IBTP,IASM,IBSM
*. Obtain alpha strings of sym IASM and type IATP
        IDUM(1) = 0
        CALL GETSTR_TOTSM_SPGP(      1,   IATP,   IASM,   NAEL, NASTR1,
     &                           IASTR,   NORB,      0,   IDUM,   IDUM)
*. Occupation of orb KORB
        DO JSTR = 1, NASTR1
          KOCC = 0
          DO JAEL = 1, NAEL
            IF(IASTR(JAEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKAOCC(JSTR) = KOCC
        END DO
C?      WRITE(6,*) ' IKAOCC array '
C?      CALL IWRTMA(IKAOCC,1,NASTR1,1,NASTR1)


*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM(1) = 0
        CALL GETSTR_TOTSM_SPGP(      2,   IBTP,   IBSM,   NBEL, NBSTR1,
     &                           IBSTR,   NORB,      0,   IDUM,   IDUM)
C?      WRITE(6,*) ' After GETSTR, NBSTR1=',NBSTR1
*. Occupation of orb KORB
        DO JSTR = 1, NBSTR1
C?        write(6,*) ' JSTR = ', JSTR
          KOCC = 0
          DO JBEL = 1, NBEL
C?          write(6,*) JBEL, IBSTR(JBEL,JSTR)
            IF(IBSTR(JBEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKBOCC(JSTR) = KOCC
        END DO
C?      WRITE(6,*) ' IKBOCC array '
C?      CALL IWRTMA(IKBOCC,1,NBSTR1,1,NBSTR1)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
C?      WRITE(6,*) ' IBLTP ', IBLTP(IASM)
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
C?      WRITE(6,*) ' NIA NIB ', NIA,NIB
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
          CALL IDAFILE(LUCIN,2,IDUM,1,IDISK(LUCIN))
          LDET=IDUM(1)
          CALL IDAFILE(LUCIN,2,IDUM,1,IDISK(LUCIN))
          CALL FRMDSC(        C,     LDET,       -1,    LUCIN,   IMZERO,
     &                  IAMPACK)
        END IF
        IF(IMZERO.NE.1) THEN
*
          IDET = 0
          DO  IB = 1,NIB
            IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
              MINIA = IB
            ELSE
              MINIA = 1
            END IF
            DO  IA = MINIA,NIA
*
              IDET = IDET + 1
C?            WRITE(6,*) ' IA IB IDET',IA,IB,IDET
              KABOCC = IKAOCC(IA)+IKBOCC(IB)
              IF(KABOCC.EQ.1) THEN
                C(IDET) = T*C(IDET)
              ELSE IF(KABOCC.EQ.2) THEN
                C(IDET) = T2 *C(IDET)
              END IF
            END DO
*           ^ End of loop over alpha strings
          END DO
*         ^ End of loop over beta strings
*
        END IF
*       ^ End of if statement for nonvanishing blocks
*. Save result on LUCOUT
        CALL ITODS([LDET],1,-1,LUCOUT)
        CALL TODSC(C,LDET,-1,LUCOUT)
      END DO
*     ^ End of loop over blocks
*. This is the end, the end of every file my friend, the end
       CALL ITODS([-1],1,-1,LUCOUT)

      RETURN
      END
*
