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
      SUBROUTINE SPSPCLS_GAS( NOCTPA, NOCTPB,   IOCA,   IOCB, NELFGP,
     &                       MXPNGAS,   NGAS,ISPSPCLS,  ICLS,  NCLS,
     &                         IPRNT)
*
* Obtain mapping a-supergroup X b-supergroup => class
*
*
* =====
*.Input
* =====
*
* NOCTPA : Number of alpha types
* NOCTPB : Number of beta types
*
* IOCA(IGAS,ISTR) occupation of AS IGAS for alpha string type ISTR
* IOCB(IGAS,ISTR) occupation of AS IGAS for beta  string type ISTR
*
* MXPNGAS : Largest allowed number of gas spaces
* NGAS    : Actual number of gas spaces

*
* ======
*.Output
* ======
*
* ISPSPCLS(IATP,IBTP) => Class of this block of determinants
*                        =0 indicates unallowed(class less) combination
*
*.Input
      INTEGER IOCA(MXPNGAS,NOCTPA),IOCB(MXPNGAS,NOCTPB)
      INTEGER NELFGP(*)
      INTEGER ICLS(NGAS,NCLS)
*.Output
      INTEGER ISPSPCLS(NOCTPA,NOCTPB)

*
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ISPSPCLS_GAS entered '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        WRITE(6,*) ' IOCA and IOCB '
        CALL IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NGAS)
        CALL IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' ICLS '
        CALL IWRTMA(ICLS,NGAS,NCLS,NGAS,NCLS)
      END IF
*
      DO 100 IATP = 1, NOCTPA
        DO 90 IBTP = 1, NOCTPB
          IICLS = 0
          DO KCLS = 1, NCLS
            IAMOKAY = 1
            DO IGAS = 1, NGAS
              IEL = NELFGP(IOCA(IGAS,IATP))+NELFGP(IOCB(IGAS,IBTP))
              IF(IEL.NE.ICLS(IGAS,KCLS)) IAMOKAY = 0
            END DO
            IF(IAMOKAY.EQ.1) IICLS=KCLS
          END DO
          ISPSPCLS(IATP,IBTP) = IICLS
   90   CONTINUE
  100 CONTINUE
*
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Matrix giving classes for alpha-beta supergroups'
        WRITE(6,*)
        CALL IWRTMA(ISPSPCLS,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
C     BLKCLS(WORK(KLCIBT),NBLOCKS,WORK(KLBLKCLS),WORK(KLSPSPCL),
C    &            NOCTPA,NOCTPB)
