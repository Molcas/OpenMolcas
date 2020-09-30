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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GETSTR_TOTSM_SPGP(ISTRTP,ISPGRP,ISPGRPSM, NEL,NSTR,
     &                               ISTR, NORBT,IDOREO,    IZ,  IREO)
*
* Obtain all super-strings of given total symmetry and given
* occupation in each GAS space
*
*.If  IDOREO .NE. 0 THEN reordering array : lexical => actual order is obtained
*
* Nomenclature of the day : superstring : string in complete
*                           orbital space, product of strings in
*                           each GAS space
* =====
* Input
* =====
*
* ISTRTP  : Type of of superstrings ( alpha => 1, beta => 2 )
* ISPGRP :  supergroup number, (relative to start of this type )
* ISPGRPSM : Total symmetry of superstrings
* NEL : Number of electrons
* IZ  : Reverse lexical ordering array for this supergroup
*
*
* ======
* Output
* ======
*
* NSTR : Number of superstrings generated
* ISTR : Occupation of superstring
* IREO : Reorder array ( if IDOREO.NE.0)
*
*
* Jeppe Olsen, July 1995
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "strbas.fh"
#include "csm.fh"
      INTEGER IZ(NORBT,NEL)
*. output
      INTEGER ISTR(*), IREO(*)
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================== '
        WRITE(6,*) ' Welcome to GETSTR_TOTSM_SPGP '
        WRITE(6,*) ' ============================== '
        WRITE(6,*)
        WRITE(6,'(A,3I3)')
     & ' Strings to be obtained : Type, supergroup, symmetry ',
     &   ISTRTP,ISPGRP,ISPGRPSM
        WRITE(6,*)
      END IF
*. Absolut number of this supergroup
      ISPGRPA = IBSPGPFTP(ISTRTP) - 1 + ISPGRP
*. Occupation per gasspace
*. Largest occupied space
      NGASL = 0
*. Largest and lowest symmetries active in each GAS space
      DO IGAS = 1, NGAS
        ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPA)
        NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
        IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
      IF(NGASL.EQ.0) NGASL = 1
*. Number of strings per GAS space and offsets for strings of given sym
      DO IGAS = 1, NGAS
        CALL ICOPVE2(iWORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(iWORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST =1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
* Largest and lowest active symmetries for each GAS space
      IF(NTEST.GE.200) THEN
         WRITE(6,*) ' Type of each GAS space '
         CALL IWRTMA(ITPFGS,1,NGAS,1,NGAS)
         WRITE(6,*) ' Number of elecs per GAS space '
         CALL IWRTMA(NELFGS,1,NGAS,1,NGAS)
      END IF
*
*. Loop over symmetries of each GAS
*
      MAXLEX = 0
      IFIRST = 1
      ISTRBS = 1
 1000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
C        NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
C        CALL NXTNUM2(ISMFGS,NGASL-1,1,MAXVAL,NONEW)
C        CALL NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        IFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
C         SYMCOM(ITASK,IOBJ,I1,I2,I12)
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
C?        write(6,*) ' ISTSMM1 : ', ISTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Obtain all strings of this symmetry
        CALL GETSTRN_GASSM_SPGP( ISMFGS,
     &                           ITPFGS,
     &                          ISTR(1+NEL*(ISTRBS-1)),
     &                             NSTR,
     &                              NEL,
*
     &                          NNSTSGP,
     &                          IISTSGP)
*. Reorder Info : Lexical => actual number
        IF(IDOREO.NE.0) THEN
*. Lexical number of NEL electrons
*. Can be made smart by using common factor for first NGAS-1 spaces
          DO JSTR = ISTRBS, ISTRBS+NSTR-1
            LEX = 1
            DO IEL = 1, NEL
              LEX = LEX + IZ(ISTR(IEL+NEL*(JSTR-1)),IEL)
            END DO
C           WRITE(6,*) ' string '
C           CALL IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
C           WRITE(6,*) ' JSTR and LEX ', JSTR,LEX
*
            MAXLEX = MAX(MAXLEX,LEX)
            IREO(LEX) = JSTR
          END DO
        END IF
*
        ISTRBS = ISTRBS + NSTR
*. ready for next symmetry distribution
        IF(NGAS-1.NE.0) GOTO 1000
 1001 CONTINUE
*. End of loop over symmetry distributions
      NSTR = ISTRBS - 1
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of strings generated ', NSTR
        WRITE(6,*)
        WRITE(6,*) ' Strings : '
        WRITE(6,*)
        CALL PRTSTR(ISTR,NEL,NSTR)
*
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) 'Largest Lexical number obtained ', MAXLEX
          WRITE(6,*) ' Reorder array '
          CALL IWRTMA(IREO,1,NSTR,1,NSTR)
        END IF
      END IF
*
      RETURN
      END
