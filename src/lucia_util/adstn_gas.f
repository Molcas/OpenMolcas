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
* Copyright (C) 1991,1994-1996, Jeppe Olsen                            *
************************************************************************
      SUBROUTINE ADSTN_GAS(  KLOFFI,   IOBSM,   IOBTP,   ISPGP, ISPGPSM,
     &                      ISPGPTP,      I1,    XI1S,   NKSTR,    IEND,
     &                        IFRST,   KFRST,    KACT,  SCLFAC)
*
*
* Obtain mappings
* a+IORB !KSTR> = +/-!ISTR> for orbitals of symmetry IOBSM and type IOBTP
* and I strings belonging to supergroup ISPGP wih symmetry ISPGPSM
* and type ISPGPTP(=1=>alpha,=2=>beta)
*
* The results are given in the form
* I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR>
* (numbering relative to TS start)
* Above +/- is stored in XI1S
*
* if some nonvanishing excitations were found, KACT is set to 1,
* else it is zero
*
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*               August 95    : GAS version
*               October 96   : Improved version
*
* ======
*. Input
* ======
*
*./BIGGY
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "strbas.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "csm.fh"
#include "lucinp.fh"
*. Local scratch
      INTEGER NELFGS(MXPNGAS), ISMFGS(MXPNGAS),ITPFGS(MXPNGAS)
      INTEGER MAXVAL(MXPNGAS),MINVAL(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
*
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
C??      DIMENSION IOFFI(LOFFI)
      PARAMETER(MXLNGAS=20)
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Will be stored as an matrix of dimension
* (NKSTR,*), Where NKSTR is the number of K-strings of
*  correct symmetry . Nk is provided by this routine.
*
* PAM Mars-2006: Allocation moved to outside this subroutine.
*      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'ADSTN ')
*      CALL MEMMAN(KLOFFI,LOFFI,'ADDL  ',2,'KLOFFI')
*
      IF(NGAS.GT.MXLNGAS) THEN
        WRITE(6,*) ' Ad hoc programming in ADSTN (IOFFI)'
        WRITE(6,*) ' Must be changed - or redimensioned '
*        STOP'ADST : IOFFI problem '
        CALL SYSABENDMSG('lucia_util/adstn_gas',
     &                          'Internal error',' ')
      END IF
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADSTN_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*) '  IOBTP IOBSM : ', IOBTP,IOBSM
        WRITE(6,*) '  ISPGP ISPGPSM ISPGPTP :  ',
     &                ISPGP,ISPGPSM,ISPGPTP
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Problemo : ADSTN_GAS'
C?      WRITE(6,*) ' SCLFAC .ne. 1 '
C?    END IF
*
*. Supergroup and symmetry of K strings
*
      ISPGRPABS = IBSPGPFTP(ISPGPTP)-1+ISPGP
      CALL NEWTYP(ISPGRPABS,1,IOBTP,KSPGRPABS)
      CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
      NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
      IF(NTEST.GE.200) WRITE(6,*)
     & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
      IF(NKSTR.EQ.0) GOTO 9999
*
      NORBTS= NOBPTS(IOBTP,IOBSM)
      ZERO =0.0D0
      CALL SETVEC(XI1S,ZERO,NORBTS*NKSTR)
      IZERO = 0
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
*. First orbital of given GASSpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1
*. First orbital of fiven GASSPace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM)


*
*. Information about I strings
* =============================
*
*. structure of group of strings defining I strings
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Number of electrons before active type
      NELB = 0
      DO IGAS = 1, IOBTP -1
        NELB = NELB + NELFGS(IGAS)
      END DO
*. Number of electrons in active space
      NACGSOB = NOBPT(IOBTP)

*. Number of strings per symmetry for each symmetry
      DO IGAS = 1, NGAS
        CALL ICOPVE2(iWORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
      END DO
*. Offset and dimension for active group in I strings
      CALL ICOPVE2(iWORK(KISTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               IACIST)
      CALL ICOPVE2(iWORK(KNSTSGP(1)),(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,
     &               NACIST)
C?     WRITE(6,*) ' IACIST and NACIST arrays '
C?     CALL IWRTMA(IACIST,1,NSMST,1,NSMST)
C?     CALL IWRTMA(NACIST,1,NSMST,1,NSMST)
*
*. Generate offsets for I strings with given symmetry in
*  each space
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
      IFIRST = 1
      NSTRINT = 0
 2000 CONTINUE
        IF(IFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 2001
        END IF
        IFIRST = 0
*. Symmetry of NGASL -1 spaces given, symmetry of full space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*.  sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,ISPGPSM)
        ISMFGS(NGASL) = ISMGSN
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL spaces '
          CALL IWRTMA(ISMFGS,1,NGASL,1,NGASL)
        END IF
*. Number of strings with this symmetry combination
        NSTRII = 1
        DO IGAS = 1, NGASL
          NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for this symmetry distribution in IOFFI
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGASL
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
*
      IF(NTEST.GE.1) THEN !SJS
        WRITE(6,*)
        WRITE(6,*) ' ============================ '
        WRITE(6,*) ' If program is crashing here,'
        WRITE(6,*) ' LOFFI needs to be increased. '
        WRITE(6,*) ' ============================ '
        WRITE(6,*)
      END IF
        WORK(KLOFFI+IOFF-1) = DBLE(NSTRINT) + 1.001D0
        NSTRINT = NSTRINT + NSTRII
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' IOFF, IOFFI(IOFF) NSTRII ',
     &                 IOFF, WORK(KLOFFI+IOFF-1),NSTRII
        END IF
*
      IF(NGASL-1.GT.0) GOTO 2000
 2001 CONTINUE


*
*. Supergroup and symmetry of K strings
*
CM    CALL NEWTYP(ISPGRPABS,1,IOBTP,KSPGRPABS)
CM    CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
CM    NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
CM    IF(NTEST.GE.200) WRITE(6,*)
CM   & ' KSM, KSPGPRABS, NKSTR : ', KSM,KSPGRPABS, NKSTR
*
*. Gas structure of K strings
*
      NGASL = 1
      DO IGAS = 1, NGAS
       ITPFGS(IGAS) = ISPGPFTP(IGAS,KSPGRPABS)
       NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
       IF(NELFGS(IGAS).GT.0) NGASL = IGAS
      END DO
*. Active group of K-strings
      KACGRP = ITPFGS(IOBTP)
*. Number of strings per symmetry distribution
      DO IGAS = 1, NGAS
        CALL ICOPVE2(iWORK(KNSTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               NNSTSGP(1,IGAS))
        CALL ICOPVE2(iWORK(KISTSGP(1)),(ITPFGS(IGAS)-1)*NSMST+1,NSMST,
     &               IISTSGP(1,IGAS))
      END DO
*
      DO IGAS = 1, NGAS
        DO ISMST = 1, NSMST
          IF(NNSTSGP(ISMST,IGAS).GT.0) MAXVAL(IGAS) = ISMST
        END DO
        DO ISMST = NSMST,1,-1
          IF(NNSTSGP(ISMST,IGAS).GT.0) MINVAL(IGAS) = ISMST
        END DO
      END DO
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
 1000 CONTINUE
        IF(KFIRST .EQ. 1 ) THEN
          DO IGAS = 1, NGASL - 1
            ISMFGS(IGAS) = MINVAL(IGAS)
          END DO
        ELSE
*. Next distribution of symmetries in NGAS -1
         CALL NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
         IF(NONEW.NE.0) GOTO 1001
        END IF
        KFIRST = 0
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' next symmetry of NGASL-1 spaces '
          CALL IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
        END IF
*. Symmetry of NGASL -1 spaces given, symmetry of total space
        ISTSMM1 = 1
        DO IGAS = 1, NGASL -1
          CALL  SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
          ISTSMM1 = JSTSMM1
        END DO
*. required sym of SPACE NGASL
        CALL SYMCOM(2,1,ISTSMM1,ISMGSN,KSM)
C?      write(6,*) ' after  SYMCOM '
C?      write(6,*) ' ngasl istsmm1 ksm',ngasl,istsmm1,ksm
        ISMFGS(NGASL) = ISMGSN
*
        DO IGAS = NGASL+1,NGAS
          ISMFGS(IGAS) = 1
        END DO
        IF(NTEST.GE.200) THEN
          WRITE(6,*) ' Next symmetry distribution '
          CALL IWRTMA(ISMFGS,1,NGAS,1,NGAS)
        END IF
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IOBTP)
        CALL  SYMCOM(3,1,IOBSM,ISMFGS(IOBTP),IACSM)
        ISMFGS(IOBTP) = IACSM
        IOFF = 1
        MULT = 1
        DO IGAS = 1, NGAS
          IOFF = IOFF + (ISMFGS(IGAS)-1)*MULT
          MULT = MULT * NSMST
        END DO
        ISMFGS(IOBTP) = ISAVE
        IBSTRINI = INT(WORK(KLOFFI+IOFF-1))
C?      WRITE(6,*) ' IOFF IBSTRINI ', IOFF,IBSTRINI
*. Number of strings before active GAS space
        NSTB = 1
        DO IGAS = 1, IOBTP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings before active GAS space
        NSTA = 1
        DO IGAS =  IOBTP+1, NGAS
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group
C?      write(6,*) ' IACSM = ', IACSM
        NIAC  = NACIST(IACSM)
        IIAC =  IACIST(IACSM)
*
        NKAC = NNSTSGP(ISMFGS(IOBTP),IOBTP)
        IKAC = IISTSGP(ISMFGS(IOBTP),IOBTP)
*. I and K strings of given symmetry distribution
        NKSD = NSTB*NKAC*NSTA
C?      write(6,*) ' nstb nsta niac nkac ',
C?   &               nstb,nsta,niac,nkac
*. Obtain annihilation n mapping for all strings of this type
*
        NORBTS= NOBPTS(IOBTP,IOBSM)
*
        NKACT = NSTFGP(KACGRP)
C?      write(6,*) ' KACGRP ', KACGRP
        CALL ADSTN_GASSM(    NSTB,    NSTA,    IKAC,    IIAC,IBSTRINI,
     &                     KSTRBS,
     &                   IWORK(KSTSTM(KACGRP,1)),
     &                   IWORK(KSTSTM(KACGRP,2)),
     &                   IBORBSPS, IBORBSP,  NORBTS,    NKAC,   NKACT,
*
     &                       NIAC,   NKSTR, KBSTRIN,    NELB, NACGSOB,
     &                         I1,    XI1S,  SCLFAC)
        KSTRBS = KSTRBS + NKSD
        IF(NGASL-1.GT.0) GOTO 1000
 1001 CONTINUE
*
 9999 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADSTN_GAS '
        WRITE(6,*) ' ===================== '
        WRITE(6,*) ' Total number of K strings ', NKSTR
        IF(NKSTR.NE.0) THEN
          DO IORB = IBORBSPS,IBORBSPS + NORBTS  - 1
            IORBR = IORB-IBORBSPS +1
            WRITE(6,*) ' Info for orbital ', IORB
            WRITE(6,*) ' Excited strings and sign '
            CALL IWRTMA(  I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
            CALL WRTMAT(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
          END DO
        END IF
      END IF
*
* PAM Mars-2006: This flush moved outside of this subroutine
*      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'ADSTN ')

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(IEND)
        CALL Unused_integer(IFRST)
        CALL Unused_integer(KFRST)
        CALL Unused_integer(KACT)
      END IF
      END
