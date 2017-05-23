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
* Copyright (C) 1991,1994-1997, Jeppe Olsen                            *
*               2012, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE ADAST_GAS(   IOBSM,   IOBTP,   NIGRP,    IGRP, ISPGPSM,
     &                           I1,    XI1S,   NKSTR,    IEND,   IFRST,
     &                        KFRST,    KACT,  SCLFAC,     IAC)
*
*
* Obtain creation or annihilation mapping
*
* IAC = 2 : Creation map
* a+IORB !KSTR> = +/-!ISTR>
*
* IAC = 1 : Annihilation map
* a IORB !KSTR> = +/-!ISTR>
*
* for orbitals of symmetry IOBSM and type IOBTP
* and Istrings defined by the NIGRP groups IGRP and symmetry ISPGPSM
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
*               September 97 : annihilation mappings added
*                              I groups defined by IGRP
*
* Giovanni Li Manni, February 2012
* Smart Loop over symmetry distributions
* in order to make the code faster
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
#include "distsym.fh"
#include "loff.fh"
*. Input
      INTEGER IGRP(NIGRP)
*. Local scratch
      INTEGER ISMFGS(MXPNGAS)
      INTEGER MXVLI(MXPNGAS),MNVLI(MXPNGAS)
*      INTEGER MXVLK(MXPNGAS),MNVLK(MXPNGAS)
      INTEGER NNSTSGP(MXPNSMST,MXPNGAS)
      INTEGER IISTSGP(MXPNSMST,MXPNGAS)
      INTEGER KGRP(MXPNGAS)
      INTEGER IACIST(MXPNSMST), NACIST(MXPNSMST)
*. Temporary solution ( for once )
cSJS      PARAMETER(LOFFI=8*8*8*8*8)
cSJS * Declaring later so it can be in terms of NGAS and NIRREP
cSJS      PARAMETER(LOFFI=NGAS**NIRREP) !SJS
      DIMENSION IOFFI(LOFFI)
*
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
      CALL QENTER('ADAST ')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' ADAST_GAS in service '
        WRITE(6,*) ' ==================== '
        WRITE(6,*)
        WRITE(6,*)'GAS space (IOBTP), Symm of it (IOBSM) :',IOBTP,IOBSM
        WRITE(6,*) ' Supergroup in action : '
        WRITE(6,'(A,I3  )') ' Number of active spaces ', NIGRP
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        WRITE(6,*) '  Symmetry of supergroup : ', ISPGPSM
        WRITE(6,*) ' SCLFAC = ', SCLFAC
        IF(IAC.EQ.1) THEN
          WRITE(6,*) ' Annihilation mapping '
        ELSE IF(IAC.EQ.2) THEN
          WRITE(6,*) ' Creation mapping '
        ELSE
          WRITE(6,*) ' Unknown IAC parameter in ADAST ',IAC
        CALL SYSABENDMSG('lucia_util/adast_gas',
     &                    'Internal error',' ')
        END IF
      END IF
*. A few preparations
      NORBTS= NOBPTS(IOBTP,IOBSM)
      NORBT= NOBPT(IOBTP)
      IACGAS = IOBTP
*. First orbital of given GASpace
       IBORBSP = IELSUM(NOBPT,IOBTP-1)+1
*. First orbital of given GASpace and Symmetry
       IBORBSPS = IOBPTS(IOBTP,IOBSM)
      IF(NTEST.GE.100) THEN
        write(6,*) ' NORBTS per GAS and sym          :', NORBTS
        write(6,*) ' NORBT per GAS                   :', NORBT
        write(6,*) ' IACGAS GAS involved             :', IACGAS
        write(6,*) ' IBORBSP 1st orb per GAS         :', IBORBSP
        write(6,*) ' IBORBSPS 1st orb per GAS and sym:', IBORBSPS
      END IF
*
*====================================================
*. K strings : Supergroup, symmetry and distributions
*====================================================
      IF(IAC.EQ.1) THEN
       IDELTA = +1
      ELSE
       IDELTA = -1
      END IF
*. Is required mapping contained within current set of maps?
*. a:) Is active GASpace included in IGRP - must be
      IACGRP = 0
      DO JGRP = 1, NIGRP
       IF(IGSFGP(IGRP(JGRP)).EQ. IACGAS) IACGRP = JGRP
      END DO
*. Note : IACGRP is not the actual active group, it is the address of the
*         active group in IGRP
      IF(IACGRP.EQ.0) THEN
        WRITE(6,*) ' ADAST in problems '
        WRITE(6,*) ' Active GASpace not included in IGRP '
        WRITE(6,*) ' Active GASpace : ', IACGAS
        WRITE(6,'(A,20I3)') ' The active groups       ',
     &                      (IGRP(I),I=1,NIGRP)
        CALL SYSABENDMSG('lucia_util/adast_gas',
     &                    'Internal error',' ')
      END IF
*. b:) active group in K strings
      NIEL = NELFGP(IGRP(IACGRP))
      NKEL = NIEL + IDELTA
      IF(NTEST.GE.1000) WRITE(6,*) ' NIEL and NKEL ',NIEL,NKEL
      IF(NKEL.EQ.-1.OR.NKEL.EQ.NOBPT(IACGAS)+1) THEN
*. No strings with this number of elecs - be happy : No work
        NKSTR = 0
        KACT = 0
        KACGRP = 0
        GOTO 9999
      ELSE
*. Find group with NKEL electrons in IACGAS
        KACGRP = 0
        DO JGRP = IBGPSTR(IACGAS),IBGPSTR(IACGAS)+NGPSTR(IACGAS)-1
          IF(NELFGP(JGRP).EQ.NKEL) KACGRP = JGRP
        END DO
        IF(NTEST.GE.1000) WRITE(6,*) ' KACGRP = ',KACGRP
*. KACGRP is the Active group itself
        IF(KACGRP.EQ.0) THEN
          WRITE(6,*)' ADAST : cul de sac, active K group not found'
          WRITE(6,*)' GAS space and number of electrons ',
     &               IACGAS,NKEL
          CALL SYSABENDMSG('lucia_util/adast_gas',
     &                    'Internal error',' ')
        END IF
      END IF
*. Okay active K group was found and is nontrivial
      CALL SYMCOM(2,0,IOBSM,KSM,ISPGPSM)
*. The K supergroup
      CALL ICOPVE(IGRP,KGRP,NIGRP)
      KGRP(IACGRP) = KACGRP
*. Number of strings and symmetry distributions of K strings
      CALL NST_SPGRP( NIGRP,   KGRP,      KSM,iWORK(KNSTSGP(1)),NSMST,
     &                NKSTR,   NKDIST)
      IF(NTEST.GE.1000) WRITE(6,*) 'KSM,NKSTR,NKDIST:', KSM,NKSTR,NKDIST
      IF(NKSTR.EQ.0) GOTO 9999
*. Last active space in K strings and number of strings per group and sym
      NGASL = 1
      DO JGRP = 1, NIGRP
       IF(NELFGP(KGRP(JGRP)).GT.0) NGASL = JGRP
       CALL ICOPVE2(iWORK(KNSTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &              NNSTSGP(1,JGRP))
       CALL ICOPVE2(iWORK(KISTSGP(1)),(KGRP(JGRP)-1)*NSMST+1,NSMST,
     &              IISTSGP(1,JGRP))
      END DO
      IF(NTEST.GE.100) WRITE(6,*) 'NGASL', NGASL
*. MIN/MAX for Kstrings
*      CALL MINMAX_FOR_SYM_DIST(NIGRP,KGRP,MNVLK,MXVLK,NKDIST_TOT)
*      IF(NTEST.GE.100) THEN
*        write(6,*) 'MNVLK and MXVLK '
*        CALL IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
*        CALL IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
*      END IF
*. (NKDIST_TOT is number of distributions, all symmetries )
* ==============
*. I Strings
* ==============
*. Generate symmetry distributions of I strings with given symmetry
      CALL TS_SYM_PNT2(    IGRP,   NIGRP,   MXVLI,   MNVLI, ISPGPSM,
     &                    IOFFI,   LOFFI)
*. Offset and dimension for active group in I strings
      CALL ICOPVE2(iWORK(KISTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               IACIST)
*
      CALL ICOPVE2(iWORK(KNSTSGP(1)),(IGRP(IACGRP)-1)*NSMST+1,NSMST,
     &               NACIST)
*. Last entry in IGRP with a nonvanisking number of strings
      NIGASL = 1
      DO JGRP = 1, NIGRP
        IF(NELFGP(IGRP(JGRP)).GT.0) NIGASL = JGRP
      END DO
*. Number of electrons before active space
      NELB = 0
      DO JGRP = 1, IACGRP-1
        NELB = NELB + NELFGP(IGRP(JGRP))
      END DO
      IF(NTEST.GE.1000) WRITE(6,*) ' NELB = ', NELB
      ZERO =0.0D0
      IZERO = 0
      CALL ISETVC(I1,IZERO,NORBTS*NKSTR)
*
* Loop over symmetry distribtions of K strings
*
      KFIRST = 1
      KSTRBS = 1
      DO IGAS = 1, NIGRP
        ISMFGS(IGAS) = 1
      END DO
 1000 CONTINUE
*. Next distribution
*        CALL NEXT_SYM_DISTR(  NGASL,  MNVLK,  MXVLK, ISMFGS,    KSM,
*     &                       KFIRST,  NONEW)
         CALL NEXT_SYM_DISTR_NEW(NSMST,NGRP,KGRP,NIGRP,
     &                           ISMFGS,KSM,KFIRST,NONEW,
     &                  iWork(ISMDFGP),iWork(NACTSYM),iWork(ISMSCR))
        IF(NTEST.GE.1000) THEN
          write(6,*) ' Symmetry distribution '
          call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
        END IF
        IF(NONEW.EQ.1) GOTO 9999
        KFIRST = 0
*. Number of strings of this symmetry distribution
        NSTRIK = 1
        DO IGAS = 1, NGASL
          NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Offset for corresponding I strings
        ISAVE = ISMFGS(IACGRP)
        CALL  SYMCOM(3,1,IOBSM,ISMFGS(IACGRP),IACSM)
        ISMFGS(IACGRP) = IACSM
        IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,IOFFI,MXVLI,MNVLI)
        ISMFGS(IACGRP) = ISAVE
*. Number of strings before active GAS space
        NSTB = 1
C       DO IGAS = 1, IOBTP-1
        DO IGAS = 1, IACGRP-1
          NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number of strings After active GAS space
        NSTA = 1
C       DO IGAS =  IOBTP +1, NIGRP
        DO IGAS =  IACGRP+1, NIGRP
          NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
        END DO
*. Number and offset for active group
        NIAC  = NACIST(IACSM)
        IIAC =  IACIST(IACSM)
        NKAC = NNSTSGP(ISMFGS(IACGRP),IACGRP)
        IKAC = IISTSGP(ISMFGS(IACGRP),IACGRP)
*. I and K strings of given symmetry distribution
        NISD = NSTB*NIAC*NSTA
        NKSD = NSTB*NKAC*NSTA
        IF(NTEST.GE.1000) THEN
        write(6,*) ' nstb nsta niac nkac ',
     &               nstb,nsta,niac,nkac
        END IF
*. Obtain annihilation/creation mapping for all strings of this type
*. Are group mappings in expanded or compact form
        IF(IAC.EQ.1.AND.ISTAC(KACGRP,2).EQ.0) THEN
          IEC = 2
          LROW_IN = NKEL
        ELSE
          IEC = 1
          LROW_IN = NORBT
        END IF
        NKACT = NSTFGP(KACGRP)
*
        IF(NSTA*NSTB*NIAC*NKAC.NE.0)
     &  CALL ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,
     &                 iWORK(KSTSTM(KACGRP,1)),iWORK(KSTSTM(KACGRP,2)),
     &                 IBORBSPS,IBORBSP,NORBTS,NKAC,NKACT,NIAC,
     &                 NKSTR,KBSTRIN,NELB,NACGSOB,I1,XI1S,SCLFAC,IAC,
     &                 LROW_IN,IEC)
        KSTRBS = KSTRBS + NKSD
        GOTO 1000
*
 9999 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from ADAST_GAS '
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
      CALL QEXIT('ADAST ')
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(IEND)
        CALL Unused_integer(IFRST)
        CALL Unused_integer(KFRST)
      END IF
      END
