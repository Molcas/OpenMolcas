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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE STRTYP_GAS(IPRNT)
*
* Find groups of strings in each GA space
*
* Output : /GASSTR/
*
* Jeppe Olsen, Oct 1994
*
      IMPLICIT REAL*8(A-H,O-Z)
*
*
#include "mxpdim.fh"
#include "cgas.fh"
#include "lucinp.fh"
#include "orbinp.fh"
#include "cstate.fh"
#include "gasstr.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "crun.fh"
*. Local scratch
      DIMENSION IOCTYP(MXPSTT),IREOSPGP(MXPSTT),ISCR(MXPSTT)
      INTEGER IOCCLS(1),IBASSPC(1)
*
      NTESTL = 00
      NTEST = MAX(IPRNT,NTESTL)
*. As input NCISPC GAS spaces IGSOCCX are given.
* Obtain space that cantains all these as special cases
*
C?    WRITE(6,*) ' NCISPC ', NCISPC
      DO IGAS = 1, NGAS
       MINI = IGSOCCX(IGAS,1,1)
       MAXI = IGSOCCX(IGAS,2,1)
C?     WRITE(6,*) ' MINI and MAXI for ISPC = 1 ',MINI,MAXI
       DO ICISPC = 2, NCISPC
        MINI = MIN(MINI,IGSOCCX(IGAS,1,ICISPC))
        MAXI = MAX(MAXI,IGSOCCX(IGAS,2,ICISPC))
C?     WRITE(6,*) ' MINI and MAXI for ISPC =  ',ICISPC,MINI,MAXI
       END DO
       IGSOCC(IGAS,1) = MINI
       IGSOCC(IGAS,2) = MAXI
      END DO
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Compound GAS space : '
        WRITE(6,*) ' ====================='
        WRITE(6,'(A)')
        WRITE(6,'(A)') '         Min. occ    Max. occ '
        WRITE(6,'(A)') '         ========    ======== '
        DO IGAS = 1, NGAS
          WRITE(6,'(A,I2,3X,I3,9X,I3)')
     &    '   GAS',IGAS,IGSOCC(IGAS,1),IGSOCC(IGAS,2)
        END DO
      END IF

*
*. Find min and max number of elecs in each subspace
*
      DO IGAS = 1, NGAS
        IF(IGAS.EQ.1) THEN
          MNGSOC(IGAS) = IGSOCC(IGAS,1)
          MXGSOC(IGAS) = IGSOCC(IGAS,2)
        ELSE
          MXGSOC(IGAS) = IGSOCC(IGAS,2)-IGSOCC(IGAS-1,1)
          MNGSOC(IGAS) = MAX(0,IGSOCC(IGAS,1)-IGSOCC(IGAS-1,2))
        END IF
      END DO
*
*. Particle and hole spaces  :
*  Hole spaces are always more than half occupied
*
      IPHGASL = 0
      NPHGAS = 0
      DO IGAS = 1, NGAS
c        IF(IUSE_PH.EQ.1) THEN
c*. P/H separation for compound space
c          IF(  MNGSOC(IGAS).GT.NOBPT(IGAS)) THEN
c            IPHGAS(IGAS) = 2
c            IPHGASL = IGAS
c            NPHGAS = NPHGAS + 1
c          ELSE
c             IPHGAS(IGAS) = 1
c          END IF
c*. P/H separation for initial  space
c          IF(IGAS.EQ.1) THEN
c            MIN_OC1 = IGSOCCX(IGAS,1,1)
c          ELSE
c            MIN_OC1 = IGSOCCX(IGAS,1,1)-IGSOCCX(IGAS-1,2,1)
c          END IF
c          IF(MIN_OC1.GT.NOBPT(IGAS)) THEN
c            IPHGAS1(IGAS) = 2
c          ELSE
c             IPHGAS1(IGAS) = 1
c          END IF
c        ELSE IF(IUSE_PH.EQ.0) THEN
          IPHGAS(IGAS) = 1
          IPHGAS1(IGAS) = 1
c        END IF
      END DO
*. Large number of particle and hole orbitals
      MXTSOB_P = 0
      MXTSOB_H = 0
      DO IGAS = 1, NGAS
        IF(IPHGAS1(IGAS).EQ.1) THEN
          MXTSOB_P = MAX(MXTSOB_P,NOBPT(IGAS))
        ELSE
          MXTSOB_H = MAX(MXTSOB_H,NOBPT(IGAS))
        END IF
      END DO
      IF (NTEST .GT. 0)
     &   WRITE(6,*) ' MXTSOB_H, MXTSOB_P = ', MXTSOB_H, MXTSOB_P
*
C?    IF(IUSE_PH.EQ.1) THEN
C?    IPHGAS(1) = 2
C?    DO I = 1, 100
C?      WRITE(6,*) ' First space enforced to hole spaces'
C?    END DO
C?    END IF
*
*. In the following I assume that the hole spaces are the first NPGAS spaces
* ( only used when calculating min number of electrons, so it can be modified
*   easily )
*. Min number of electrons in hole spaces
C      IF (NTEST .GT. 0)
C     &   WRITE(6,*) ' IPHGASL, NPHGAS ',  IPHGASL, NPHGAS
C_Jesper      IF(IPHGASL.NE.NPHGAS) THEN
C_Jesper        WRITE(6,*) ' The hole spaces are not the first orbital spaces'
C_Jesper        STOP       ' The hole spaces are not the first orbital spaces'
C_Jesper      END IF
C      MNHL = IGSOCC(IPHGASL,1)
C     MNHL  = 0
C     DO IGAS = 1, NGAS
C       IF(IPHGAS(IGAS).EQ.2) THEN
C         MNHL = MNHL + MNGSOC(IGAS)
C       END IF
C     END DO
*
      IF (NTEST .GT. 0) THEN
        WRITE(6,*)' IPHGASL,NPHGAS:',IPHGASL,NPHGAS
        IF (IPHGASL.GT.0)
     &          WRITE(6,*)'IGSOCC(IPHGASL,1):',IGSOCC(IPHGASL,1)
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Min and Max occupation in each GAS space: '
        WRITE(6,'(A)') ' ========================================= '
        WRITE(6,*)
        DO IGAS = 1,  NGAS
          WRITE(6,'(A,I2,4X,2I3)')
     &    '  GAS',IGAS,MNGSOC(IGAS),MXGSOC(IGAS)
        END DO
*
        WRITE(6,*)' Particle(1) or hole(2) spaces (for compound space)'
        CALL IWRTMA(IPHGAS,1,NGAS,1,NGAS)
        WRITE(6,*)' Particle(1) or hole(2) spaces (for initial space)'
        CALL IWRTMA(IPHGAS1,1,NGAS,1,NGAS)
       END IF
*.
*. Occupation classes corresponding to largest CI space
*
      IBASSPC(1)=0
      CALL OCCLS(         1,    NOCCLS,    IOCCLS,    NACTEL,      NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),       0,    IBASSPC,   NOBPT)
      NMXOCCLS = NOCCLS
*
* Split into alpha and beta parts
*
*. Number of alpha. and beta electrons
*
      NAEL = (MS2 + NACTEL ) / 2
      NBEL = (NACTEL - MS2 ) / 2
*
      IF(NAEL + NBEL .NE. NACTEL ) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I4)')   MS2,NACTEL,NAEL,NBEL
        WRITE(6,*)
     &  ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT '
*        STOP ' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY '
        CALL SYSABENDMSG('lucia_util/strtyp_gas','Internal error',' ')
      END IF
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I6)')   MS2,NACTEL,NAEL,NBEL
      END IF
*
      IF(NAEL + NBEL .NE. NACTEL ) THEN
        WRITE(6,*) '  MS2 NACTEL NAEL NBEL '
        WRITE(6,'(5I4)')   MS2,NACTEL,NAEL,NBEL
        WRITE(6,*)
     &  ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT '
*          STOP ' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY '
        CALL SYSABENDMSG('lucia_util/strtyp_gas','Internal error',' ')
      END IF

*. Number of electrons to be subtracted or added
*
      MAXSUB = 2
      MAXADD = 2
*. electrons are only added for systems that atleast have halffilled
*. shells
      IGRP = 0
      MXAL = NAEL
      MNAL = NAEL
      MXBL = NBEL
      MNBL = NBEL
      NORBL = NTOOB
      DO IGAS = 1, NGAS
*. occupation constraints 1
       MXA1 = MIN(MXGSOC(IGAS),NOBPT(IGAS),MXAL)
       MXB1 = MIN(MXGSOC(IGAS),NOBPT(IGAS),MXBL)
       MNA1 = MAX(0,MNGSOC(IGAS)-MXA1)
       MNB1 = MAX(0,MNGSOC(IGAS)-MXB1)
*. Additional checks can be made here
       MXA = MXA1
       MXB = MXB1
       MNA = MNA1
       MNB = MNB1
*
       MXAL = MXAL - MNA
       MNAL = MAX(0,MNAL-MXA)
       MXBL = MXBL - MNB
       MNBL = MAX(0,MNBL-MXB)
*
       IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation numbers for IGAS = ', IGAS
        WRITE(6,*) ' MXAL MNAL MXBL MNBL ',MXAL,MNAL,MXBL,MNBL
        WRITE(6,*) ' MXA MNA MXB MNB ',MXA,MNA,MXB,MNB
       END IF
*
       MNAB = MIN(MNA,MNB)
       MXAB = MAX(MXA,MXB)
*. Additional holes only allowed in particle spaces
       IF(IPHGAS(IGAS).EQ.1) THEN
         MNAB = MAX(0,MNAB-MAXSUB)
       ELSE IF(IPHGAS(IGAS).EQ.2) THEN
         MNAB = MNAB
       END IF
C. For coupled cluster- could be refined ...
       MNAB = 0
*. Additional electrons allowed in hole spaces
       IF(IPHGAS(IGAS).EQ.2)  MXAB = MIN(MXAB + 2,NOBPT(IGAS))
*
       IF(NTEST.GE.100) WRITE(6,*) ' MNAB,MXAB',MNAB,MXAB
       NGPSTR(IGAS) = MXAB-MNAB+1
       IBGPSTR(IGAS) = IGRP + 1
       MNELFGP(IGAS) = MNAB
       MXELFGP(IGAS) = MXAB
*
       IADD = 0
       DO JGRP = IGRP+1,IGRP+NGPSTR(IGAS)
         IF(JGRP.GT.MXPSTT) THEN
           WRITE(6,*) ' Too many string groups '
           WRITE(6,*) ' Current limit ', MXPSTT
           WRITE(6,*) ' STOP : GASSTR, Too many string groups'
*                        STOP ' GASSTR, Too many string groups'
           CALL SYSABENDMSG('lucia_util/gasstr','Internal error',' ')
          END IF
*
         IADD = IADD + 1
         IEL = MNAB-1+IADD
         NELFGP(JGRP) = IEL
         IGSFGP(JGRP) = IGAS
         NSTFGP(JGRP) = IBION_LUCIA(NOBPT(IGAS),IEL)
       END DO
       IGRP = IGRP + NGPSTR(IGAS)
      END DO
      NGRP = IGRP
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,'(A)') ' Information about Groups of strings '
        WRITE(6,'(A)') ' =================================== '
        WRITE(6,*)
        WRITE(6,*) '     GAS  MNEL  MXEL IBGRP  NGRP'
        WRITE(6,*) '    ============================'
        DO IGAS = 1, NGAS
          WRITE(6,'(5(2X,I4))') IGAS,MNELFGP(IGAS),
     &          MXELFGP(IGAS),IBGPSTR(IGAS),NGPSTR(IGAS)
        END DO
        WRITE(6,'(A,I3)')
     &  ' Total number of groups generated ', NGRP
*
        WRITE(6,'(A)') ' Information about each string group '
        WRITE(6,'(A)') ' ===================================='
        WRITE(6,*)
        IITYPE = 0
        WRITE(6,'(A)') ' GROUP  GAS   NEL      NSTR '
        WRITE(6,'(A)') ' ==========================='
        DO IGRP = 1, NGRP
          IITYPE = IITYPE + 1
          WRITE(6,'(3(2X,I4),2X,I8)')
     &    IITYPE,IGSFGP(IGRP),NELFGP(IGRP),NSTFGP(IGRP)
        END DO
      END IF
*
*. Creation-annihilation connections between groups
*
      DO IGRP = 1, NGRP
        ISTAC(IGRP,1) = 0
        ISTAC(IGRP,2) = 0
        DO JGRP = 1, NGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)-1) ISTAC(IGRP,2) = JGRP
          IF(IGSFGP(IGRP).EQ.IGSFGP(JGRP).AND.
     &       NELFGP(IGRP).EQ.NELFGP(JGRP)+1) ISTAC(IGRP,1) = JGRP
        END DO
      END DO
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================================'
        WRITE(6,*) ' Annihilation / Creation connections'
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        CALL IWRTMA(ISTAC,NGRP,2,MXPSTT,2)
      END IF
*.
*
* Construct number of type ( combinations of groups ) with nael and nbel strings
*
*
* Type 1 : NAEL electrons
*      2 : NBEL ELECTRONS
*      3 : NAEL -1 ELECTRONS
*      4 : NBEL -1 ELECTRONS
*      5 : NAEL -2 ELECTRONS
*      6 : NBEL -2 ELECTRONS
*
      NSTTYP = 6
      NSTTP = 6
*
c      IF(IUSE_PH.EQ.1) THEN
c*. allow N+1,N+2 resolution string
c        NSTTYP = 10
c        NSTTP  = 10
c      END IF
*. alpha
      NELEC(1) = NAEL
      NELFTP(1) = NAEL
*. beta
      NELEC(2) = NBEL
      NELFTP(2) = NBEL
*. alpha -1
      NELEC(3) = NAEL-1
      NELFTP(3) = NAEL-1
*. beta  -1
      NELEC(4) = NBEL-1
      NELFTP(4) = NBEL-1
*. alpha -2
      NELEC(5) = NAEL-2
      NELFTP(5) = NAEL-2
*. beta  -2
      NELEC(6) = NBEL-2
      NELFTP(6) = NBEL-2
*
c      IF(IUSE_PH.EQ.1) THEN
c*. Alpha + 1
c        NELEC(7) = NAEL+1
c        NELFTP(7) = NAEL+1
c*. beta  + 1
c        NELEC(8) = NBEL+1
c        NELFTP(8) = NBEL+1
c*. Alpha + 2
c        NELEC(9) = NAEL+2
c        NELFTP(9) = NAEL+2
c*. beta  + 2
c        NELEC(10) = NBEL+2
c        NELFTP(10) = NBEL+2
c      END IF
*. Can easily be extended to relativistic case !!
      DO ITP = 1, NSTTYP
        NOCTYP(ITP) = 0
        NSPGPFTP(ITP) = 0
      END DO
*
* Loop over types, i.e.  given number of electrons
*
      IOFF = 1
      NABEL = NAEL + NBEL
      NSPGP_TOT = 0
      DO 2000 ITYP = 1, NSTTYP
*. Number of electrons in reference space ( alpha or beta )
        IF(MOD(ITYP,2).EQ.1) THEN
*. alpha type
          NELEC_REF = NELEC(1)
        ELSE
*. Beta type
          NELEC_REF = NELEC(2)
        END IF
*. If we are studying beta type, and number of alpha and beta
* electrons are identical, just refer to alpha
        IF(NAEL.EQ.NBEL.AND.MOD(ITYP,2).EQ.0) THEN
          IBSPGPFTP(ITYP) =  IBSPGPFTP(ITYP-1)
          NOCTYP(ITYP) =   NOCTYP(ITYP-1)
          NSPGPFTP(ITYP) =  NSPGPFTP(ITYP-1)
        ELSE
*
*
*. Number of electrons removed compared to reference
        IDEL = NELEC(ITYP) - NELEC_REF
C?      WRITE(6,*) '  GASSPC : ITYP IDEL ', ITYP,IDEL
*. Initial type of strings, relative to offset for given group
        DO IGAS = 1, NGAS
          IOCTYP(IGAS) = 1
        END DO
        NSPGP = 0
        IBSPGPFTP(ITYP) = IOFF
        IF(NELEC(ITYP).LT.0) THEN
          NOCTYP(ITYP) = 0
          NSPGPFTP(ITYP) =  0
          GOTO 2000
        END IF
*. Number of electrons in present type
*. Loop over  SUPER GROUPS with current nomenclature!
*. Temp max for loop
        MXLOOP = 10000
        IONE = 1
        NLOOP = 0
 1000   CONTINUE
*. Number of electrons in present supergroup
          NEL = 0
          DO IGAS = 1, NGAS
            NEL = NEL + NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
          END DO
*
          IF(NEL.GT.NELEC(ITYP)) THEN
*. If the number of electrons is to large find next number that
* can be correct.
* The following uses that within a given GAS space
* the number of elecs increases as the type number increases
*
*. First integer  that can be reduced
            IRED = 0
            DO IGAS = 1, NGAS
              IF(IOCTYP(IGAS).NE.1) THEN
                IRED = IGAS
                GOTO 888
              END IF
            END DO
  888       CONTINUE
            IF(IRED.EQ.NGAS) THEN
              NONEW = 1
            ELSE IF(IRED.LT.NGAS) THEN
              IOCTYP(IRED) = 1
*. Increase remanining part
              CALL NXTNUM2(
     &        IOCTYP(IRED+1),NGAS-IRED,IONE,NGPSTR(IRED+1),NONEW)
            END IF
            GOTO 2803
          END IF

          IF(NEL.EQ.NELEC(ITYP)) THEN
*. test 1 has been passed, check additional occupation constraints
*
           I_AM_OKAY = 1
*. Number of extra holes in hole spaces
CE         IF(IUSE_PH.EQ.1) THEN
CE           IDELP = 0
CE           IDELM = 0
CE           DO IGAS = 1, NGAS
CE             IF(IPHGAS(IGAS).EQ.2) THEN
CE               NELH =  NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
CE               IF(NELH.LT.MNGSOC(IGAS))
CE    &          IDELM = IDELM +MNGSOC(IGAS)-NELH
CE               IF(NELH.GT.MXGSOC(IGAS))
CE    &          IDELP = IDELP + NELH-MXGSOC(IGAS)
CE             END IF
CE           END DO
CE           IF(IDELM.GT.0.OR.IDELP.GT.MAX(0,IDEL)) THEN
CE             I_AM_OKAY = 0
CE             WRITE(6,*) ' P/H rejected supergroup '
CE             CALL IWRTMA(IOCTYP,1,NGAS,1,NGAS)
CE             WRITE(6,*) ' IDELM, IDELP ', IDELM, IDELP
CE           END IF
CE         END IF
*
*. Check from above
*
           DO IGAS = NGAS, 1, -1
*. Number of electrons when all electrons of AS IGAS have been added
             IF(IGAS.EQ.NGAS ) THEN
               IEL = MAX(NABEL,NABEL+IDEL)
             ELSE
               IEL = IEL-NELFGP(IOCTYP(IGAS+1)+IBGPSTR(IGAS+1)-1)
               IF(IEL.LT.MAX(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL))
     &         I_AM_OKAY = 0
             END IF
           END DO
*
* Check from below
*
           IEL = 0
           IOELMX = 0
           DO IGAS = 1, NGAS
             IEL = IEL + NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
             IOELMX = IOELMX+NOBPT(IGAS)
             IF(IEL+IOELMX.LT.MIN(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL))
     &       I_AM_OKAY = 0
           END DO
*
           IF(I_AM_OKAY.EQ.1) THEN
*. passed !!!
             NSPGP = NSPGP + 1
*. Copy supergroup to ISPGPFTP with absolute group numbers
             DO IGAS = 1, NGAS
               ISPGPFTP(IGAS,IOFF-1+NSPGP)
     &       = IOCTYP(IGAS)+IBGPSTR(IGAS)-1
             END DO
           END IF
*
          END IF
*. Next type of strings
          IONE = 1
          CALL NXTNUM2(IOCTYP,NGAS,IONE,NGPSTR,NONEW)
 2803   CONTINUE
        IF(NONEW.EQ.0) GOTO 1000
*. End of loop over possible supergroups, save information about current type
        IOFF = IOFF + NSPGP
        NOCTYP(ITYP) = NSPGP
        NSPGPFTP(ITYP) =  NSPGP
        NSPGP_TOT =  NSPGP_TOT +  NSPGP
      END IF
 2000 CONTINUE
      NTSPGP = NSPGP_TOT
*
      IF(NSPGP_TOT .GT. MXPSTT ) THEN
        WRITE(6,*) ' Too many super groups = ', NSPGP_TOT
        WRITE(6,*) ' Increase MXPSTT to this value'
        WRITE(6,*) ' See you later '
        WRITE(6,*)
        WRITE(6,*) ' STOP Increase MXPSTT '
*        STOP' Increase MXPSTT'
        CALL SYSABENDMSG('lucia_util/strtyp_gas','Internal error',' ')
      END IF
*
*. Reorder supergroups according to dimension
*
      DO ITYP= 1, NSTTP
        IBTYP = IBSPGPFTP(ITYP)
        NSPGP = NSPGPFTP(ITYP)
*.Dimension of supergroups
        DO ISPGP = 1, NSPGP
          IDIM = 1
          DO JGAS = 1, NGAS
            IDIM = IDIM * NSTFGP(ISPGPFTP(JGAS,ISPGP+IBTYP-1))
          END DO
          IOCTYP(ISPGP) = IDIM
        END DO
*. Reorder
C            ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
        CALL ORDINT(IOCTYP,ISCR,NSPGP,IREOSPGP,NTEST)
C?      WRITE(6,*) ' IREO array '
C?      CALL IWRTMA(IREOSPGP,1,NSPGP,1,NSPGP)
*. And reorder the definition of supergroups
        DO ISPGP = 1, NSPGP
          CALL ICOPVE(ISPGPFTP(1,ISPGP+IBTYP-1),NELFSPGP(1,ISPGP),NGAS)
        END DO
        DO ISPGP_N = 1, NSPGP
          ISPGP_O = IREOSPGP(ISPGP_N)
          CALL ICOPVE(NELFSPGP(1,ISPGP_O),ISPGPFTP(1,ISPGP_N+IBTYP-1),
     &                NGAS)
        END DO
*
      END DO
*     ^ End of loop over types
*
      IF(NTEST.GE.2) THEN
       WRITE(6,*) ' Total number of super groups ', NTSPGP
       WRITE(6,*) ' Number of alpha supergroups  ', NSPGPFTP(1)
       WRITE(6,*) ' Number of beta  supergroups  ', NSPGPFTP(2)
       WRITE(6,*)
       WRITE(6,*)
      END IF

*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' Information about types of strings'
        WRITE(6,*) ' =================================='
        WRITE(6,*)
        DO ITYP = 1, NSTTYP
          WRITE(6,*)
          WRITE(6,*) '      Type : ', ITYP
          WRITE(6,*) '      ==============='
          WRITE(6,*) '      Number of electrons  ',NELFTP(ITYP)
          WRITE(6,*) '      Number of super groups ', NSPGPFTP(ITYP)
          WRITE(6,*) '      Supergroups '
          DO ISPGP = 1, NSPGPFTP(ITYP)
            IOFF = IBSPGPFTP(ITYP)
            CALL IWRTMA(ISPGPFTP(1,IOFF-1+ISPGP),1,NGAS,1,NGAS)
          END DO
        END DO
      END IF
*
*
      RETURN
      END
