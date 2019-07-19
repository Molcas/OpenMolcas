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
* Copyright (C) 1991,1995,1998, Jeppe Olsen                            *
************************************************************************
      SUBROUTINE GSBBD1_LUCIA(   RHO1,  NACOB,  ISCSM,  ISCTP,  ICCSM,
     &                          ICCTP,   IGRP,   NROW,   NGAS,   ISEL,
     &                           ICEL,     SB,     CB,  ADSXA, SXSTST,
     &                         STSTSX,MXPNGAS, NOBPTS, IOBPTS,  ITSOB,
     &                           MAXI,   MAXK,   SSCR,   CSCR,     I1,
*
     &                           XI1S,     I2,   XI2S,      H,  NSMOB,
     &                          NSMST,  NSMSX, MXPOBS,  RHO1S, SCLFAC,
     &                        IUSE_PH, IPHGAS,IDOSRHO1, SRHO1,   IAB)
*
* SUBROUTINE GSBBD1_LUCIA --> 40
*
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen
* Particle-Hole version of Jan. 98
*
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1 : Updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "para_info.fh"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*), ITSOB(*)
      INTEGER IPHGAS(*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*), SRHO1(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(*),XI1S(*)
      DIMENSION I2(*),XI2S(*)
*.Local arrays ( assume MPNGAS = 16 ) !!!
      DIMENSION ITP(16*16),JTP(16*16)
*
      DIMENSION IJ_REO(2),IJ_DIM(2),IJ_SM(2),IJ_TP(2),IJ_AC(2)
      DIMENSION IJ_OFF(2)
      DIMENSION ISCR(2)
      DIMENSION ICGRP(16),ISGRP(16)
      DIMENSION H(*)
*. Add or subtract for spindensity
      IF(IAB.EQ.1) THEN
        XAB = 1.0D0
      ELSE
        XAB = -1.0D0
      END IF
*.Local arrays
      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================='
        WRITE(6,*) ' GSBBD1 in action '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' ISCSM, ICCSM = ', ISCSM, ICCSM
*
        WRITE(6,*) ' GSBBD1, sclfac ',SCLFAC
      END IF
*
      IFRST = 1
      JFRST = 1
*. Number of partitionings over column strings
CSVC: determine optimum number of partions as the lowest multiple of
C     NPROCS that satisfies a block size smaller than MAXI:
      NIPART=0
      DO
        NIPART=NIPART+NPROCS
        NIPARTSZ=MAX(NROW-1,0)/NIPART+1
        IF (NIPARTSZ.LE.MAXI) EXIT
      END DO

*. Groups defining supergroups
C          GET_SPGP_INF(ISPGP,ITP,IGRP)
      CALL GET_SPGP_INF(ICCTP,IGRP,ICGRP)
      CALL GET_SPGP_INF(ISCTP,IGRP,ISGRP)

* Type of single excitations that connects the two column strings
      CALL SXTYP2_GAS(    NSXTP,      ITP,      JTP,     NGAS,     ISEL,
     &                     ICEL,   IPHGAS)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(NTEST.GE.1000)
     &WRITE(6,*) ' ISCSM,ICCSM IJSM ', ISCSM,ICCSM,IJSM
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
*. Hvilken vej skal vi valge,
*. Mi pojdem drugim putem (C)
*. VV: the code below confuses Absoft compiler and was rewritten.
        NOP = 2
        IJ_AC(1) = 2
        IJ_AC(2) = 1
        IJ_TP(1) = ITYP
        IJ_TP(2) = JTYP
c        IF(IUSE_PH.EQ.1) THEN
c          CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TP,IJ_AC,IJ_REO,SIGNIJ)
c        ELSE
          IJ_REO(1) = 1
          IJ_REO(2) = 2
          SIGNIJ = 1.0D0
c        END IF
*
        if (IJ_REO(1).eq.1) THEN

*
        IJ_TP(1) = ITYP
        IJ_TP(2) = JTYP
        else
        IXXX = IJ_AC(1)
        IJ_AC(1) = IJ_AC(2)
        IJ_AC(2) = IXXX
*
        ISCR(1) = ITYP
        ISCR(2) = JTYP
        IJ_TP(1) = JTYP
        IJ_TP(2) = ITYP
        endif

c        ISCR(1) = IJ_AC(1)
c        ISCR(2) = IJ_AC(2)
c        IJ_AC(1) = ISCR(IJ_REO(1))
c        IJ_AC(2) = ISCR(IJ_REO(2))
*
c nasty code to avoid optimization
c        if(iscr(1).eq.-1000) print *,IJ_TP,IJ_REO
c        ISCR(1) = ITYP
c        ISCR(2) = JTYP
c        IJ_TP(1) = ISCR(IJ_REO(1))
c        IJ_TP(2) = ISCR(IJ_REO(2))

        DO 800 ISM = 1, NSMOB
*. new i and j so new intermediate strings
          KFRST = 1
*
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.1000) write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
          IBIORB = IOBPTS(ITYP,ISM)
          IBJORB = IOBPTS(JTYP,JSM)
*. Reorder
*          ISCR(1) = ISM
*          ISCR(2) = JSM
*          IJ_SM(1) = ISCR(IJ_REO(1))
*          IJ_SM(2) = ISCR(IJ_REO(2))
**
*          ISCR(1) = NIORB
*          ISCR(2) = NJORB
*          IJ_DIM(1) = ISCR(IJ_REO(1))
*          IJ_DIM(2) = ISCR(IJ_REO(2))
**
*          ISCR(1) = IBIORB
*          ISCR(2) = IBJORB
*          IJ_OFF(1) = ISCR(IJ_REO(1))
*          IJ_OFF(2) = ISCR(IJ_REO(2))
**
**
      IF(IJ_REO(1).EQ.1) THEN
        IJ_SM(1)=ISM
        IJ_SM(2)=JSM
        IJ_DIM(1)=NIORB
        IJ_DIM(2)=NJORB
        IJ_OFF(1)=IBIORB
        IJ_OFF(2)=IBJORB
      ELSE
        IJ_SM(1)=JSM
        IJ_SM(2)=ISM
        IJ_DIM(1)=NJORB
        IJ_DIM(2)=NIORB
        IJ_OFF(1)=IBJORB
        IJ_OFF(2)=IBIORB
      END IF

          IF(NTEST.GE.2000)
     &    WRITE(6,*) ' NIORB NJORB ', NIORB,NJORB
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          SCLFACS = SCLFAC*SIGNIJ
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
          CALL ADAST_GAS(IJ_SM(2),IJ_TP(2),    NGAS,   ICGRP,   ICCSM,
     &                         I1,    XI1S,  NKASTR,    IEND,   IFRST,
     &                      KFRST,    KACT, SCLFACS,IJ_AC(1))
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
      IF (NKASTR.EQ.0) GOTO 800
          ONE = 1.0D0
          CALL ADAST_GAS(IJ_SM(1),IJ_TP(1),    NGAS,   ISGRP,   ISCSM,
     &                         I2,    XI2S,  NKASTR,    IEND,   IFRST,
     &                      KFRST,    KACT,     ONE,IJ_AC(1))
      IF (NKASTR.EQ.0) GOTO 800
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(     I1,   XI1S,IJ_DIM(2),   I2, XI2S,
     &                        IJ_DIM(1),NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF
C         WRITE(6,*) ' NKAEFF NKASTR', NKAEFF,NKASTR

*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = MIN(KTOP + MAXK,NKAEFF)
              IF(KTOP.EQ.NKAEFF) THEN
                KEND = 1
              ELSE
                KEND = 0
              END IF
              LKABTC = KTOP - KBOT +1

*. This is the place to start over partitioning of I strings
              DO 701 IIPART = 1+MYRANK, NIPART, NPROCS
                IBOT = (IIPART-1)*NIPARTSZ+1
                ITOP = MIN(IBOT+NIPARTSZ-1,NROW)
                NIBTC = ITOP - IBOT + 1
                IF (NIBTC.LE.0) EXIT
* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
                DO JJORB = 1,IJ_DIM(2)
                  ICGOFF = 1 + (JJORB-1)*LKABTC*NIBTC
                  CALL MATCG(     CB,CSCR(ICGOFF),NROW,NIBTC, IBOT,
     &                        LKABTC,
     &                       I1(KBOT+(JJORB-1)*NKASTR),
     &                       XI1S(KBOT+(JJORB-1)*NKASTR) )
                END DO
* Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
                DO IIORB = 1,IJ_DIM(1)
*.Gather S Block
                  ISGOFF = 1 + (IIORB-1)*LKABTC*NIBTC
                  CALL MATCG(     SB,SSCR(ISGOFF),NROW,NIBTC, IBOT,
     &                        LKABTC,
     &                       I2(KBOT+(IIORB-1)*NKASTR),
     &                       XI2S(KBOT+(IIORB-1)*NKASTR) )
                END DO
*
                IF(NTEST.GE.1000) THEN
                 WRITE(6,*) ' CSCR and SSCR '
                 CALL WRTMAT(CSCR,IJ_DIM(2),NKI,IJ_DIM(2),NKI)
                 CALL WRTMAT(SSCR,IJ_DIM(1),NKI,IJ_DIM(1),NKI)
                END IF
*
*. And then the hard work
                NKI = LKABTC*NIBTC
                FACTORC = 0.0D0
                FACTORAB = 1.0D0
                CALL MATML7(  RHO1S,   SSCR,   CSCR,IJ_DIM(1),IJ_DIM(2),
     &                          NKI,IJ_DIM(1),  NKI,IJ_DIM(2),FACTORC,
     &                      FACTORAB,     1)
*
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' Block to one-body density '
                  CALL WRTMAT(RHO1S,IJ_DIM(1),IJ_DIM(2),
     &                              IJ_DIM(1),IJ_DIM(2))
                END IF
*. Scatter out to complete matrix
                DO JJORB = 1, IJ_DIM(2)
                  JORB = IJ_OFF(2)-1+JJORB
                  DO IIORB = 1, IJ_DIM(1)
                    IORB = IJ_OFF(1)-1+IIORB
                    RHO1((JORB-1)*NACOB+IORB) =
     &              RHO1((JORB-1)*NACOB+IORB) +
     &              RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)
                    IF(IDOSRHO1.EQ.1) THEN
                      SRHO1((JORB-1)*NACOB+IORB) =
     &                SRHO1((JORB-1)*NACOB+IORB) +
     &                XAB*RHO1S((JJORB-1)*IJ_DIM(1)+IIORB)
                    END IF
                  END DO
                END DO
*               /\ End of hard work

  701     CONTINUE
*. /\ end of this I partitioning
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
*. End of loop over I partitionings
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
C!    stop ' enforced stop in RSBBD1 '
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(SXSTST)
        CALL Unused_integer_array(ITSOB)
        CALL Unused_real_array(H)
        CALL Unused_integer(IUSE_PH)
      END IF
      END
