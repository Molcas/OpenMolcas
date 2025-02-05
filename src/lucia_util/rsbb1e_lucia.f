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
* Copyright (C) 1991,1997, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE RSBB1E_LUCIA(  ISCSM,  ISCTP,  ICCSM,  ICCTP,   IGRP,
     &                           NROW,   NGAS,   ISEL,   ICEL,     SB,
     &                             CB,  ADSXA, STSTSX, NOBPTS,   MAXI,
     &                           MAXK,   SSCR,   CSCR,     I1,   XI1S,
     &                             I2,   XI2S,      H,  NSMOB,  NSMST,
*
     &                          NSMSX,    MOC, MXSXST, IH2TRM, SCLFAC,
     &                        IUSE_PH, IPHGAS, NTESTG)
*
* SUBROUTINE RSBB1E_LUCIA --> 33
*
*
* One electron excitations on column strings
*. If IH2TRM .ne. 0 then the diagonal and one-electron
*  excitations arising from the two body operator is also
* included
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ISBCTP : Base for sigma column types
* ICCSM,ICCTP : Symmetry and type of C     columns
* ICBCTP : Base for C     column types
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active sets
* ISEL : Occupation in each active set for sigma block
* ICEL : Occupation in each active set for C     block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* MOC  : Use MOC method ( instead of N-1 resolution method )
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for one electron integrals
*
* Jeppe Olsen, Winter of 1991
*              IUSE_PH added winter of 97
*
      use Constants, only: One
      USE Para_Info, ONLY: MyRank, nProcs
      use lucia_data, only: MXPOBS,MXPNGAS,MXPTSOB
      IMPLICIT NONE
      INTEGER ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NGAS,MAXI,MAXK,NSMOB,
     &        NSMST,NSMSX, MOC,MXSXST,IH2TRM,IUSE_PH,NTESTG
      REAL*8  SCLFAC
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
      INTEGER IPHGAS(NGAS)
*.Specific Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      REAL*8 CB(*)
*.Output
      REAL*8 SB(*)
*.Scatch
      REAL*8 SSCR(*),CSCR(*),XI1S(*),H(*),XI2S(*)
      INTEGER I1(*),I2(*)
*.Local arrays ( assume MPNGAS = 16 ) !!!
      INTEGER ITP(16),JTP(16)
      INTEGER ISGRP(16),ICGRP(16)
*. For transposing integral block
      REAL*8 HSCR(MXPTSOB*MXPTSOB)
*
      INTEGER IJ_REO(2),IJ_DIM(2),IJ_SM(2),IJ_TP(2),IJ_AC(2)
* Type of single excitations that connects the two column strings
      INTEGER NTESTL,NTEST,NIPART,NIPARTSZ,IFRST,IJSM,IJTP,NSXTP,ITYP,
     &        JTYP,IXXX,ISM,JSM,KFRST,NIORB,NJORB,IDOCOMP,NKAEFF,NKASTR,
     &        KBOT,KTOP,KEND,LKABTC,IIPART,IBOT,ITOP,NIBTC,JJORB,ICGOFF,
     &        NIK,IIORB,ISBOFF,IEND,KACT
      REAL*8 SIGNIJ,SCLFACS,FACTORC,FACTORAB

C     MOC = 1
      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500)THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Information from RSBB1E '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' RSBB1E : MOC,IH2TRM,IUSE_PH ', MOC,IH2TRM,IUSE_PH
        WRITE(6,*) ' ISEL : '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' ICEL : '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF

*. Number of partitionings over column strings
CSVC: determine optimum number of partitions as the lowest multiple of
C     NPROCS that satisfies a block size smaller than MAXI:
      NIPART=0
      DO
        NIPART=NIPART+NPROCS
        NIPARTSZ=MAX(NROW-1,0)/NIPART+1
        IF (NIPARTSZ.LE.MAXI) EXIT
      END DO

*. Obtain groups
C     GET_SPGP_INF(ISPGP,ITP,IGRP)
      CALL GET_SPGP_INF(ICCTP,IGRP,ICGRP)
      CALL GET_SPGP_INF(ISCTP,IGRP,ISGRP)
*
      IFRST = 1
*. Types of single excitations that connect ISEL and ICEL
      CALL SXTYP2_GAS(    NSXTP,      ITP,      JTP,     NGAS,     ISEL,
     &                     ICEL,   IPHGAS)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.2000)
     &  write(6,*) ' ITYP JTYP ', ITYP,JTYP
*. Is this combination of types allowed
c        IJ_ACT = 1
c        IF(IJ_ACT.EQ.0) GOTO 900
*. Hvilken vej skal vi valge,
c        NOP = 2
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
c        ISCR(1) = ITYP
c        ISCR(2) = JTYP
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
*
c        ISCR(1) = ITYP
c        ISCR(2) = JTYP
c        IJ_TP(1) = ISCR(IJ_REO(1))
c        IJ_TP(2) = ISCR(IJ_REO(2))
*
        DO 800 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
*. New intermediate strings will be accessed so
          KFRST = 1
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.2000)
     &    write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
*. Reorder
*
*          ISCR(1) = ISM
*          ISCR(2) = JSM
*          IJ_SM(1) = ISCR(IJ_REO(1))
*          IJ_SM(2) = ISCR(IJ_REO(2))
*
*          ISCR(1) = NIORB
*          ISCR(2) = NJORB
*          IJ_DIM(1) = ISCR(IJ_REO(1))
*          IJ_DIM(2) = ISCR(IJ_REO(2))
*
      IF(IJ_REO(1).EQ.1) THEN
        IJ_SM(1)=ISM
        IJ_SM(2)=JSM
        IJ_DIM(1)=NIORB
        IJ_DIM(2)=NJORB
      ELSE
        IJ_SM(1)=JSM
        IJ_SM(2)=ISM
        IJ_DIM(1)=NJORB
        IJ_DIM(2)=NIORB
      END IF

          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*. Fetch integrals : For CI-transformations using RSBB1E
*. most of the blocks vanishes
*.Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
           IF(IJ_REO(1).EQ.1) THEN
*. obtain integrals h(j,i)
             CALL GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
             CALL TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
           ELSE
*. Obtain integrals h(i,j)
             CALL GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
           END IF
COLD       XNORM = INPROD(H,H,IJ_DIM(1)*IJ_DIM(2))
COLD       IF(XNORM.EQ.0) GOTO 800
          IF(MOC.EQ.0) THEN
*
*
* ======================================================================
*.                   Use N-1 resolution method
* ======================================================================
*
*
*. Obtain annihilation/creation maps for all K strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          SCLFACS = SIGNIJ*SCLFAC
          IF(NTEST.GE.1000)
     &    WRITE(6,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
          CALL ADAST_GAS(IJ_SM(2),IJ_TP(2),    NGAS,   ICGRP,   ICCSM,
     &                         I1,    XI1S,  NKASTR,    IEND,   IFRST,
     &                      KFRST,    KACT, SCLFACS,IJ_AC(1))
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          CALL ADAST_GAS(IJ_SM(1),IJ_TP(1),    NGAS,   ISGRP,   ISCSM,
     &                         I2,    XI2S,  NKASTR,    IEND,   IFRST,
     &                      KFRST,    KACT,     ONE,IJ_AC(1))
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(     I1,   XI1S,IJ_DIM(2),   I2, XI2S,
     &                        IJ_DIM(1),NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF
*. Loop over partitionings of the row strings
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
                IBOT = (IIPART-1)*MAXI+1
                ITOP = MIN(IBOT+MAXI-1,NROW)
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
*.Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
C               CALL GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
C               CALL TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
*. Problems when HOLE switches blocks around ?
C               CALL GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' RSBB1E H BLOCK '
                  CALL WRTMAT(H,IJ_DIM(2),IJ_DIM(1),IJ_DIM(2),IJ_DIM(1))
                END IF
*.Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
                NIK = NIBTC*LKABTC
                FACTORC = 0.0D0
                FACTORAB = 1.0D0
                IF(NTEST.GE.2000) THEN
                  WRITE(6,*) ' CSCR array,NIK X NJORB array '
                  CALL WRTMAT(CSCR,NIK,IJ_DIM(2),NIK,IJ_DIM(2))
                END IF
                CALL MATML7(   SSCR,   CSCR,      H,    NIK,IJ_DIM(1),
     &                        NIK,IJ_DIM(2),IJ_DIM(2),IJ_DIM(1),FACTORC,
     &                      FACTORAB,     0)
                IF(NTEST.GE.2000) THEN
                  WRITE(6,*) ' SSCR array,NIK X NIORB array '
                  CALL WRTMAT(SSCR,NIK,IJ_DIM(1),NIK,IJ_DIM(1))
                END IF
*.S(I,a+ K) =  S(I, a+ K) + sgn*Sscr(I,K,i)
                DO IIORB = 1,IJ_DIM(1)
                  ISBOFF = 1+(IIORB-1)*LKABTC*NIBTC
                  CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                         LKABTC,
     &                        I2(KBOT+(IIORB-1)*NKASTR),
     &                        XI2S(KBOT+(IIORB-1)*NKASTR))
                END DO
*
  701       CONTINUE
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
*. End of loop over I partitioninigs
          END IF
*.(End of algorithm switch)
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NSMSX)
        CALL Unused_integer(MXSXST)
      END IF
      END SUBROUTINE RSBB1E_LUCIA
