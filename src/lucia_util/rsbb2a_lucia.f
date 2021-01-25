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
* Copyright (C) 1991, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE RSBB2A_LUCIA(  ISCSM,  ISCTP,  ICCSM,  ICCTP,   IGRP,
     &                           NROW,  NSCOL,   NGAS,   ISOC,   ICOC,
     &                             SB,     CB,  ADSXA,  STSTDX,
     &                         SXDXSX,NOBPTS,
     &                           MAXI,   MAXK,   SSCR,   CSCR,     I1,
*
     &                           XI1S,          XINT,  NSMOB,
     &                          NSMST,    NSMDX,
     &                           MXSXST, MXSXBL,   IMOC, SCLFAC,
     &                         NTESTG, IPHGAS)
*
*
* SUBROUTINE RSBB2A_LUCIA --> 46
*
*
* two electron excitations on column strings
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NSCOL : Number of columns in S block
* ISEL1(3) : Number of electrons in RAS1(3) for S block
* ICEL1(3) : Number of electrons in RAS1(3) for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* DXSTST : Sym of sx,!st> => sym of sx !st>
* STSTDX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* SXDXSX : Symmetry of SX1,SX1*SX2 => symmetry of SX2
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals, strings,
*                           single excitations, double excitations
* MAXI   : Largest number of 'spectator strings' treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
*
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
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
      USE Para_Info, ONLY: MyRank, nProcs
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*)
      INTEGER IPHGAS(NGAS)
*
*.Input
      DIMENSION CB(*)
      INTEGER ISOC(NGAS),ICOC(NGAS)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),XINT(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*)
c      DIMENSION RIKSX(MXSXBL,4),RKJSX(MXSXBL,4)
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
C-jwk-cleanup      INTEGER I4_DIM(4),I4_SM(4)
      DIMENSION I4_TP(4),I4_REO(4),ISCR(4)
      INTEGER I4_AC(4)
*
      INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
*     DIMENSION SCR(MXPTSOB*MXPTSOB*MXPTSOB*MXPTSOB)
*

#include "oper.fh"
*
C-jwk-cleanup      DIMENSION IACAR(2),ITPAR(2)
      Call Allocate_Work(ipSCR,MXPTSOB**4)
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' RSBB2A speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' MXSXST MXSXBL = ', MXSXST,MXSXBL
        WRITE(6,*) ' RSBB2A : IMOC ', IMOC
        WRITE(6,*) ' ISOC and ICOC : '
        CALL IWRTMA(ISOC,1,NGAS,1,NGAS)
        CALL IWRTMA(ICOC,1,NGAS,1,NGAS)
*
*
      END IF
      IFRST = 1
      JFRST = 1
*
*.Types of DX that connects the two strings
*
      IDXSM = STSTDX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0)  GOTO 2001
*. Connecting double excitations
      CALL DXTYP2_GAS(   NDXTYP,      ITP,      JTP,      KTP,      LTP,
     &                     NGAS,     ISOC,     ICOC,   IPHGAS)
      DO 2000 IDXTYP = 1, NDXTYP
        ITYP = ITP(IDXTYP)
        JTYP = JTP(IDXTYP)
        KTYP = KTP(IDXTYP)
        LTYP = LTP(IDXTYP)
*. Is this combination of types allowed
c         IJKL_ACT = 1
c         IF(IJKL_ACT.EQ.0) GOTO 2000
*
C?      write(6,*) ' test inserted in RSBB2AN'
C?      NPTOT = 0
C?      IF(ITYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(JTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(KTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(LTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(NPTOT.EQ.3) GOTO 2000
*
        ITYP_ORIG = ITYP
        JTYP_ORIG = JTYP
        KTYP_ORIG = KTYP
        LTYP_ORIG = LTYP
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' ITYP_ORIG, JTYP_ORIG, KTYP_ORIG, LTYP_ORIG',
     &                 ITYP_ORIG, JTYP_ORIG, KTYP_ORIG, LTYP_ORIG
        END IF
        NIJKL1 = 0
        IF(ITYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(JTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(KTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(LTYP.EQ.1) NIJKL1 = NIJKL1+1
*. Optimal ordering of operators
        I4_AC(1) = 2
        I4_AC(2) = 2
        I4_AC(3) = 1
        I4_AC(4) = 1
        I4_TP(1) = ITYP
        I4_TP(2) = KTYP
        I4_TP(3) = LTYP
        I4_TP(4) = JTYP
c        IF(IUSE_PH.EQ.1) THEN
c          NOP = 4
c          CALL ALG_ROUTERX(ISOC,JSOC,NOP,I4_TP,I4_AC,I4_REO,SIGN4)
c        ELSE
          DO IJKL = 1, 4
            I4_REO(IJKL) = IJKL
          END DO
c        END IF
*. Type of operators : TP and AC
        DO IJKL = 1, 4
C         ISCR( I4_REO(IJKL) ) = I4_TP(IJKL)
          ISCR(IJKL) = I4_TP( I4_REO(IJKL) )
        END DO
        DO IJKL = 1, 4
          I4_TP(IJKL) = ISCR(IJKL)
        END DO
        ITYP = I4_TP(1)
        KTYP = I4_TP(2)
        LTYP = I4_TP(3)
        JTYP = I4_TP(4)
        DO IJKL = 1, 4
          ISCR(IJKL) = I4_AC( I4_REO(IJKL) )
        END DO
        DO IJKL = 1, 4
          I4_AC(IJKL) = ISCR(IJKL)
        END DO
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' I4_AC, IT_TP  defined '
          WRITE(6,*) ' I4_AC, I4_TP '
          CALL IWRTMA(I4_AC,1,4,1,4)
          CALL IWRTMA(I4_TP,1,4,1,4)
        END IF
*
CSVC: determine optimum number of partions as the lowest multiple of
C     NPROCS that satisfies a block size smaller than MAXI:
        NPART=0
        DO
          NPART=NPART+NPROCS
          NPARTSZ=MAX(NROW-1,0)/NPART+1
          IF (NPARTSZ.LE.MAXI) EXIT
        END DO
*
*      ==================================
        IF(I4_AC(1).EQ.I4_AC(2) ) THEN
*      ==================================
*
*. a+ a+ a a or a a a+ a+
*. Largest possible number of orbital pairs
          MI = 0
          MJ = 0
          MK = 0
          ML = 0
          DO IOBSM = 1, NSMST
            MI = MAX(MI,NOBPTS(ITYP,IOBSM))
            MJ = MAX(MJ,NOBPTS(JTYP,IOBSM))
            MK = MAX(MK,NOBPTS(KTYP,IOBSM))
            ML = MAX(ML,NOBPTS(LTYP,IOBSM))
          END DO
          MXPAIR = MAX(MI*MK,MJ*ML)
*. Largest posssible
*. Symmetry of allowed Double excitation,loop over excitations
          DO 1950 IKOBSM = 1, NSMOB
            JLOBSM = SXDXSX(IKOBSM,IDXSM)
            IF(NTEST.GE.100) WRITE(6,*) ' IKOBSM,JLOBSM', IKOBSM,JLOBSM
            IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
            KFRST = 1
*
*. Number of batchs of symmetry pairs IK
*
            LENGTH = 0
            NIKBT = 0
            NBLK = 0
            NBLKT = 0
            DO ISM = 1, NSMOB
              KSM = ADSXA(ISM,IKOBSM)
              NI = NOBPTS(ITYP,ISM)
              NK = NOBPTS(KTYP,KSM)
              IF(NTEST.GE.100) write(6,*) ' NI, NK' , NI,NK
*
              IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                NIK = NI*(NI+1)/2
              ELSE IF(ITYP.GT.KTYP.OR.(ITYP.EQ.KTYP.AND.ISM.GT.KSM))THEN
                NIK = NI*NK
              ELSE
                NIK = 0
              END IF
              IF(NIK.NE.0) THEN
                NBLKT = NBLKT + 1
                IF(LENGTH+NIK .GT. MXPAIR) THEN
*. The present batch is complete
                  NIKBT = NIKBT+1
                  IKBT(1,NIKBT) = NBLKT - NBLK
                  IKBT(2,NIKBT) = NBLK
                  IKBT(3,NIKBT) = LENGTH
                  LENGTH = 0
                  NBLK = 0
                END IF
                NBLK = NBLK + 1
                LENGTH = LENGTH + NIK
                IKSMBT(1,NBLKT) = ISM
                IKSMBT(2,NBLKT) = KSM
              END IF
            END DO
*. The last batch
            IF(NBLK.NE.0) THEN
              NIKBT = NIKBT+1
              IKBT(1,NIKBT) = NBLKT - NBLK + 1
              IKBT(2,NIKBT) = NBLK
              IKBT(3,NIKBT) = LENGTH
            END IF

*.
            IF(NTEST.GE.2000) THEN
              WRITE(6,*) ' ITYP, KTYP, IKOBSM,  NIKBT = ',
     &                     ITYP, KTYP, IKOBSM,  NIKBT
              WRITE(6,*) ' IKBT : Offset, number, length '
              DO JIKBT = 1, NIKBT
                WRITE(6,'(3i3)') (IKBT(II,JIKBT), II = 1, 3)
              END DO
              WRITE(6,*) ' IKSMBT '
              CALL IWRTMA(IKSMBT,2,NBLKT,2,8)
            END IF
*
*. Number of batchs of symmetry pairs JL
*
            LENGTH = 0
            NJLBT = 0
            NBLK = 0
            NBLKT = 0
            DO JSM = 1, NSMOB
              LSM = ADSXA(JSM,JLOBSM)
              NJ = NOBPTS(JTYP,JSM)
              NL = NOBPTS(LTYP,LSM)
*
              IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                NJL = NJ*(NJ+1)/2
              ELSE IF(JTYP.GT.LTYP.OR.(JTYP.EQ.LTYP.AND.JSM.GT.LSM))THEN
                NJL = NJ*NL
              ELSE
                NJL = 0
              END IF
              IF(NJL.NE.0) THEN
                NBLKT = NBLKT + 1
                IF(LENGTH+NJL .GT. MXPAIR) THEN
*. The present batch is complete
                  NJLBT = NJLBT+1
                  JLBT(1,NJLBT) = NBLKT - NBLK
                  JLBT(2,NJLBT) = NBLK
                  JLBT(3,NJLBT) = LENGTH
                  LENGTH = 0
                  NBLK = 0
                END IF
                NBLK = NBLK + 1
                LENGTH = LENGTH + NJL
                JLSMBT(1,NBLKT) = JSM
                JLSMBT(2,NBLKT) = LSM
              END IF
            END DO
*. The last batch
            IF(NBLK.NE.0) THEN
              NJLBT = NJLBT+1
              JLBT(1,NJLBT) = NBLKT - NBLK + 1
              JLBT(2,NJLBT) = NBLK
              JLBT(3,NJLBT) = LENGTH
            END IF
*.
            IF(NTEST.GE.2000) THEN
              WRITE(6,*) ' JTYP, LTYP, JLOBSM,  NJLBT = ',
     &                     JTYP, LTYP, JLOBSM,  NJLBT
              WRITE(6,*) ' JLBT : Offset, number, length '
              DO JJLBT = 1, NJLBT
                WRITE(6,'(3i3)') (JLBT(II,JJLBT), II = 1, 3)
              END DO
              WRITE(6,*) ' JLSMBT '
              CALL IWRTMA(JLSMBT,2,NBLKT,2,8)
            END IF
*
*. Loop over batches of IK strings
            DO 1940 IKBTC = 1, NIKBT
              IF(NTEST.GE.1000) WRITE(6,*) ' IKBTC = ', IKBTC
*. Loop over batches of JL strings
              DO 1930 JLBTC = 1, NJLBT
                IFIRST = 1
*. Loop over batches of I strings
                DO 1801 IIPART = 1+MYRANK, NPART, NPROCS
                  IBOT = 1+(IIPART-1)*NPARTSZ
                  ITOP = MIN(IBOT+NPARTSZ-1,NROW)
                  NIBTC = ITOP-IBOT+1
                  IF (NIBTC.LE.0) EXIT
*.Loop over batches of intermediate strings
                  KBOT = 1- MAXK
                  KTOP = 0
 1800             CONTINUE
                    KBOT = KBOT + MAXK
                    KTOP = KTOP + MAXK
*
                    IONE = 1
                    JLBOFF = 1
                    NJLT = JLBT(3,JLBTC)
                    DO JLPAIR = 1, JLBT(2,JLBTC)
                      JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                      LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
                      NJ = NOBPTS(JTYP,JSM)
                      NL = NOBPTS(LTYP,LSM)
                      IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                        NJL = NJ*(NJ+1)/2
                        JLSM = 1
                      ELSE
                        NJL = NJ * NL
                        JLSM = 0
                      END IF

*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
                      II12 = 1
                      K12 = 1
                      IONE = 1
C?       write(6,*) ' Before ADAADAST '
*. Creation / annihilation maps , conjugated of above
                      IF(I4_AC(4).EQ.1) THEN
                        JAC = 2
                      ELSE
                        JAC = 1
                      END IF
                      IF(I4_AC(3).EQ.1) THEN
                        LAC = 2
                      ELSE
                        LAC = 1
                      END IF
                      CALL ADAADAST_GAS( IONE,  JSM, JTYP,   NJ,  JAC,
     &                                   IONE,  LSM, LTYP,   NL,  LAC,
     &                                  ICCTP,ICCSM, IGRP, KBOT, KTOP,
     &                                     I1, XI1S, MAXK,NKBTC, KEND,
     &                                  JFRST,KFRST, II12,  K12,SCLFAC)
*
                      JFRST = 0
                      KFRST = 0
*
                      IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                      J = 0
                      L = 1
*
                      DO  IJL = 1, NJL
                        CALL NXTIJ(     J,     L,    NJ,    NL,  JLSM,
     &                              NONEW)
                        I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                        JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                        IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                          ZERO = 0.0D0
                          CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                        ELSE
                          CALL MATCG(    CB,CSCR(JLOFF),NROW,NIBTC,IBOT,
     &                                NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                        END IF
                      END DO
*
                      JLBOFF = JLBOFF + NJL
                    END DO
*
*. ( End of loop over jlpair in batch )
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals xint(ik,jl) = (ij!kl)-(il!kj)
                    IF(IFIRST.EQ.1) THEN
                      IXCHNG = 1
* Obtain integrals in ik batch
                      NIKT = IKBT(3,IKBTC)
                      NJLT = JLBT(3,JLBTC)
                      JLOFF = 1
                      DO JLPAIR = 1, JLBT(2,JLBTC)
                      IKOFF = 1
                      DO IKPAIR = 1, IKBT(2,IKBTC)
*
                        ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                        KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                        JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                        LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
*
                        IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                          IKSM = 1
                          NIK =
     &                    NOBPTS(ITYP,ISM)*(NOBPTS(ITYP,ISM)+1)/2
                        ELSE
                          IKSM = 0
                          NIK =
     &                    NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
                        END IF
*
                        IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                          JLSM = 1
                          NJL =
     &                    NOBPTS(JTYP,JSM)*(NOBPTS(JTYP,JSM)+1)/2
                        ELSE
                          JLSM = 0
                          NJL =
     &                    NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM)
                        END IF
* ================================================================
*. Required form of integrals : Coulomb - Exchange of just Coulomb
* ================================================================
                        ICOUL = 0
*. Use coulomb - exchange
                        IXCHNG = 1
*. fetch integrals
*. Full conjugation symmetry, do do not worry
*                       CALL GETINT(   SCR,  ITYP,   ISM,  JTYP,   JSM,
                        CALL GETINT(Work(ipSCR),
     &                                       ITYP,   ISM,  JTYP,   JSM,
     &                                KTYP,   KSM,  LTYP,   LSM,IXCHNG,
     &                                IKSM,  JLSM, ICOUL)
*                       ^ End if similarity transformed Hamiltonian is
*                         used
                        DO JL = 1, NJL
*                         CALL COPVEC(SCR((JL-1)*NIK+1),
                          CALL COPVEC(Work((JL-1)*NIK+ipSCR),
     &                         XINT((JLOFF-1+JL-1)*NIKT+IKOFF),NIK)
                        END DO
                        IKOFF = IKOFF + NIK
                      END DO
                      JLOFF = JLOFF + NJL
                      END DO
                    END IF
*                   ^ End if integrals should be fetched
                    IFIRST = 0
*.and now ,to the work
                    LIKB = NIBTC*NKBTC
                    IF(NTEST.GE.3000) THEN
                     WRITE(6,*) ' Integral block '
                     CALL WRTMAT(XINT,NIKT,NJLT,NIKT,NJLT)
                    END IF
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' CSCR matrix '
                      CALL WRTMAT(CSCR,LIKB,NJLT,LIKB,NJLT)
                    END IF
*
C?                  MXACIJO = MXACIJ
c                    MXACIJ = MAX(MXACIJ,LIKB*NJLT,LIKB*NIKT)
C?                  IF(MXACIJ.GT.MXACIJO) THEN
C?                    write(6,*) ' New max MXACIJ = ', MXACIJ
C?                    write(6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
C?                    WRITE(6,*) ' ITYP,JTYP,KTYP,LTYP',
C?   &                             ITYP,JTYP,KTYP,LTYP
C?                    WRITE(6,*)'NIJT, NJLT, NIBTC NKBTC',
C?   &                           NIJT, NJLT,NIBTC,NKBTC
C?                  END IF
*
                    FACTORC = 0.0D0
                    FACTORAB = 1.0D0
                    CALL MATML7(   SSCR,   CSCR,   XINT,   LIKB,   NIKT,
     &                             LIKB,   NJLT,   NIKT,   NJLT,FACTORC,
     &                          FACTORAB,     2)
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' SSCR matrix '
                      CALL WRTMAT(SSCR,LIKB,NIKT,LIKB,NIKT)
                    END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                    II12 = 2
*
                    IONE = 1
                    IKBOFF = 1
                    DO IKPAIR = 1, IKBT(2,IKBTC)
                      ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                      KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                      NI = NOBPTS(ITYP,ISM)
                      NK = NOBPTS(KTYP,KSM)
                      IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                        NIK = NI*(NI+1)/2
                        IKSM = 1
                      ELSE
                        NIK = NI * NK
                        IKSM = 0
                      END IF
                      IF(IFRST.EQ.1) KFRST = 1
                      ONE = 1.0D0
*
                      IAC = I4_AC(1)
                      KAC = I4_AC(2)
*
                      CALL ADAADAST_GAS( IONE,  ISM, ITYP,   NI,  IAC,
     &                                   IONE,  KSM, KTYP,   NK,  KAC,
     &                                  ISCTP,ISCSM, IGRP, KBOT, KTOP,
     &                                     I1, XI1S, MAXK,NKBTC, KEND,
     &                                  IFRST,
*
     &                                  KFRST,
     &                                   II12,
     &                                    K12,
     &                                  ONE         )
*
                      IFRST = 0
                      KFRST = 0
*
                      I = 0
                      K = 1
                      DO IK = 1, NIK
                        CALL NXTIJ(     I,     K,    NI,    NK,  IKSM,
     &                              NONEW)
                        IKOFF = (K-1)*NI + I
                        ISBOFF = 1+(IKBOFF-1+IK-1)*NIBTC*NKBTC
                        IF(IKSM.EQ.1.AND.I.EQ.k) THEN
* a+ i a+i gives trivially zero
                        ELSE
                         CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                                NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                        END IF
                      END DO
                      IKBOFF = IKBOFF + NIK
*
                    END DO
*                   ^ End of loop over IKPAIRS in batch
*
                  IF(KEND.EQ.0) GOTO 1800
*.                ^ End of loop over partitionings of resolution strings
 1801           CONTINUE
*               ^ End of loop over partionings of I strings
 1930         CONTINUE
*             ^ End of loop over batches of JL
 1940       CONTINUE
*           ^ End of loop over batches of IK
 1950     CONTINUE
*         ^ End of loop over IKOBSM
*
*      ==============================================
        ELSE IF(.NOT.( I4_AC(1).EQ. I4_AC(2)) ) THEN
*      ==============================================
*
*
* Three types of operators :
* a+ a  a+ a
* a+ a  a  a+
* a  a+ a+ a
*
* The first end up with
* -a+ i ak a+l aj X2(ik,jl)
*
* Number two and three end up with
* -a i a k a l aj XC(ik,jl)  ( In coulomb form)
*
          JLSM = 0
          IKSM = 0
*. Symmetry of allowed Double excitation,loop over excitations
          DO 2950 IKOBSM = 1, NSMOB
            JLOBSM = SXDXSX(IKOBSM,IDXSM)
            IF(JLOBSM.EQ.0) GOTO 2950
*. types + symmetries defined => K strings are defined
            KFRST = 1
            DO ISM = 1, NSMOB
              KSM = ADSXA(ISM,IKOBSM)
              DO JSM = 1, NSMOB
                LSM = ADSXA(JSM,JLOBSM)
                IF(NTEST.GE.2000) WRITE(6,*) ' ISM KSM LSM JSM',
     &          ISM,KSM,LSM,JSM
                ISCR(I4_REO(1)) = ISM
                ISCR(I4_REO(2)) = KSM
                ISCR(I4_REO(3)) = LSM
                ISCR(I4_REO(4)) = JSM
*
                ISM_ORIG = ISCR(1)
                KSM_ORIG = ISCR(2)
                LSM_ORIG = ISCR(3)
                JSM_ORIG = ISCR(4)
*
C           DO ISM_ORIG = 1, NSMOB
C             KSM_ORIG = ADSXA(ISM_ORIG,IKOBSM)
C             DO JSM_ORIG = 1, NSMOB
C               LSM_ORIG = ADSXA(JSM_ORIG,JLOBSM)
*
C               ISCR(1) = ISM_ORIG
C               ISCR(2) = KSM_ORIG
C               ISCR(3) = LSM_ORIG
C               ISCR(4) = JSM_ORIG
*
C               ISM = ISCR(I4_REO(1))
C               KSM = ISCR(I4_REO(2))
C               LSM = ISCR(I4_REO(3))
C               JSM = ISCR(I4_REO(4))
*
                NI = NOBPTS(ITYP,ISM)
                NJ = NOBPTS(JTYP,JSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
                NIK = NI*NK
                NJL = NJ*NL
                IF(NIK.EQ.0.OR.NJL.EQ.0) GOTO 2803
*
                ITPSM_ORIG = (ITYP_ORIG-1)*NSMOB + ISM_ORIG
                JTPSM_ORIG = (JTYP_ORIG-1)*NSMOB + JSM_ORIG
                KTPSM_ORIG = (KTYP_ORIG-1)*NSMOB + KSM_ORIG
                LTPSM_ORIG = (LTYP_ORIG-1)*NSMOB + LSM_ORIG
*
                IF(ITPSM_ORIG.GE.KTPSM_ORIG.AND.
     &             JTPSM_ORIG.GE.LTPSM_ORIG) THEN
*
                IFIRST = 1
*. Loop over batches of I strings
                DO 2801 IIPART = 1+MYRANK, NPART, NPROCS
                  IBOT = 1+(IIPART-1)*NPARTSZ
                  ITOP = MIN(IBOT+NPARTSZ-1,NROW)
                  NIBTC = ITOP-IBOT+1
                  IF (NIBTC.LE.0) EXIT
*.Loop over batches of intermediate strings
                  KBOT = 1- MAXK
                  KTOP = 0
 2800             CONTINUE
                    KBOT = KBOT + MAXK
                    KTOP = KTOP + MAXK
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
                    II12 = 1
                    K12 = 1
                    IONE = 1
*. Creation / annihilation maps , conjugated of above
                    IF(I4_AC(4).EQ.1) THEN
                      JAC = 2
                    ELSE
                      JAC = 1
                    END IF
                    IF(I4_AC(3).EQ.1) THEN
                      LAC = 2
                    ELSE
                      LAC = 1
                    END IF
C                   KFRST = 1
                    CALL ADAADAST_GAS( IONE,  JSM, JTYP,   NJ,  JAC,
     &                                 IONE,  LSM, LTYP,   NL,  LAC,
     &                                ICCTP,ICCSM, IGRP, KBOT, KTOP,
     &                                   I1, XI1S, MAXK,NKBTC, KEND,
     &                                JFRST,KFRST, II12,  K12,SCLFAC)
*
                    JFRST = 0
                    KFRST = 0
*
                    IF(NKBTC.EQ.0) GOTO 2801
*. Loop over jl in TS classes and gather
                    J = 0
                    L = 1
                    DO  IJL = 1, NJL
                      CALL NXTIJ(     J,     L,    NJ,    NL,  JLSM,
     &                            NONEW)
                      I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                      JLOFF = (IJL-1)*NKBTC*NIBTC+1
                      CALL MATCG(    CB,CSCR(JLOFF),NROW,NIBTC,IBOT,
     &                            NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                    END DO
*
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals as xint(ik,jl) = (ij!kl)-(il!kj)
                    IKSM = 0
                    JLSM = 0
                    IF(IFIRST.EQ.1) THEN
                      IF(I4_AC(1).EQ.I4_AC(3)) THEN
* a+ a a+ a
                        ICOUL = 2
                      ELSE
* a+ a a a+ or a+ a a a+
                        ICOUL = 1
                      END IF
*. Use coulomb - exchange or just coulomb integrals ?
                      IF(ITPSM_ORIG.EQ.KTPSM_ORIG
     &                .AND.JTPSM_ORIG.EQ.LTPSM_ORIG)THEN
*. No use of exchange
                        IXCHNG = 0
                        FACX = -0.5D0
                      ELSE IF(ITPSM_ORIG.NE.KTPSM_ORIG
     &                .OR.JTPSM_ORIG.NE.LTPSM_ORIG) THEN
*. Exchange used, combines two terms
                        IXCHNG = 1
                        FACX = -0.5D0
                      END IF
                      IF(ITPSM_ORIG.NE.KTPSM_ORIG
     &                .AND.JTPSM_ORIG.NE.LTPSM_ORIG)THEN
*. Exchange used, combines four terms
                        IXCHNG = 1
                        FACX = -1.0D0
                      END IF
           IF( NTEST.GE.1000) WRITE(6,*)
     &   ' ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX',
     &     ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX
*. fetch integrals
* we want the operator in the form a+i ak a+l aj ((ij!lk)-(ik!lj))
                      IF(ICOUL.EQ.2) THEN
*. Obtain X2(ik,lj) = (ij!lk)
                      CALL GETINT(  XINT,  ITYP,   ISM,  JTYP,   JSM,
     &                              LTYP,   LSM,  KTYP,   KSM,IXCHNG,
     &                              IKSM,  JLSM, ICOUL)
                      ELSE IF (ICOUL.EQ.1) THEN
                        CALL GETINT(  XINT,  ITYP,   ISM,  KTYP,   KSM,
     &                                JTYP,   JSM,  LTYP,   LSM,IXCHNG,
     &                                IKSM,  JLSM, ICOUL)
                      END IF
*
                    END IF
*                   ^ End if integrals should be fetched
                    IFIRST = 0
*.and now ,to the work
                    LIKB = NIBTC*NKBTC
                    IF(NTEST.GE.3000) THEN
                     WRITE(6,*) ' Integral block '
                     CALL WRTMAT(XINT,NIK,NJL,NIK,NJL)
                    END IF
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' CSCR matrix '
                      CALL WRTMAT(CSCR,LIKB,NJL,LIKB,NJL)
                    END IF
*
C?                  MXACIJO = MXACIJ
c                    MXACIJ = MAX(MXACIJ,LIKB*NJL,LIKB*NIK)
C?                  IF(MXACIJ.GT.MXACIJO) THEN
C?                    write(6,*) ' New max MXACIJ = ', MXACIJ
C?                    write(6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
C?                    WRITE(6,*) ' ITYP,JTYP,KTYP,LTYP',
C?   &                             ITYP,JTYP,KTYP,LTYP
C?                    WRITE(6,*)'NIJ NJL NIBTC NKBTC',
C?   &                           NIJ,NJL,NIBTC,NKBTC
C?                  END IF
*
                    FACTORC = 0.0D0
                    FACTORAB = FACX
                    CALL MATML7(   SSCR,   CSCR,   XINT,   LIKB,    NIK,
     &                             LIKB,    NJL,    NIK,    NJL,FACTORC,
     &                          FACTORAB,     2)
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' SSCR matrix '
                      CALL WRTMAT(SSCR,LIKB,NIK,LIKB,NIK)
                    END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                    II12 = 2
*
                    IONE = 1
                    IF(IFRST.EQ.1) KFRST = 1
                    ONE = 1.0D0
*
                    IAC = I4_AC(1)
                    KAC = I4_AC(2)
*
C                   KFRST = 1
                    CALL ADAADAST_GAS( IONE,  ISM, ITYP,   NI,  IAC,
     &                                 IONE,  KSM, KTYP,   NK,  KAC,
     &                                ISCTP,ISCSM, IGRP, KBOT, KTOP,
     &                                   I1, XI1S, MAXK,NKBTC, KEND,
     &                                IFRST,
*
     &                                KFRST,
     &                                 II12,
     &                                  K12,
     &                                ONE          )
*
                    IFRST = 0
                    KFRST = 0
*
                    I = 0
                    K = 1
                    DO IK = 1, NIK
                      CALL NXTIJ(     I,     K,    NI,    NK,  IKSM,
     &                            NONEW)
                      IKOFF = (K-1)*NI + I
                      ISBOFF = 1+(IK-1)*NIBTC*NKBTC
                      CALL MATCAS(SSCR(ISBOFF),SB,NIBTC, NROW, IBOT,
     &                             NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                    END DO
C                   write(6,*) ' first element of updated SB', SB(1)
*
                  IF(KEND.EQ.0) GOTO 2800
*. End of loop over partitionings of resolution strings
 2801           CONTINUE
*               ^ End of loop over batches of I strings
              END IF
*             ^ End of if I. ge. K, J.ge. L
 2803         CONTINUE
              END DO
*             ^ End of loop over KSM
            END DO
*           ^ End of loop over ISM
 2950     CONTINUE
        END IF
*       ^ End of a+ a+ a a/a a a+ a+ versus a+ a a+ a switch

 2000 CONTINUE
*
 2001 CONTINUE
*
C?      WRITE(6,*) ' Memcheck at end of RSBB2A '
C?      CALL MEMCHK
C?      WRITE(6,*) ' Memcheck passed '
*
      Call Free_Work(ipSCR)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NSCOL)
        CALL Unused_integer(NSMDX)
      END IF
      END
