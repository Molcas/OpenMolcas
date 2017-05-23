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
      SUBROUTINE GSBBD2A_LUCIA(   RHO2,  RHO2S,  RHO2A,  NACOB,  ISCSM,
     &                           ISCTP,  ICCSM,  ICCTP,   IGRP,   NROW,
     &                            NGAS,   ISEL,   ICEL,     SB,     CB,
     &                           ADSXA, SXSTST, STSTSX, SXDXSX,MXPNGAS,
     &                          NOBPTS, IOBPTS,   MAXI,   MAXK,   SSCR,
*
     &                            CSCR,     I1,   XI1S,     I2,   XI2S,
     &                               X,  NSMOB,  NSMST,  NSMSX, MXPOBS,
     &                          SCLFAC,  IPACK)
*
* SUBROUTINE GSBBD2A_LUCIA --> 37
*
*
* Contributions to two-electron density matrix from column excitations
*
* GAS version, '96 , Jeppe Olsen
*
* =====
* Input
* =====
* RHO2  : two body density matrix to be updated
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
* IPACK  : Should we pack the density?
*
* ======
* Output
* ======
* RHO2, RHO2S, RHO2A : Updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S, I2,XI2S : For holding creations/annihilations
*              type and symmetry
*
* Jeppe Olsen, Fall of 96
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "para_info.fh"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST), SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*)
      LOGICAL IPACK
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*),X(*)
*.Output
      DIMENSION RHO2(*),RHO2S(*),RHO2A(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
C     INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
      REAL*8, PARAMETER :: ZERO = 0.0D0
*

      NTEST = 000
      IF(NTEST.GE.1000) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' GSBBD2A in action '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of active left strings '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' Occupation of active Right strings '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
      IFRST = 1
      JFRST = 1
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
* Type of single excitations that connects the two column strings
      CALL DXTYP_GAS(    NDXTP,      ITP,      JTP,      KTP,      LTP,
     &                    NGAS,     ISEL,     ICEL)
*.Symmetry of Double excitation that connects IBSM and JBSM
*. For general use : STSTSX => STSTDX
      IDXSM = STSTSX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0) GOTO 2001
      IF(NTEST.GE.1000)
     &WRITE(6,*) ' ISCSM,ICCSM ', ISCSM,ICCSM
      DO 2000 IDXTP =  1, NDXTP
        ITYP = ITP(IDXTP)
        JTYP = JTP(IDXTP)
        KTYP = KTP(IDXTP)
        LTYP = LTP(IDXTP)
        IF(NTEST.GE.1000)
     &  write(6,*) ' ITYP JTYP KTYP LTYP ', ITYP,JTYP,KTYP,LTYP
        DO 1950 IKOBSM = 1, NSMOB
          JLOBSM = SXDXSX(IKOBSM,IDXSM)
          IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
          KFRST = 1
*. Loop over of symmetry of i orbitals
          DO 1940 ISM = 1, NSMOB
          KSM = ADSXA(ISM,IKOBSM)
          NI = NOBPTS(ITYP,ISM)
          NK = NOBPTS(KTYP,KSM)
          IF(NI.EQ.0.OR.NK.EQ.0) GOTO 1940
*. Loop over batches of j orbitals
          DO 1930 JSM = 1, NSMOB
          LSM = ADSXA(JSM,JLOBSM)
          NJ = NOBPTS(JTYP,JSM)
          NL = NOBPTS(LTYP,LSM)
          IF(NJ.EQ.0.OR.NL.EQ.0) GOTO 1930
*
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          KOFF = IOBPTS(KTYP,KSM)
          LOFF = IOBPTS(LTYP,LSM)
*
          IF(IOFF.LT.KOFF) GOTO 1930
          IF(JOFF.LT.LOFF) GOTO 1930
*
*
* =========================================================================
*                    Use N-2 projection method
* =========================================================================
*
              IFIRST = 1
*. Loop over batches of I strings
              CALL SETVEC(X,ZERO,NI*NJ*NK*NL)
              DO 1801 IIPART = 1+MYRANK, NPART, NPROCS
                IBOT = 1+(IIPART-1)*NPARTSZ
                ITOP = MIN(IBOT+NPARTSZ-1,NROW)
                NIBTC = ITOP-IBOT+1
                IF (NIBTC.LE.0) EXIT
*.Loop over batches of intermediate strings
                KBOT = 1- MAXK
                KTOP = 0
 1800           CONTINUE
                  KBOT = KBOT + MAXK
                  KTOP = KTOP + MAXK
*
* =========================================================
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
* =========================================================
*
                  IONE = 1
                  JLBOFF = 1
                  IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                    NJL = NJ*(NJ+1)/2
                    JLSM = 1
                  ELSE
                    NJL = NJ * NL
                    JLSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
                  CALL QENTER('ADADS')
                  II12 = 1
                  K12 = 1
                  IONE = 1
                  CALL ADADST_GAS(  IONE,   JSM,  JTYP,    NJ,  IONE,
     &                               LSM,  LTYP,    NL, ICCTP, ICCSM,
     &                              IGRP,  KBOT,  KTOP,    I1,  XI1S,
     &                              MAXK, NKBTC,  KEND, JFRST, KFRST,
     &                              II12,   K12,SCLFAC)
*
                  JFRST = 0
                  KFRST = 0
*
                  CALL QEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
                  CALL QENTER('MATCG')
                  DO  IJL = 1, NJL
                    CALL NXTIJ(      J,      L,     NJ,     NL,   JLSM,
     &                           NONEW)
                    I1JL = (L-1)*NJ+J
*. JAN28
                    IF(JLSM.NE.0) THEN
                      IJLE = J*(J-1)/2+L
                    ELSE
                      IJLE = IJL
                    END IF
*. JAN28
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
C                   JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                    JLOFF = (JLBOFF-1+IJLE-1)*NKBTC*NIBTC+1
                    IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                      CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(    CB,CSCR(JLOFF),NROW,NIBTC,IBOT,
     &                            NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                    END IF
                  END DO
                  CALL QEXIT ('MATCG')
*
*
* =========================================================
*
*. obtain sb(KB,IA,ik) = sum(IB)<KB!a kb a ib !IB>S(IA,IB)
*
* =========================================================
*
                  IONE = 1
                  IKBOFF = 1
                  IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                    NIK = NI*(NI+1)/2
                    IKSM = 1
                  ELSE
                    NIK = NI * NK
                    IKSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
CT                CALL QENTER('ADADS')
                  II12 = 2
                  K12 = 1
                  IONE = 1
                  IF(IFRST.EQ.1) KFRST = 1
                  ONE = 1.0D0
                  CALL ADADST_GAS(  IONE,   ISM,  ITYP,    NI,  IONE,
     &                               KSM,  KTYP,    NK, ISCTP, ISCSM,
     &                              IGRP,  KBOT,  KTOP,    I1,  XI1S,
     &                              MAXK, NKBTC,  KEND, IFRST, KFRST,
     &                              II12,   K12,ONE   )
*
                  IFRST = 0
                  KFRST = 0
*
CT                CALL QEXIT('ADADS')
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  I = 0
                  K = 1
*
CT                CALL QENTER('MATCG')
                  DO  IIK = 1, NIK
                    CALL NXTIJ(      I,      K,     NI,     NK,   IKSM,
     &                           NONEW)
                    I1IK = (K-1)*NI+I
*. JAN28
                    IF(IKSM.NE.0) THEN
                      IIKE = I*(I-1)/2+K
                    ELSE
                      IIKE = IIK
                    END IF
*. JAN28
*.SB(IA,KB,ik) = +/-S(IA,a+ka+iIA)
C                   IKOFF = (IKBOFF-1+IIK-1)*NKBTC*NIBTC+1
                    IKOFF = (IKBOFF-1+IIKE-1)*NKBTC*NIBTC+1
                    IF(IKSM.EQ.1.AND.I.EQ.K) THEN
*. a+j a+j gives trivially zero
                      CALL SETVEC(SSCR(IKOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(    SB,SSCR(IKOFF),NROW,NIBTC,IBOT,
     &                            NKBTC,I1(1,I1IK),XI1S(1,I1IK))
                    END IF
                  END DO
CT                CALL QEXIT ('MATCG')
*
*
* =================================================================
*
* RHO2C(ik,jl)  = RHO2C(ik,jl) - sum(Ia,Kb)SB(Ia,Kb,ik)*CB(Ia,Kb,jl)
*
* =================================================================
*
* The minus ??
*
* Well, the density matrices are constructed as

* <I!a+i a+k aj al!> = -sum(K) <I!a+ia+k!K><J!aj al!K>, and
* the latter matrices are the ones we are constructing
*
              IOFF = IOBPTS(ITYP,ISM)
              JOFF = IOBPTS(JTYP,JSM)
              KOFF = IOBPTS(KTYP,KSM)
              LOFF = IOBPTS(LTYP,LSM)
              NTESTO = NTEST
C?            IF(IOFF.EQ.3.AND.JOFF.EQ.3.AND.KOFF.EQ.4.AND.LOFF.EQ.4)
C?   &            NTEST = 5000
                  LDUMMY = NKBTC*NIBTC
                  IF(NTEST.GE.2000) THEN
                    WRITE(6,*) ' CSCR matrix '
                    CALL WRTMAT(CSCR,LDUMMY,NJL,LDUMMY,NJL)
                    WRITE(6,*) ' SSCR matrix '
                    CALL WRTMAT(SSCR,LDUMMY,NIK,LDUMMY,NIK)
                  END IF

                  IF(IFIRST.EQ.1) THEN
                    FACTOR = 0.0D0
                  ELSE
                    FACTOR = 1.0D0
                  END IF
C                 MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
                  LDUMMY = NKBTC*NIBTC
                  ONEM = -1.0D0
                  CALL MATML7(      X,   SSCR,   CSCR,    NIK,    NJL,
     &                         LDUMMY,    NIK, LDUMMY,    NJL, FACTOR,
     &                           ONEM,      1)
                  IFIRST = 0
                  IF(NTEST.GE.2000) THEN
      WRITE(6,*) ' Updated X matrix IK,JL,IK,JL',NIK,NJL,NIK,NJL
                    CALL WRTMAT(X,NIK,NJL,NIK,NJL)
                  END IF

*
                IF(KEND.EQ.0) GOTO 1800
*. End of loop over partitionings of resolution strings
 1801         CONTINUE
*. Rho2(ik,jl) has been constructed for ik,jl belonging to
*. Scatter out to density matrix
              IOFF = IOBPTS(ITYP,ISM)
              JOFF = IOBPTS(JTYP,JSM)
              KOFF = IOBPTS(KTYP,KSM)
              LOFF = IOBPTS(LTYP,LSM)
!      write(*,*)"I, J, K, L offsets & IPACK :",IOFF,JOFF,KOFF,LOFF,Ipack
              CALL ADTOR2(    RHO2,   RHO2S,   RHO2A,       X,       1,
     &                          NI,    IOFF,      NJ,    JOFF,      NK,
     &                        KOFF,      NL,    LOFF,   NACOB,   IPACK)
C                  ADTOR2(RHO2,RHO2T,ITYPE,
C    &                  NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)

!              write(6,*) ' updated density matrix A'," norb = 4 ",norb
!         write(6,*) ' offset ',"IOFF,JOFF,KOFF,LOFF",IOFF,JOFF,KOFF,LOFF
!              call prsym(rho2s,NORB*(NORB+1)/2)

 1930       CONTINUE
 1940     CONTINUE
 1950   CONTINUE
 2000 CONTINUE
 2001 CONTINUE
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(SXSTST)
        CALL Unused_integer_array(I2)
        CALL Unused_real_array(XI2S)
      END IF
      END
