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

      SUBROUTINE GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,SB,CB,
     &                  ADSXA,SXSTST,STSTSX,SXDXSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,X,
     &                  NSMOB,NSMST,NSMSX,MXPOBS)
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
*
* ======
* Output
* ======
* RHO2 : Updated density block
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
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST), SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO2(*),X(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*.Local arrays
      DIMENSION ITP(3*3),JTP(3*3),KTP(3*3),LTP(3*3)
*
*
* Type of single excitations that connects the two column strings
      NTEST=0
      CALL DXTYP_GAS(NDXTP,ITP,JTP,KTP,LTP,3,ISEL,ICEL)
*.Symmetry of Double excitation that connects IBSM and JBSM
*. For general use : STSTSX => STSTDX
      IDXSM = STSTSX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0) GOTO 2001
      DO 2000 IDXTP =  1, NDXTP
        ITYP = ITP(IDXTP)
        JTYP = JTP(IDXTP)
        KTYP = KTP(IDXTP)
        LTYP = LTP(IDXTP)
        DO 1950 IKOBSM = 1, NSMOB
          JLOBSM = SXDXSX(IKOBSM,IDXSM)
          IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
*         KFRST = 1
*. Loop over of symmetry of i orbitals
          DO 1940 ISM = 1, NSMOB
          KSM = ADSXA(ISM,IKOBSM)
          NI = NOBPTS(ITYP,ISM)
          NK = NOBPTS(KTYP,KSM)
          IOFF = IOBPTS(ITYP,ISM)
          KOFF = IOBPTS(KTYP,KSM)
          IF(NI.EQ.0.OR.NK.EQ.0) GOTO 1940
*. Loop over batches of j orbitals
          DO 1930 JSM = 1, NSMOB
          IFIRST = 1
          LSM = ADSXA(JSM,JLOBSM)
          NJ = NOBPTS(JTYP,JSM)
          NL = NOBPTS(LTYP,LSM)
          JOFF = IOBPTS(JTYP,JSM)
          LOFF = IOBPTS(LTYP,LSM)
          IF(IOFF.LT.KOFF) GOTO 1930
          IF(JOFF.LT.LOFF) GOTO 1930
          IF(NJ.EQ.0.OR.NL.EQ.0) GOTO 1930
*
* =========================================================================
*                    Use N-2 projection method
* =========================================================================
*
*. Loop over batches of I strings
              NPART = NROW/MAXI
              IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
              IF(NTEST.GE.2000)
     &        write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
              DO 1801 IPART = 1, NPART
                IBOT = 1+(IPART-1)*MAXI
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP-IBOT+1
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
                  nkStref=maxk  ! ????????
                  JLBOFF = 1
                  IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                    NJL = NJ*(NJ+1)/2
                    JLSM = 1
                  ELSE
                    NJL = NJ * NL
                    JLSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
              lOFF = IOBPTS(lTYP,lSM)
              jOFF = IOBPTS(jTYP,jSM)
                  CALL ADADST(JTYP,JSM,JOFF,NJ,
     &                        LTYP,LSM,LOFF,NL,jlsm,
     &                        ICCTP,ICCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,
     &                        NKBTC,nkstref,KEND)

*
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
                  DO  IJL = 1, NJL
                    CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
*                   I1JL = (J-1)*NJ+L
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
*. JAN28
                    IF(JLSM.NE.0) THEN
                      IJLE = J*(J-1)/2+L
                    ELSE
                      IJLE = IJL
                    END IF
*
                    JLOFF = (JLBOFF-1+IJLE-1)*NKBTC*NIBTC+1
                    IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                    ELSE
*EAW BEGIN 970407
*                     CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
*    &                            I1(1,I1JL),XI1S(1,I1JL))
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,IJL),XI1S(1,IJL))
*EAW END
                    END IF
                  END DO
*
*
* =========================================================
*
*. obtain sb(KB,IA,ik) = sum(IB)<KB!a kb a ib !IB>S(IA,IB)
*
* =========================================================
*
                  IKBOFF = 1
                  IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                    NIK = NI*(NI+1)/2
                    IKSM = 1
                  ELSE
                    NIK = NI * NK
                    IKSM = 0
                  END IF
*. Obtain all double excitations from this group of K strings
              KOFF = IOBPTS(KTYP,KSM)
              iOFF = IOBPTS(iTYP,iSM)
*                 IF(IFRST.EQ.1) KFRST = 1
                   CALL ADADST(ITYP,ISM,IOFF,NI,
     &                        KTYP,KSM,KOFF,NK,IKSM,
     &                        ISCTP,ISCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,
     &                        NKBTC,NKSTREF,KEND)


*
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  I = 0
                  K = 1
*
                  DO  IIK = 1, NIK
                    CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
*                   I1IK = (K-1)*NI+I
*. JAN28
                    IF(IKSM.NE.0) THEN
                      IIKE = I*(I-1)/2+K
                    ELSE
                      IIKE = IIK
                    END IF
*. JAN28
*.SB(IA,KB,ik) = +/-S(IA,a+ka+iIA)
                    IKOFF = (IKBOFF-1+IIKE-1)*NKBTC*NIBTC+1
                    IF(IKSM.EQ.1.AND.I.EQ.K) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(SSCR(IKOFF),ZERO,NKBTC*NIBTC)
                    ELSE
*EAW-BEGIN 970407
*                     CALL MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,
*    &                            I1(1,I1IK),XI1S(1,I1IK))
                      CALL MATCG(SB,SSCR(IKOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,IIK),XI1S(1,IIK))
*EAW-END
                    END IF
                  END DO
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
                  LDUMMY = NKBTC*NIBTC

                 IF(IFIRST.EQ.1) THEN
                    FACTOR = 0.0D0
                  ELSE
                    FACTOR = 1.0D0
                  END IF
                  LDUMMY = NKBTC*NIBTC
                  ONEM = -1.0D0
                  CALL  DGEMM_('T','N',NIK,NJL,
     &                        LDUMMY,ONEM,SSCR,max(1,LDUMMY),
     &                        CSCR,max(1,LDUMMY),
     &                        FACTOR,X,max(1,NIK))
                  IFIRST = 0

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
              CALL ADTOR2_MCLR(RHO2,X,1,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,
     &                    NACOB)
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
