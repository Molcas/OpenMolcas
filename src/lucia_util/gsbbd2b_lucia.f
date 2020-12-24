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
      SUBROUTINE GSBBD2B_LUCIA(   RHO2,  RHO2S,  RHO2A,   IASM,   IATP,
     &                            IBSM,   IBTP,    NIA,    NIB,   JASM,
     &                            JATP,   JBSM,   JBTP,    NJA,    NJB,
     &                           IAGRP,  IBGRP,   NGAS,   IAOC,   IBOC,
     &                            JAOC,   JBOC,     SB,     CB,  ADSXA,
*
     &                          STSTSX,MXPNGAS, NOBPTS, IOBPTS,   MAXK,
     &                              I1,   XI1S,     I2,   XI2S,     I3,
     &                            XI3S,     I4,   XI4S,      X,  NSMOB,
     &                           NSMST,  NSMSX,  NSMDX, MXPOBS, IUSEAB,
     &                           CJRES,  SIRES,   NORB, NTESTG, SCLFAC,
*
     &                         S2_TERM1, IPACK)
*
* SUBROUTINE GSBBD2B_LUCIA --> 52
*
*
* alpha-beta contribution to two-particle density matrix
* from given c-block and s-block.
*
* S2_TERM1 = - <L!a+i alpha a+jbeta a i beta a j alpha !R>
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
* IPACK  : Should we pack the density?
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
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* X : Space for block of two-electron integrals
*
* Jeppe Olsen, Fall of 1996
*
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "para_info.fh"
#include "WrkSpc.fh"
#include "loff.fh"
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
*.Input
      DIMENSION CB(*),SB(*)
      LOGICAL IPACK
*. Output
      DIMENSION RHO2(*),RHO2S(*),RHO2A(*)
*.Scratch
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION CJRES(*),SIRES(*)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)

      DIMENSION IAOC(*),JAOC(*),IBOC(*),JBOC(*)
*

      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' GSBBD2B speaking '
        WRITE(6,*) ' ================== '
      END IF
      !WRITE(6,*) ' NJAS NJB = ',NJA,NJB
      !WRITE(6,*) ' IAGRP IBGRP = ', IAGRP,IBGRP
      !WRITE(6,*) ' MXPNGAS = ', MXPNGAS
      !WRITE(6,*) ' NSMOB = ', NSMOB

*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP_GAS(   NKLTYP,      KTP,      LTP,     NGAS,     IBOC,
     &                    JBOC)
      CALL SXTYP_GAS(   NIJTYP,      ITP,      JTP,     NGAS,     IAOC,
     &                    JAOC)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
* Repeated allocation/deallocation inside ADSTN_GAS has been
* outerlooped to here. KLOFFI added to call parameters of
* ADSTN_GAS. PAM March 2006.
      CALL GETMEM('KLOFFI','ALLO','REAL',KLOFFI,LOFFI)
      Call FZero(Work(KLOFFI),LOFFI)

      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          ntest=00 !yjma
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all Ka strings
*. a+j!ka> = +/-/0 * !Ja>
          CALL ADSTN_GAS(  KLOFFI,     JSM,    JTYP,    JATP,    JASM,
     &                      IAGRP,      I1,    XI1S,  NKASTR,    IEND,
     &                      IFRST,   KFRST,    KACT,  SCLFAC)
          IF (NKASTR.EQ.0) GOTO 1940
*. a+i!ka> = +/-/0 * !Ia>
          ONE    = 1.0D0
          CALL ADSTN_GAS(  KLOFFI,     ISM,    ITYP,    IATP,    IASM,
     &                      IAGRP,      I3,    XI3S,  NKASTR,    IEND,
     &                      IFRST,   KFRST,    KACT,  ONE   )
          IF (NKASTR.EQ.0) GOTO 1940
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
C             COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
              CALL COMPRS2LST(     I1,   XI1S,     NJ,     I3,   XI3S,
     &                             NI, NKASTR, NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF

*. Loop over batches of KA strings
          NKABTC=0
          DO
            NKABTC=NKABTC+NPROCS
            NKABTCSZ=MAX(NKAEFF-1,0)/NKABTC+1
            IF (NKABTCSZ.LE.MAXK) EXIT
          END DO
*
          DO 1801 IKABTC = 1+MYRANK, NKABTC, NPROCS
            KABOT = (IKABTC-1)*NKABTCSZ + 1
            KATOP = MIN(KABOT+NKABTCSZ-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
            IF (LKABTC.LE.0) EXIT
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, NJ
              CALL GET_CKAJJB(     CB,     NJ,    NJA,  CJRES, LKABTC,
     &                            NJB,
     &                             JJ,
     &                        I1(KABOT+(JJ-1)*NKASTR),
     &                        XI1S(KABOT+(JJ-1)*NKASTR))
*
            END DO
*. Obtain S(ka,i,Ib) for Ka in batch
            DO II = 1, NI
              CALL GET_CKAJJB(     SB,     NI,    NIA,  SIRES, LKABTC,
     &                            NIB,
     &                             II,
     &                        I3(KABOT+(II-1)*NKASTR),
     &                        XI3S(KABOT+(II-1)*NKASTR))
*
            END DO
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*
              DO 1930 KSM = 1, NSMOB
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
                KOFF = IOBPTS(KTYP,KSM)
                LOFF = IOBPTS(LTYP,LSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP) IKORD=1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
                ONE = 1.0D0
                CALL ADSTN_GAS( KLOFFI,    LSM,   LTYP,   JBTP,   JBSM,
     &                           IBGRP,     I2,   XI2S, NKBSTR,   IEND,
     &                           IFRST,  KFRST,   KACT, ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADSTN_GAS( KLOFFI,    KSM,   KTYP,   IBTP,   IBSM,
     &                           IBGRP,     I4,   XI4S, NKBSTR,   IEND,
     &                           IFRST,  KFRST,   KACT,    ONE)
                IF(NKBSTR.EQ.0) GOTO 1930
*
*. Update two-electron density matrix
*  Rho2b(ij,kl) =  Sum(ka)S(Ka,i,Ib)<Ib!Eb(kl)!Jb>C(Ka,j,Jb)
*
                ZERO = 0.0D0
                CALL SETVEC(X,ZERO,NI*NJ*NK*NL)
*
                CALL ABTOR2(  SIRES,  CJRES, LKABTC,    NIB,    NJB,
     &                       NKBSTR,      X,     NI,     NJ,     NK,
     &                           NL, NKBSTR,     I4,   XI4S,     I2,
     &                         XI2S,  IKORD)
*. contributions to Rho2(ij,kl) has been obtained, scatter out
              !call wrtmat(x,ni*nj,nk*nl,ni*nj,nk*nl)
*. Contribution to S2
                IF(KTYP.EQ.JTYP.AND.KSM.EQ.JSM.AND.
     &            ITYP.EQ.LTYP.AND.ISM.EQ.LSM) THEN
                  DO I = 1, NI
                    DO J = 1, NJ
                      IJ = (J-1)*NI+I
                      JI = (I-1)*NJ+J
                      NIJ = NI*NJ
                      S2_TERM1 = S2_TERM1-X((JI-1)*NIJ+IJ)
                    END DO
                  END DO
                END IF

                CALL ADTOR2(   RHO2,  RHO2S,  RHO2A,      X,      2,
     &                           NI,   IOFF,     NJ,   JOFF,     NK,
     &                         KOFF,     NL,   LOFF,   NORB,  IPACK)

!              write(6,*) ' updated density matrix B ',"norb = ",norb
!         write(6,*) ' offset ',"IOFF,JOFF,KOFF,LOFF",IOFF,JOFF,KOFF,LOFF
!              call prsym(rho2s,NORB*(NORB+1)/2)

 1930         CONTINUE
 2000       CONTINUE
 1801     CONTINUE
*. End of loop over partitioning of alpha strings
 1940   CONTINUE
 2001 CONTINUE
* This 'flush' outerlooped here. Was previously inside ADSTN_GAS.
      CALL GETMEM('KLOFFI','FREE','REAL',KLOFFI,LOFFI)
*
 9999 CONTINUE
*
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NSMSX)
        CALL Unused_integer(NSMDX)
      END IF
      END
