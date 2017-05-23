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
* Copyright (C) 1991-1994, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                  JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,
     &                  IAEL1,IAEL3,JAEL1,JAEL3,
     &                  IBEL1,IBEL3,JBEL1,JBEL3,
     &                  SB,CB,
     &                  ADSXA,STSTSX,
     &                  NTSOB,IBTSOB,ITSOB,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  ICJKAIB,CJRES,SIRES,S2,NTEST,ISIGN,ieaw,
     &                  TimeDep)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
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
*
* ICJKAIB =1 =>  construct C(Ka,Jb,j) and S(Ka,Ib,i) as intermediate
*                 matrices in order to reduce overhead
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
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + ICJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      Logical TimeDep
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*),I1(MAXK,*),XI1S(MAXK,*)
      DIMENSION I2(MAXK,*),XI2S(MAXK,*)
      DIMENSION I3(MAXK,*),XI3S(MAXK,*)
      DIMENSION I4(MAXK,*),XI4S(MAXK,*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION S2(*)
*.Local arrays
      DIMENSION ITP(3),JTP(3),KTP(3),LTP(3)
*
      ZERO = 0.0D0
      ONEM = -1.0D0
*     IUSEAB = 0
*     ICJKAIB = 1
      IROUTE = 1
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
*.Types of SX that connects the two strings
      CALL SXTYP(NKLTYP,KTP,LTP,IBEL1,IBEL3,JBEL1,JBEL3)
      CALL SXTYP(NIJTYP,ITP,JTP,IAEL1,IAEL3,JAEL1,JAEL3)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
****. TESTTING
        N1IND = 0
        N2IND = 0
        N3IND = 0
        IF(ITYP.EQ.1) N1IND = N1IND + 1
        IF(ITYP.EQ.2) N2IND = N2IND + 1
        IF(ITYP.EQ.3) N3IND = N3IND + 1
*
        IF(JTYP.EQ.1) N1IND = N1IND + 1
        IF(JTYP.EQ.2) N2IND = N2IND + 1
        IF(JTYP.EQ.3) N3IND = N3IND + 1
*
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          IOFF = IBTSOB(ITYP,ISM)
          JOFF = IBTSOB(JTYP,JSM)
          NI = NTSOB(ITYP,ISM)
          NJ = NTSOB(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Loop over batches of KA strings
          KABOT = 1 - MAXK
          KATOP = 0
 1801     CONTINUE
            KABOT = KABOT + MAXK
            KATOP = KATOP + MAXK
*. Find Ka strings that connect with Ja strings for given group of Jorbs
           CALL ADST(JOFF,NJ,JATP,JASM,IAGRP,KABOT,KATOP,
     &              I1(1,1),XI1S(1,1),MAXK,NKABTC,KAEND)
           CALL ADST(IOFF,NI,IATP,IASM,IAGRP,KABOT,KATOP,
     &               I3(1,1),XI3S(1,1),MAXK,NKABTC,KAEND)
           IF(NKABTC.EQ.0) GOTO 1940
*. Generate - if required C(Ka,Jb,j)
           IF(ICJKAIB.NE.0) THEN
             IIOFF = 1
             LCJ = NJB*NKABTC
             DO JJ = 1, NJ
               IF(JJ.EQ.1) THEN
                 IIOFF = 1
               ELSE
                 IIOFF = IIOFF + LCJ
               END IF
                CALL GATRMT(CB,NJA,NJB,CJRES(IIOFF),NKABTC,NJB,
     &                      I1(1,JJ),XI1S(1,JJ) )
             End Do
*
*. We have now C gathered in the form C(Ka,Jb,j).
*. If Ka is small, say 1, it can be advantageous to switch
* around to C(j,Ka,Jb). This is mediated by the switch IROUTE
* IROUTE = 1 : Normal (i.e. old) route,
* IROUTE = 2 : New route with j first
* IROUTE = 3 : C(Ka,j,Jb)
*
             IF(NJ.GE.NKABTC.AND.NI.GE.NKABTC ) THEN
               IROUTE = 2
             ELSE
                IROUTE = 3
             END IF
* DOES THIS WORK?
C9805EAW     IROUTE = 1
             If (TimeDep) IROUTE=1
*
             IF(IROUTE.EQ.2) THEN
*. C(Ka,Jb,j) => C(j,Ka,Jb)
               CALL TRPMAT(CJRES,LCJ,NJ,SIRES)
               CALL DCOPY_(NJ*LCJ,SIRES,1,CJRES,1)
             END IF
             IF(IROUTE.EQ.3) THEN
*. C(Ka,Jb,j) => C(Ka,j,JB)
               DO JB = 1,NJB
                 DO J = 1, NJ
                   IOFFIN = (J-1)*NJB*NKABTC+(JB-1)*NKABTC + 1
                   IOFFOUT =(JB-1)*NKABTC*NJ + (J-1)*NKABTC + 1
                   DO KA = 1, NKABTC
                     SIRES(IOFFOUT-1+KA) = CJRES(IOFFIN-1+KA)
                   END DO
                 END DO
               END DO
               CALL DCOPY_(NJ*LCJ,SIRES,1,CJRES,1)
             END IF
*
             CALL DCOPY_(NIB*NKABTC*NI,ZERO,0,SIRES,1)

           END IF
*
           DO 2000 KLTYP = 1, NKLTYP
             KTYP = KTP(KLTYP)
             LTYP = LTP(KLTYP)
*. Testing
             N1IND = 0
             N2IND = 0
             N3IND = 0
             IF(ITYP.EQ.1) N1IND = N1IND + 1
             IF(ITYP.EQ.2) N2IND = N2IND + 1
             IF(ITYP.EQ.3) N3IND = N3IND + 1

             IF(JTYP.EQ.1) N1IND = N1IND + 1
             IF(JTYP.EQ.2) N2IND = N2IND + 1
             IF(JTYP.EQ.3) N3IND = N3IND + 1
*
             IF(KTYP.EQ.1) N1IND = N1IND + 1
             IF(KTYP.EQ.2) N2IND = N2IND + 1
             IF(KTYP.EQ.3) N3IND = N3IND + 1
*
             IF(LTYP.EQ.1) N1IND = N1IND + 1
             IF(LTYP.EQ.2) N2IND = N2IND + 1
             IF(LTYP.EQ.3) N3IND = N3IND + 1
             DO 1930 KSM = 1, NSMOB
              IFIRST = 1
              LSM = ADSXA(KSM,KLSM)
              IF(LSM.EQ.0) GOTO 1930
              KOFF = IBTSOB(KTYP,KSM)
              LOFF = IBTSOB(LTYP,LSM)
              NK = NTSOB(KTYP,KSM)
              NL = NTSOB(LTYP,LSM)
*. If IUSEAB  is used, only terms with i.ge.k will be generated so
              IKORD = 0
              IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
              IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.GT.KTYP)
     &          GOTO 1930
              IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &           IKORD = 1
*
              IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*.Loop over batches of KB strings
              KBBOT = 1- MAXK
              KBTOP = 0
 1800         CONTINUE
               KBBOT = KBBOT + MAXK
               KBTOP = KBTOP + MAXK
*. obtain cb(KA,KB,jl) =  sum(JA,JB)<KA!a la!JA><KB!a jb !JB>C(JA,JB)
*
               CALL ADST(LOFF,NL,JBTP,JBSM,IBGRP,KBBOT,KBTOP,
     &                      I2(1,1),XI2S(1,1),MAXK,NKBBTC,KBEND)
               CALL ADST(KOFF,NK,IBTP,IBSM,IBGRP,KBBOT,KBTOP,
     &                   I4(1,1),XI4S(1,1),MAXK,NKBBTC,KBEND)
               IF(NKBBTC.EQ.0) GOTO 1930
*
*. Modern low copy version
               IF(IFIRST.EQ.1) THEN
                 IXCHNG = 0
                 IF(IROUTE.EQ.1) THEN
*. Integrals stored as (j l i k )
*                   Write(*,*)'Timedep in rsbb2bn;',TimeDep
                   If (TimeDep) Then
                      CALL GETINT_td(XINT,JTYP,JSM,ITYP,ISM,LTYP,LSM,
     &                       KTYP,KSM,0,0,iroute,ieaw)
                   Else
                      CALL GETINT_MCLR(XINT,JTYP,JSM,ITYP,ISM,LTYP,LSM,
     &                       KTYP,KSM,IXCHNG,0,0,0,ieaw)
                   End If
                 ELSE IF(IROUTE.EQ.2 ) THEN
*. Integrals stored as (i j k l )
                   If (TimeDep) Then
                      CALL GETINT_td(XINT,ITYP,ISM,JTYP,JSM,KTYP,KSM,
     &                         LTYP,LSM,0,0,iroute,ieaw)
                   Else
                       CALL GETINT_MCLR(XINT,ITYP,ISM,JTYP,JSM,KTYP,KSM,
     &                         LTYP,LSM,IXCHNG,0,0,1,ieaw)
                   End If
                 ELSE IF(IROUTE.EQ.3) THEN
*. Integrals stored as (j i k l )
                   If (TimeDep) Then
                        CALL GETINT_td(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,
     &                         LTYP,LSM,0,0,iroute,ieaw)
                   Else
                       CALL GETINT_MCLR(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,
     &                         LTYP,LSM,IXCHNG,0,0,1,ieaw)
                   End If
                 END IF
                 IF(ISIGN.EQ.-1) THEN
                   CALL DSCAL_(NI*NJ*NK*NL,ONEM,XINT,1)
                 END IF
                 IFIRST = 0
               END IF
               CALL SKICKJ(SIRES,CJRES,NKABTC,NIB,NJB,
     &                     NKBBTC,XINT,NI,NJ,NK,NL,MAXK,
     &                     I4,XI4S,I2,XI2S,IKORD,
     &                     IDUM,iXDUM,XDUM,IROUTE,NTEST )
*
               IF(KBEND.EQ.0) GOTO 1800
*. End of loop over partitioning of beta strings
 1930        CONTINUE
 2000      CONTINUE
*. Scatter out from s(Ka,Ib,i)
*. Restore order !!
           IF(IROUTE.EQ.2) THEN
             CALL TRPMAT(SIRES,NI,NIB*NKABTC,CJRES)
             CALL DCOPY_(NI*NIB*NKABTC,CJRES,1,SIRES,1)
           END IF
           IF(IROUTE.EQ.3) THEN
              DO JB = 1,NIB
                 DO J = 1, NI
                  IOFFIN = (J-1)*NIB*NKABTC+(JB-1)*NKABTC + 1
                  IOFFOUT =(JB-1)*NKABTC*NI + (J-1)*NKABTC + 1
                  DO KA = 1, NKABTC
                    CJRES(IOFFIN-1+KA) = SIRES(IOFFOUT-1+KA)
                  END DO
                 END DO
              END DO
              CALL DCOPY_(NI*NIB*NKABTC,CJRES,1,SIRES,1)
           END IF
           IF(ICJKAIB.EQ.1) THEN
              DO II = 1, NI
                CALL SCARMT(SIRES((II-1)*NKABTC*NIB+1),NKABTC,NIB,
     &                      SB,NIA,NIB,I3(1,II),XI3S(1,II))

              End Do
           END IF
           IF(KAEND.EQ.0) GOTO 1801
*. End of loop over partitioning of alpha strings
 1940     CONTINUE
 2001 CONTINUE
*
 9999 CONTINUE
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(ITSOB)
        CALL Unused_real_array(SSCR)
        CALL Unused_real_array(CSCR)
        CALL Unused_integer(NSMSX)
        CALL Unused_integer(NSMDX)
        CALL Unused_real_array(S2)
      END IF
      END
