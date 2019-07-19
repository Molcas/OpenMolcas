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
* Copyright (C) 1991,1995, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE GSBBD1(RHO1,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,RHO1S)
*


* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen
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
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(3,*), IOBPTS(3,*), ITSOB(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(*),XI1S(*)
      DIMENSION I2(*),XI2S(*)
*
*.Local arrays
      DIMENSION ITP(3*3),JTP(3*3)
*
      IFRST = 1
      JFRST = 1
*
* Type of single excitations that connects the two column strings
      CALL SXTYP_GAS(NSXTP,ITP,JTP,3,ISEL,ICEL)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        DO 800 ISM = 1, NSMOB
*. new i and j so new intermediate strings
          KFRST = 1
*
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 800
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
          IBIORB = IOBPTS(ITYP,ISM)
          IBJORB = IOBPTS(JTYP,JSM)
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
COLD. Loop over partitionings of the row strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = KTOP + MAXK
*EAWBEGIN970207
*             -1 -> MAXK
*             KTOP -> -1
*             KTOP=-1
*. Single excitation information independent of I strings
*
*.set up I1(K) =  XI1S(K) a JORB !J STRING >
              CALL ADST(IBJORB,NJORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &             I1,XI1S,MAXK,NKBTC,KEND)
*.set up I2(K) =  XI1S(K) a JORB !J STRING >
              CALL ADST(IBIORB,NIORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &             I2,XI2S,MAXK,NKBTC,KEND)
*EAWEND
*. Appropriate place to start partitioning over I strings
*. Loop over partitionings of the row strings
          DO 701 IPART = 1, NIPART
            IBOT = (IPART-1)*MAXI+1
            ITOP = MIN(IBOT+MAXI-1,NROW)
            NIBTC = ITOP - IBOT + 1

* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
*.Gather  C Block
              DO JJORB = 1,NJORB
                ICGOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I1(1+(JJORB-1)*MAXK),
     &                         XI1S(1+(JJORB-1)*MAXK)     )
              END DO
*
*
* Obtain SSCR(I,K,IORB) = SUM(I)<K!A IORB!J>S(I,J)
              DO IIORB = 1,NIORB
*.Gather S Block
                ISGOFF = 1 + (IIORB-1)*NKBTC*NIBTC
                CALL MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I2(1+(IIORB-1)*MAXK),
     &                     XI2S(1+(IIORB-1)*MAXK) )
              END DO
              NKI = NKBTC*NIBTC
              If (NKI*NIORB*NJORB.ne.0)  Then
              CALL DGEMM_('T','N',NIORB,NJORB,NKI,1.0D0,SSCR,NKI,CSCR,
     &                    NKI,0.0D0,RHO1S,NIORB)
              Else
               call dcopy_(NIORB*NJORB,[0.0d0],0,RHO1S,1)
              End IF
*. Scatter out to complete matrix
              DO 610 JJORB = 1, NJORB
                JORB = IBJORB-1+JJORB
                DO 605 IIORB = 1, NIORB
                  IORB = IBIORB-1+IIORB
                  RHO1((JORB-1)*NACOB+IORB) =
     &            RHO1((JORB-1)*NACOB+IORB) +
     &            RHO1S((JJORB-1)*NIORB+IIORB)
 605            CONTINUE
 610         CONTINUE
*
*
  701     CONTINUE
*. /\ end of this I partitioning
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
*. End of loop over I partitioninigs
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(SXSTST)
        CALL Unused_integer(MXPNGAS)
        CALL Unused_integer_array(ITSOB)
        CALL Unused_real(H)
      END IF
      END
