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
      SUBROUTINE RSBB1E(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  ISEL1,ISEL3,ICEL1,ICEL3,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,
     &                  NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,SIGN)
*
* One electron excitations on column strings
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* ISEL1(3) : Number of electrons in RAS1(3) for S block
* ICEL1(3) : Number of electrons in RAS1(3) for C block
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
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST,
     &        STSTSX(NSMST,NSMST)
      INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),H(*)
*.Local arrays
      DIMENSION ITP(3),JTP(3)
      ZERO = 0.0D0
      ONE = 1.0D0
      ONEM = -1.0D0
*
* Type of single excitations that connects the two column strings
*
      CALL SXTYP(NSXTP,ITP,JTP,ISEL1,ISEL3,ICEL1,ICEL3)
*
*.Symmetry of single excitation that connects IBSM and JBSM
*
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(IJSM.EQ.0) GOTO 1001
*
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        DO 800 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 800
          NIORB = NTSOB(ITYP,ISM)
          NJORB = NTSOB(JTYP,JSM)
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*. Loop over partitionings of the row strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
          DO 701 IPART = 1, NIPART
            IBOT = (IPART-1)*MAXI+1
            ITOP = MIN(IBOT+MAXI-1,NROW)
            NIBTC = ITOP - IBOT + 1
*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = KTOP + MAXK
              JBORB = IBTSOB(JTYP,JSM)
              NJORB = NTSOB(JTYP,JSM)
* Obtain CSCR(JORB,I,K) = SUM(J)<K!A JORB!J>C(I,J)
* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
              DO 500 JJORB = 1,NJORB
*               IF(JJORB.EQ.1) THEN
*                 JOFF = 1
*               ELSE
*                 JOFF = IOFF + NK
*. Note : NK is not not known until first call to ADST
*               END IF
                JORB = ITSOB(JBORB-1+JJORB)
*.set up I1(K) =  XI1S(K) a JORB !J STRING >
                CALL ADST(JORB,1,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &               I1,XI1S,MAXK,NKBTC,KEND)
*.Gather  C Block
*. First index : JORB, second index : JaKb
                ICGOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I1,XI1S)
  500         CONTINUE
*.Obtain one electron integrals (JSM,JTP,ISM,ITP)
              CALL NGETH1(SSCR,ISM,ITYP,JSM,JTYP)
              CALL TRPMAT(SSCR,NIORB,NJORB,H)
*
*.Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
                NIK = NIBTC*NKBTC
*
              CALL  DGEMM_('N','N',NIK,NIORB,NJORB,ONE ,
     &                    CSCR,max(NIK,1),H,max(1,NJORB),
     &                    ZERO,SSCR,max(1,NIK))
*.S(I,a+ K) =  S(I, a+ K) + sgn*Sscr(I,K,i)
              IBORB = IBTSOB(ITYP,ISM)
              DO 400 IIORB = 1,NIORB
                IORB = ITSOB(IBORB-1+IIORB)
*.set up I1(IORB,K) = a IORB !I STRING >
                CALL ADST(IORB,1,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &          I1,XI1S,1,NKBTC,KEND)
*. Well, someplace the minus must come in
                IF( SIGN .EQ. -1.0D0) THEN
                  CALL DSCAL_(NKBTC,ONEM,XI1S,1)
                END IF
                ISBOFF = 1+(IIORB-1)*NKBTC*NIBTC
                CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &               NKBTC,I1,XI1S)
  400         CONTINUE
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
  701     CONTINUE
*. End of loop over I partitioninigs
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(SXSTST)
        CALL Unused_integer(NSMSX)
      END IF
      END
