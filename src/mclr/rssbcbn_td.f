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
      SUBROUTINE RSSBCBN_td(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                  IAEL1,IAEL3,IBEL1,IBEL3,
     &                  JAEL1,JAEL3,JBEL1,JBEL3,
     &                  NAEL,NBEL,
     &                  IJAGRP,IJBGRP,
     &                  SB,CB,IDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  XINT,C2,NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IPRNT,IST,
     &                  CJRES,SIRES,NOPART,TimeDep)
*
* Contributions to sigma block (iasm iatp, ibsm ibtp ) from
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* IAEL1,IAEL3 : Number of elecs in RAS1(RAS3) for alpha strings in sigma
* IBEL1,IBEL3 : Number of elecs in RAS1(RAS3) for  beta strings in sigma
* JAEL1,JAEL3 : Number of elecs in RAS1(RAS3) for alpha strings in C
* JBEL1,JBEL3 : Number of elecs in RAS1(RAS3) for  beta strings in C
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* IDOH2 : = 0 => no two electron operator
* IDOH2 : = 1 =>    two electron operator
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* IST, IDOH2 : See RASSG3 input description
*
* ======
* Output
* ======
* SB : fresh sigma block
*
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
      Logical TimeDep
*. Output
      DIMENSION CB(*),SB(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*),
     &          I3(MAXK,*),XI3S(MAXK,*),I4(MAXK,*),XI4S(MAXK,*)
      DIMENSION C2(*),CJRES(*),SIRES(*)
*
      NTEST =  0
      NTEST = MAX(NTEST,IPRNT)
      icheck=3
*
      ieaw=0
      if (ist.eq.2) ieaw=1
*
*
* =============================
* Sigma beta beta contribution
* =============================
*
* Sigma aa(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*                 * ((ij!kl)-(il!kj)) C(IA,JB)
*                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*
      IF(IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
*         One electron part
*
         IF(IST.EQ.1) THEN
           SIGN = 1.0D0
         ELSE
           SIGN = -1.0D0
         END IF
         IF(NBEL.GE.1) THEN
            CALL RSBB1E(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &         IBEL1,IBEL3,JBEL1,JBEL3,
     &         SB,CB,
     &         ADSXA,SXSTST,STSTSX,NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,XINT,
     &         NSMOB,NSMST,NSMSX,MXPOBS,SIGN)
*
*               Call RECPRT('SSCR after RSBB1E',' ',SSCR,5,1)
         END IF
*
*         Two electron part
*
         IF((iand(icheck,1).eq.1).and.IDOH2.NE.0.AND.NBEL.GE.2)
     &   THEN
*         Write(*,*)'Timedep in rssbcbn_td',TimeDep
         CALL RSBB2A(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
     &                IBEL1,IBEL3,JBEL1,JBEL3,
     &                SB,CB,
     &                ADSXA,DXSTST,STSTDX,SXDXSX,
     &                NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,XINT,
     &                NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,SIGN,
     &                NOPART,TimeDep,ieaw )
*
*               Call RECPRT('SSCR after RSBB2A',' ',SSCR,5,1)
         END IF
      END IF
*
*====================================*
* Mixed alpha-beta double excitations
*====================================*
*
      IF((iand(icheck,2).eq.2).and.IDOH2.NE.0.AND.
     &    NAEL.GE.1.AND.NBEL.GE.1) THEN
*
          ieaw=0
          if (ist.eq.2) ieaw=1
*          Write(*,*)'ieaw in rssbcbn_td ',ieaw
          CALL TRPMAT(CB,NJA,NJB,C2)
          CALL DCOPY_(NJA*NJB,C2,1,CB,1)
          CALL TRPMAT(SB,NIA,NIB,C2)
          CALL DCOPY_(NIA*NIB,C2,1,SB,1)
          IIITRNS = 1
          IF(IIITRNS.EQ.1.AND.NIB.GT.NIA.AND.NJB.GT.NJA) THEN
             JJJTRNS = 1
          ELSE
             JJJTRNS = 0
          END IF
          IF(JJJTRNS.EQ.1.AND.IST.EQ.2) THEN
            IFACTOR = -1
          ELSE
            IFACTOR = 1
          END IF
          IF (JJJTRNS.EQ.0) THEN
          CALL RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                IJAGRP,IJBGRP,
     &                IAEL1,IAEL3,JAEL1,JAEL3,
     &                IBEL1,IBEL3,JBEL1,JBEL3,
     &                SB,CB,
     &                ADSXA,STSTSX,
     &                NTSOB,IBTSOB,ITSOB,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                XINT,
     &                NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,0,1,
     &                CJRES,SIRES,C2,NTEST,IFACTOR,ieaw,
     &                TimeDep )
*
*          Call RECPRT('SSCR after RSBB2BN',' ',SSCR,5,1)
*
          ELSE IF ( JJJTRNS.EQ.1) THEN
            CALL TRPMAT(SB,NIB,NIA,C2)
            CALL DCOPY_(NIA*NIB,C2,1,SB,1)
            CALL TRPMAT(CB,NJB,NJA,C2)
            CALL DCOPY_(NJA*NJB,C2,1,CB,1)
*
            CALL RSBB2BN(IBSM,IBTP,IASM,IATP,NIB,NIA,
     &                JbSM,JbTP,JaSM,JaTP,NJb,NJa,
     &                IJbGRP,IJaGRP,
     &                IbEL1,IbEL3,JbEL1,JbEL3,
     &                IaEL1,IaEL3,JaEL1,JaEL3,
     &                SB,CB,
     &                ADSXA,STSTSX,
     &                NTSOB,IBTSOB,ITSOB,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                XINT,
     &                NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,0,1,
     &                CJRES,SIRES,C2,NTEST,IFACTOR,ieaw,
     &                TimeDep)
*
*                Call RECPRT('SSCR after RSBB2BN',' ',SSCR,5,1)
*
            CALL TRPMAT(SB,NIA,NIB,C2)
            CALL DCOPY_(NIA*NIB,C2,1,SB,1)
            CALL TRPMAT(CB,NJA,NJB,C2)
            CALL DCOPY_(NJA*NJB,C2,1,CB,1)
          END IF
*.        Restore order !
          CALL TRPMAT(CB,NJB,NJA,C2)
          CALL DCOPY_(NJA*NJB,C2,1,CB,1)
          CALL TRPMAT(SB,NIB,NIA,C2)
          CALL DCOPY_(NIA*NIB,C2,1,SB,1)
      END IF
*
* =============================
* Sigma alpha contribution
* =============================
*
*. Transpose for alpha excitations
*
      IF(NAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
           CALL TRPMAT(CB,NJA,NJB,C2)
           CALL DCOPY_(NJA*NJB,C2,1,CB,1)
           CALL TRPMAT(SB,NIA,NIB,C2)
           CALL DCOPY_(NIA*NIB,C2,1,SB,1)
*
* alpha single excitation
*
           SIGN = 1.0D0
           CALL RSBB1E(IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &                IAEL1,IAEL3,JAEL1,JAEL3,
     &                SB,CB,
     &                ADSXA,SXSTST,STSTSX,
     &                NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,XINT,
     &                NSMOB,NSMST,NSMSX,MXPOBS,SIGN)
*
*                Call RECPRT('SSCR after RSBB1E',' ',SSCR,5,1)
*
* alpha double excitation
*
           IF((iand(icheck,1).eq.1).and.NAEL.GE.2.AND.IDOH2.NE.0)
     &     Then
*            Write(*,*)'Timedep in rssbcbn_td',TimeDep
            CALL RSBB2A(IASM,IATP,JASM,JATP,IJAGRP,NIB,
     &           IAEL1,IAEL3,JAEL1,JAEL3,
     &           SB,CB,
     &           ADSXA,DXSTST,STSTDX,SXDXSX,
     &           NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &           SSCR,CSCR,I1,XI1S,XINT,
     &           NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,SIGN,
     &           NOPART,TimeDep,ieaw)
*
*                Call RECPRT('SSCR after RSBB2A',' ',SSCR,5,1)
*
           END IF
* Restore order !
           CALL TRPMAT(SB,NIB,NIA,C2)
           CALL DCOPY_(NIA*NIB,C2,1,SB,1)
           CALL TRPMAT(CB,NJB,NJA,C2)
           CALL DCOPY_(NJA*NJB,C2,1,CB,1)
      END IF
*
      RETURN
      END
