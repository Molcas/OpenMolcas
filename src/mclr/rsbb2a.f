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
      SUBROUTINE RSBB2A(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
     &                  ISEL1,ISEL3,ICEL1,ICEL3,
     &                  SB,CB,
     &                  ADSXA,DXSTST,STSTDX,SXDXSX,
     &                  NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,SIGN,
     &                  NOPART,TimeDep,ieaw)
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
* XINT : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      Logical TimeDep
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),DXSTST(*)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),XINT(*)
*.Local arrays
      DIMENSION ITP(36),JTP(36),KTP(36),LTP(36)
*
*.Types of DX that connects the two strings
*
*      Write(*,*)'ieaw in rsbb2a', ieaw
      ONEM = -1.0D0
      ZERO = 0.0D0
      IDXSM = STSTDX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0)  GOTO 2001
      CALL DXTYP(NDXTYP,ITP,JTP,KTP,LTP,ISEL1,ISEL3,ICEL1,ICEL3)
      DO 2000 IDXTYP = 1, NDXTYP
        ITYP = ITP(IDXTYP)
        JTYP = JTP(IDXTYP)
        KTYP = KTP(IDXTYP)
        LTYP = LTP(IDXTYP)
*. Type of intermediate strings
        CALL NEWTYP_MCLR(IGRP,ICCTP,[1],[JTYP],1,K1GRP,K1TP)
        CALL NEWTYP_MCLR(K1GRP,K1TP,[1],[LTYP],1,K2GRP,K2TP)
        IF(K2TP.LE.0) GOTO 2000
*. Symmetry of allowed Double excitation,loop over excitations
        DO 1950 IKSM = 1, NSMSX
          JLSM = SXDXSX(IKSM,IDXSM)
          IF(JLSM.EQ.0) GOTO 1950
          DO 1940 ISM = 1, NSMOB
*. Works only for D2h
            KSM = ADSXA(ISM,IKSM)
            IF(KSM.EQ.0) GOTO 1940
*. sym of intermediate strings
*
            CALL SYMCOM_MCLR(3,0,ISM,ISCSM,K1SM)
            CALL SYMCOM_MCLR(3,0,KSM,K1SM,K2SM)
* Intermediate K strings are of type K2TP and Sym K2Sm
*
            NKSTR = NSTAGTS(K2GRP,K2TP,K2SM)
            IF(NOPART.EQ.0) THEN
              NKSTREF = MIN(NKSTR,MAXK)
            ELSE
              NKSTREF = NKSTR
            END IF
            IOFF = IBTSOB(ITYP,ISM)
            KOFF = IBTSOB(KTYP,KSM)
            NI = NTSOB(ITYP,ISM)
            NK = NTSOB(KTYP,KSM)
            IF(KOFF.GT.IOFF) GOTO 1940
            IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
              IKPSM = 1
              NIK = NI*(NI+1)/2
            ELSE
              IKPSM = 0
              NIK = NI*NK
            END IF
            IF(NOPART.EQ.1) THEN
*             CALL SETVEC(SSCR,ZERO,NKSTREF*NROW*NIK)
              call dcopy_(NKSTREF*NROW*NIK,ZERO,0,SSCR,1)
            END If
            DO 1930 JSM = 1, NSMOB
              LSM = ADSXA(JSM,JLSM)
              IF(LSM.EQ.0) GOTO 1930
              JOFF = IBTSOB(JTYP,JSM)
              LOFF = IBTSOB(LTYP,LSM)
              IF(LOFF.GT.JOFF) GOTO 1930
              NJ = NTSOB(JTYP,JSM)
              NL = NTSOB(LTYP,LSM)
              IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                JLPSM = 1
                NJL = NJ*(NJ+1)/2
              ELSE
                JLPSM = 0
                NJL = NJ * NL
              END IF
              IF(NI.EQ.0.OR.NJ.EQ.0.OR.NK.EQ.0.OR.NL.EQ.0) GOTO 1930
              IFIRST = 1
*. Loop over batches of I strings
              IF(NOPART.EQ.0) THEN
                NPART = NROW/MAXI
                IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
              ELSE
                NPART = 1
              END IF
              DO 1801 IPART = 1, NPART
                IBOT = 1+(IPART-1)*MAXI
                IF(NOPART.EQ.0) THEN
                  ITOP = MIN(IBOT+MAXI-1,NROW)
                ELSE
                  ITOP = NROW
                END IF
                NIBTC = ITOP-IBOT+1
*
*.Loop over batches of intermediate strings
*
                KBOT = 1- MAXK
                KTOP = 0
 1800           CONTINUE
                  IF(NOPART.EQ.0) THEN
                    KBOT = KBOT + MAXK
                    KTOP = KTOP + MAXK
                  ELSE
                    KBOT = 1
                    KTOP = NKSTREF
                  END IF
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Generate arrays a+ja+l |kstr> for kstr in current interval
                  CALL ADADST(JTYP,JSM,JOFF,NJ,
     &                        LTYP,LSM,LOFF,NL,JLPSM,
     &                        ICCTP,ICCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,
     &                        NKBTC,NKSTREF,KEND)
                  IF(NKBTC.EQ.0) GOTO 1930
                  J = 0
                  L = 1
                  DO  IJL = 1, NJL
                    CALL NXTIJ(J,L,NJ,NL,JLPSM,NONEW)
                    LEFF = L + LOFF - 1
                    JEFF = J + JOFF - 1
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                    JLOFF = (IJL-1)*NKBTC*NIBTC+1
                    JLOFF2 = (IJL-1)*NKSTREF + 1
                    CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
     &                          I1(JLOFF2),XI1S(JLOFF2))
                  End Do
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals (ij!kl)-(il!kj)
                  IF(IFIRST.EQ.1) THEN
                    IXCHNG = 1
*                   Write(*,*)'TimeDep in rsbb2a is:',TimeDep
                   If (TimeDep) Then
                        CALL GETINT_td(XINT,ITYP,ISM,JTYP,JSM,KTYP,
     &                              KSM,LTYP,LSM,IKPSM,JLPSM,4,ieaw)
                   Else
*                         Write(*,*)'I call getint not getint_td'
                         CALL GETINT_MCLR(XINT,ITYP,ISM,JTYP,JSM,KTYP,
     &                             KSM,LTYP,LSM,IXCHNG,IKPSM,JLPSM,0,0)
                   End If
                  END IF
                  IFIRST = 0
*.and now, to the work
                  LIKB = NIBTC*NKBTC
                  IF(NOPART.EQ.1) THEN
                    FACTORC = 1.0D0
                  ELSE
                    FACTORC = 0.0D0
                  END If
                  FACTORAB = 1.0D0
                  CALL  DGEMM_('N','T',LIKB,NIK,NJL,FACTORAB,
     &                        CSCR,LIKB,XINT,NIK,
     &                        FACTORC,SSCR,LIKB)
* ============================
* Loop over ik and scatter out
* ============================
*
* Generate arrays a+i a+k !kstr>
                  IF(NOPART.EQ.0) THEN
                   CALL ADADST(ITYP,ISM,IOFF,NI,
     &                        KTYP,KSM,KOFF,NK,IKPSM,
     &                        ISCTP,ISCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,
     &                        NKBTC,NKSTREF,KEND)
                   I = 0
                   K = 1
                   DO  IK = 1, NIK
                    CALL NXTIJ(I,K,NI,NK,IKPSM,NONEW)
                    IEFF = I + IOFF - 1
                    KEFF = K + KOFF - 1
                    ISBOFF = 1+(IK-1)*NIBTC*NKBTC
                    IKOFF = (IK-1)*NKSTREF+1
                    IF( SIGN .EQ. -1.0D0)
     &                CALL DSCAL_(NKSTREf,ONEM,XI1S(IKOFF),1)
                    CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                   NKBTC,I1(IKOFF),XI1S(IKOFF))
                   End Do
                  END If
*
                IF(KEND.EQ.0.AND.NOPART.EQ.0) GOTO 1800
*. End of loop over partitionings of resolution strings
 1801         CONTINUE
 1930       CONTINUE
*. End of loop over JSM
            IF(NOPART.EQ.1) THEN
*. processing of ik dependent terms is done indepDently of
*  jl dependEnt terms :
                  CALL ADADST(ITYP,ISM,IOFF,NI,
     &                        KTYP,KSM,KOFF,NK,IKPSM,
     &                        ISCTP,ISCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,
     &                        NKBTC,NKSTREF,KEND)
                  I = 0
                  K = 1
                  DO IK = 1, NIK
                    CALL NXTIJ(I,K,NI,NK,IKPSM,NONEW)
                    IEFF = I + IOFF - 1
                    KEFF = K + KOFF - 1
                    ISBOFF = 1+(IK-1)*NIBTC*NKBTC
                    IKOFF = (IK-1)*NKSTREF+1
*. Well, someplace the minus must come in
                    IF( SIGN .EQ. -1.0D0)
     &                CALL DSCAL_(NKSTREf,ONEM,XI1S(IKOFF),1)
                    CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                   NKBTC,I1(IKOFF),XI1S(IKOFF))
                  END DO
                  END If


 1940     CONTINUE
 1950   CONTINUE
 2000 CONTINUE
*
 2001 CONTINUE
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(DXSTST)
        CALL Unused_integer_array(ITSOB)
        CALL Unused_integer(NSMDX)
      END IF
      END
