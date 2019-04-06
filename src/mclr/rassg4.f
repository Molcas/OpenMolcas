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
*               1996, Anders Bernhardsson                              *
************************************************************************
      SUBROUTINE RASSG4(C,S,CB,SB,C2,ICOCOC,ISOCOC,ICSMOS,ISSMOS,
     &                  ICBLTP,ISBLTP,
     &                  NORB1,NORB2,NORB3,NACOB,
     &                  NSSOA,ISSOA,NSSOB,ISSOB,
     &                  NAEL,IAGRP,NBEL,IBGRP,NOCTPA,NOCTPB,
     &                  NSMST,NSMOB,NSMSX,NSMDX,NTSOB,IBTSOB,ITSOB,
     &                  MAXIJ,MAXK,MAXI,ICSMOD,IINMOD,LI,LC,LS,
     &                  XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &                  SXDXSX,ADSXA,ASXAD,
     &                  IAEL1,IAEL3,
     &                  IBEL1,IBEL3,IDC,
     &                  ISOOSC,NSOOSC,ISOOSE,NSOOSE,
     &                  ICOOSC,NCOOSC,ICOOSE,NCOOSE,
     &                  IASOOS,IACOOS,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  IDOH2,ISTRFL,PS,IPRNT,LUC,LUHC,IST,
     &                  CJRES,SIRES,NOPARt,TimeDep)
      IMPLICIT REAL*8(A-H,O-Z)
*
*      LOOP OVER SIGMA AND C VECTOR
*
#include "detdim.fh"
* Jeppe Olsen , Winter of 1991
* small modifications by eaw 96
*
* =====
* Input
* =====
*
* ICOCOC : Allowed type combinations for C
* ISOCOC : Allowed type combinations for S(igma)
* ICSMOS : Symmetry array for C
* ISSMOS : Symmetry array for S
* ICBLTP : Block types for C
* ISBLTP : Block types for S
*
* NORB1(2,3) : Number of orbitals in RAS1(2,3)
* NACOB : Number of active orbitals
* H     : Active one-body Hamiltonian with core contributions
* C     : CI vector
* CB    : Array able to hold largest STT block of C
* NSSOA : Number of strings per type and symmetry for alpha strings
* ISSOA : Offset for strings if given type and symmetry, alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* ISSOB : Offset for strings if given type and symmetry, beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* ITSOB : Orbitals of given type and symmetry
* IBTSOB: Offset for ITSOB
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* ICSMOD : 1 => Single symmetry blocks of C and S are treated
*               simultaneously
* ICSMOD : 2 => Single symmetry-occ-occ blocks of C and S are treated
*               simultaneously
* IINMOD :
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector (space for res. matr.)
* SSCR : Scratch array for S vector (space for res. matr.)
*
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* IST :  = 1 => Singlet operator
*        = 2 => Triplet operator\
*
* IDOH2 : = 1 => both one and two particle parts
*         = 0 => only one-electron operator
* A triplet one electron operator is defined as E(aa)-E(bb)
* A triplet two-electron operator is defined as (E(aa)+E(bb))(E(aa)-E(bb))
*
      Logical TimeDep
*.General input
      INTEGER SXSTST(1), DXSTST(1) ! HMMMMMM
      INTEGER ICOCOC(NOCTPA,NOCTPB),ISOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST),ISSMOS(NSMST)
      INTEGER ICBLTP(NSMST),ISBLTP(NSMST)
      INTEGER NSSOA(NOCTPA,nsmst),ISSOA(NOCTPA,nsmst)
      INTEGER NSSOB(NOCTPB,nsmst),ISSOB(NOCTPB,nsmst)
      INTEGER SXSTSM(NSMSX,NSMST)
      INTEGER STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER ADSXA(MXPOBS,2*MXPOBS),ASXAD(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NTSOB(3,NSMOB),IBTSOB(3,NSMOB),ITSOB(mxporb)
      INTEGER IAEL1(*),IAEL3(*)
      INTEGER IBEL1(*),IBEL3(*)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      INTEGER ISOOSC(NOCTPA,NOCTPB,NSMST),NSOOSC(NOCTPA,NOCTPB,NSMST)
      INTEGER ISOOSE(NOCTPA,NOCTPB,NSMST),NSOOSE(NOCTPA,NOCTPB,NSMST)
      INTEGER ICOOSC(NOCTPA,NOCTPB,NSMST),NCOOSC(NOCTPA,NOCTPB,NSMST)
      INTEGER ICOOSE(NOCTPA,NOCTPB,NSMST),NCOOSE(NOCTPA,NOCTPB,NSMST)
      INTEGER IASOOS(NOCTPA,NOCTPB,NSMST),IACOOS(NOCTPA,NOCTPB,NSMST)
      DIMENSION I1(MAXK,*),I2(MAXK,*),XI1S(MAXK,*),XI2S(MAXK,*),
     &          I3(MAXK,*),I4(MAXK,*),XI3S(MAXK,*),XI4S(MAXK,*)
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
*.
      DIMENSION C(*),S(*),ISTRFL(*)
*
      ZERO = 0.0D0
      PL=Zero
* ================================
* 1 : Arrays for accessing C and S
* ================================
*

********************************************************************
*
*.Sigma, compact form
      CALL ZOOS(ISSMOS,ISBLTP,NSMST,ISOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,IDC,ISOOSC,NSOOSC,NSCMBC,0)
*. Sigma with expanded diagonal blocks
      CALL ZOOS(ISSMOS,ISBLTP,NSMST,ISOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,IDC,ISOOSE,NSOOSE,NSCMBE,1)
*.C, compact form
      CALL ZOOS(ICSMOS,ICBLTP,NSMST,ISOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,IDC,ICOOSC,NCOOSC,NCCMBC,0)
*.C, Full determinant form
      CALL ZOOS(ICSMOS,ICBLTP,NSMST,ISOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,1  ,ICOOSE,NCOOSE,NCCMBE,1)
*
*********************************************************************
*

*. Initialize loop over batches of sigma blocks
      ISENSM = 1
      ISENTA = 1
      ISENTB = 1
      IFRSTS = 1
* Loop over batches over sigma blocks
      IF(LUHC.GT.0) REWIND LUHC

10001 CONTINUE
*. Next batch of sigma blocks
        ISSTSM = ISENSM
        ISSTTA = ISENTA
        ISSTTB = ISENTB
        CALL INCOOS(IDC,ISBLTP,NSOOSE,NOCTPA,NOCTPB,ISSTSM,ISSTTA,
     &              ISSTTB,NSMST,ISENSM,ISENTA,ISENTB,IASOOS,
     &              LS,ISFINI,NSBLK,IFRSTS,ISOCOC)
        IF(NSBLK.EQ.0.AND.ISFINI.NE.0) GOTO 10002
        IFRSTS = 0
*. Initialize sigma blocks
        IS1SM = ISSTSM
        IS1TA = ISSTTA
        IS1TB = ISSTTB
        LSBLK = 0
        DO ISBLK = 1, NSBLK
          LSBLK = LSBLK + NSOOSE(IS1TA,IS1TB,IS1SM)
          IF(ISBLK.NE.NSBLK) THEN
            CALL NXTBLK_MCLR(IS1TA,IS1TB,IS1SM,NOCTPA,NOCTPB,NSMST,
     &      ISBLTP,IDC,NONEWS,ISOCOC)
          END IF
        End Do
*       CALL SETVEC(SB,ZERO ,LSBLK)
        call dcopy_(LSBLK,[ZERO],0,SB,1)
*. Initialize loop over blocks over C vector
        ICENSM = 1
        ICENTA = 1
        ICENTB = 1
        IF(LUC.GT.0) REWIND LUC
*. Loop over blocks of C vector
        IFRSTC = 1
 9001   CONTINUE
          ICSTSM = ICENSM
          ICSTTA = ICENTA
          ICSTTB = ICENTB
*
          CALL INCOOS(IDC,ICBLTP,NCOOSE,NOCTPA,NOCTPB,ICSTSM,ICSTTA,
     &         ICSTTB,NSMST,ICENSM,ICENTA,ICENTB,IACOOS,
     &         LC,IFINIC,NCBLK,IFRSTC,ICOCOC)

*. If no more C blocks goto next batch of sigma blocks
          IF(NCBLK.EQ.0.AND.IFINIC.NE.0) GOTO 10001
          IFRSTC = 0
*. Read C blocks into core
          IC1SM = ICSTSM ! Symmetry alpha
          IC1TA = ICSTTA ! Type alpha string
          IC1TB = ICSTTB ! Type Beta string
          ICOFF = 1
          DO ICBLK = 1, NCBLK
            ICBSM = ICSMOS(IC1SM) ! Symmetry Beta string
            IF(ICOCOC(IC1TA,IC1TB).EQ.1) THEN
*
*              C CI vector in
*              CB Block of CI vector out
*

             CALL GSTTBL_MCLR(C,CB(ICOFF),IC1TA,IC1SM,IC1TB,ICBSM,
     &                  ICOCOC,NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                  PL,LUC,C2)
            END IF

            ICOFF = ICOFF + NCOOSE(IC1TA,IC1TB,IC1SM)
            IF(ICBLK.NE.NCBLK) THEN
              CALL NXTBLK_MCLR(IC1TA,IC1TB,IC1SM,NOCTPA,NOCTPB,NSMST,
     &        ICBLTP,IDC,NONEWC,ICOCOC)
            END IF
          End Do
*
*********************************************************************
*
*. Loop over sigma and C blocks in core and obtain  contribution from
* given C block to given S block
          ISOFF = 1
          IASM = ISSTSM
          IATP = ISSTTA
          IBTP = ISSTTB
          DO 10000 ISBLK = 1, NSBLK
           call xflush(6)
           IBSM = ISSMOS(IASM)
           NIA = NSSOA(IATP,IASM)
           NIB = NSSOB(IBTP,IBSM)
           IF(NIA.NE.0.AND.NIB.NE.0) THEN
            JASM = ICSTSM
            JATP = ICSTTA
            JBTP = ICSTTB
            ICOFF = 1
            DO 9000 ICBLK = 1, NCBLK
              call xflush(6)
              JBSM = ICSMOS(JASM)
              NJA = NSSOA(JATP,JASM)
              NJB = NSSOB(JBTP,JBSM)
              XNORM2 = DDot_(NJA*NJB,CB(ICOFF),1,CB(ICOFF),1)
              IF(NIA.NE.0.AND.NIB.NE.0.AND.NJA.NE.0.AND.NJB.NE.0 .AND.
     &        ISOCOC(IATP,IBTP).EQ.1.AND.ICOCOC(JATP,JBTP).EQ.1.AND.
     &         XNORM2.NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
!                write(*,*)"Other symmetry blocks that can be obtained"
!                call xflush(6)

                CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                      LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
C    &                      TimeDep)
                DO 8765 IPERM = 1, NPERM
                  LLASM = LASM(IPERM)
                  LLBSM = LBSM(IPERM)
                  LLATP = LATP(IPERM)
                  LLBTP = LBTP(IPERM)
                  NLLA = NSSOA(LLATP,LLASM)
                  NLLB = NSSOB(LLBTP,LLBSM)
                  IF(LTRP(IPERM).EQ.1) THEN
                    LROW = NSSOA(LATP(IPERM-1),LASM(IPERM-1))
                    LCOL = NSSOB(LBTP(IPERM-1),LBSM(IPERM-1))
                    CALL TRPMAT(CB(ICOFF),LROW,LCOL,C2)
                    call dcopy_(LROW*LCOL,C2,1,CB(iCOFF),1)
                  END IF
                  IF(LSGN(IPERM).EQ.-1)
     &            CALL DSCAL_(LROW*LCOL,-1.0d0,CB(ICOFF),1)
*
*
*    Generation of contribution to sigma block
*    from given CI block
*
!                 Write(*,*)'TimeDep in rassg4',TimeDep
!                 call xflush(6)
                 If (TimeDep) Then
*                      Write(*,*)'I call rssbcbn_td'
                      CALL RSSBCBN_td(IASM,IATP,IBSM,IBTP,
     &                LLASM,LLATP,LLBSM,LLBTP,
     &                IAEL1(IATP),IAEL3(IATP),
     &                IBEL1(IBTP),IBEL3(IBTP),
     &                IAEL1(LLATP),IAEL3(LLATP),
     &                IBEL1(LLBTP),IBEL3(LLBTP),
     &                NAEL,NBEL,
     &                IAGRP,IBGRP,
     &                SB(ISOFF),CB(ICOFF),IDOH2,
     &                ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,
     &                C2,NSMOB,NSMST,NSMSX,NSMDX,NIA,NIB,NLLA,NLLB,
     &                MXPOBS,IPRNT,IST,CJRES,SIRES,NOPART,TimeDep)
*
!                      Call RECPRT('SSCR in rassg4',' ',SSCR,5,1) !yma
!                      call xflush(6)
                 Else
                      CALL RSSBCBN_MCLR(IASM,IATP,IBSM,IBTP,
     &                LLASM,LLATP,LLBSM,LLBTP,
     &                IAEL1(IATP),IAEL3(IATP),
     &                IBEL1(IBTP),IBEL3(IBTP),
     &                IAEL1(LLATP),IAEL3(LLATP),
     &                IBEL1(LLBTP),IBEL3(LLBTP),
     &                NAEL,NBEL,
     &                IAGRP,IBGRP,
     &                SB(ISOFF),CB(ICOFF),IDOH2,
     &                ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                NTSOB,IBTSOB,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,
     &                C2,NSMOB,NSMST,NSMSX,NSMDX,NIA,NIB,NLLA,NLLB,
     &                MXPOBS,IPRNT,IST,CJRES,SIRES,NOPART,TimeDep)
                 End If
*
 8765           CONTINUE
*. Transpose or scale to restore order ??
                  IF(LTRP(NPERM+1).EQ.1) THEN
                    CALL TRPMAT(CB(ICOFF),NJB,NJA,C2)
                    call dcopy_(NJA*NJB,C2,1,CB(ICOFF),1)
                  END IF
                  IF(LSGN(NPERM+1).EQ.-1)
     &            CALL DSCAL_(NJA*NJB,-1.0d0,CB(ICOFF),1)
*
              END IF
              ICOFF = ICOFF + NCOOSE(JATP,JBTP,JASM)
*. NeXt C block
              IF(ICBLK.NE.NCBLK) THEN
                CALL NXTBLK_MCLR(JATP,JBTP,JASM,NOCTPA,NOCTPB,NSMST,
     &          ICBLTP,IDC,NONEWC,ICOCOC)
              END IF
 9000       CONTINUE
*. End of loop over C blocks in Batch
*. NeXt S block
           END IF
           ISOFF = ISOFF + NSOOSE(IATP,IBTP,IASM)
           IF(ISBLK.NE.NSBLK) THEN
               CALL NXTBLK_MCLR(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSMST,
     &         ISBLTP,IDC,NONEWS,ISOCOC)
            END IF
10000     CONTINUE
********************************************************************
*. End of loop over S blocks in batch
*. End of loop over batches of C blocks

          IF(IFINIC.EQ.0) GOTO 9001
*. Transfer S block to permanent storage
C
*          Call RECPRT('SB in rassg4',' ',S,5,1)
C
          I1ASM = ISSTSM
          I1TA  = ISSTTA
          I1TB  = ISSTTB
          IOFF  = 1
          DO ISBLK = 1, NSBLK
           I1BSM = ISSMOS(I1ASM)
           IF(ISOCOC(I1TA,I1TB).EQ.1) THEN
            CALL PSTTBL_MCLR(S,SB(IOFF),I1TA,I1ASM,I1TB,I1BSM,ISOCOC,
     &                  NOCTPA,NOCTPB,NSSOA,NSSOB,PS,
     &                  ISOOSC,2,IDC,LUHC,C2)
           END IF
           IOFF = IOFF + NSOOSE(I1TA,I1TB,I1ASM)
           IF(ISBLK.NE.NSBLK) THEN
            CALL NXTBLK_MCLR(I1TA,I1TB,I1ASM,NOCTPA,NOCTPB,NSMST,
     &      ISBLTP,IDC,NONEWS,ISOCOC)
           END IF
          End Do
      IF(ISFINI.EQ.0) GOTO 10001
*. End of loop over batches of S blocks
10002 CONTINUE
********************************************************************
      IF(LUHC.GT.0) CALL ITODS([-1],1,LBLK,LUHC)

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(NORB1)
        CALL Unused_integer(NORB2)
        CALL Unused_integer(NORB3)
        CALL Unused_integer(NACOB)
        CALL Unused_integer_array(ISSOA)
        CALL Unused_integer_array(ISSOB)
        CALL Unused_integer(MAXIJ)
        CALL Unused_integer(ICSMOD)
        CALL Unused_integer(IINMOD)
        CALL Unused_integer(LI)
        CALL Unused_integer_array(SXSTSM)
        CALL Unused_integer_array(ASXAD)
      END IF
      END
