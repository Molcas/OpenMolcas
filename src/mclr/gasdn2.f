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
      SUBROUTINE GASDN2(I12,RHO1,RHO2,
     &           R,L,CB,SB,C2,ICOCOC,ISOCOC,ICSMOS,ISSMOS,
     &           ICBLTP,ISBLTP,NACOB,
     &           NSSOA,ISSOA,NSSOB,ISSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,
     &           MXPNGAS,NOBPTS,IOBPTS,
     &           MAXK,MAXI,LC,LS,
     &           CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,
     &           NGAS,NELFSPGPA,NELFSPGPB,IDC,
     &           ISOOSC,NSOOSC,ISOOSE,NSOOSE,
     &           ICOOSC,NCOOSC,ICOOSE,NCOOSE,
     &           IASOOS,IACOOS,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &           MXPOBS,IPRNT,RHO1S,
     &           LUL,LUR,PSL,PSR,RHO1P,XNATO,ieaw,n1,n2 )
*
*
* Jeppe Olsen , Winter of 1991
* GAS modificatios, August 1995
*
* =====
* Input
* =====
*
* I12    : = 1 => calculate one-electron density matrix
*          = 2 => calculate one-and two-electron density matrix
* RHO1   : Initial one-electron density matrix
* RHO2   : Initial two-electron density matrix
*
* ICOCOC : Allowed type combinations for C
* ISOCOC : Allowed type combinations for S(igma)
* ICSMOS : Symmetry array for C
* ISSMOS : Symmetry array for S
* ICBLTP : Block types for C
* ISBLTP : Block types for S
*
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* ISSOA : Offset for strings if given type and symmetry, alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* ISSOB : Offset for strings if given type and symmetry, beta strings
* NBEL  : Number of active beta electrons
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
*
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* RHO1S: Scratch array for one body
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
* The L and R vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER DXSTST(1),SXSTST(1)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB),ISOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST),ISSMOS(NSMST)
      INTEGER ICBLTP(*),ISBLTP(*)
      INTEGER NSSOA(NSMST,NOCTPA),ISSOA(NSMST,NOCTPA)
      INTEGER NSSOB(NSMST,NOCTPB),ISSOB(NSMST,NOCTPB)
      INTEGER SXSTSM(NSMSX,NSMST)
      INTEGER STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER ADSXA(MXPOBS,2*MXPOBS),ASXAD(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NELFSPGPA(3,*)
      INTEGER NELFSPGPB(3,*)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION CSCR(*),SSCR(*)
      DIMENSION ISOOSC(NOCTPA,NOCTPB,NSMST),NSOOSC(NOCTPA,NOCTPB,NSMST)
      DIMENSION ISOOSE(NOCTPA,NOCTPB,NSMST),NSOOSE(NOCTPA,NOCTPB,NSMST)
      DIMENSION ICOOSC(NOCTPA,NOCTPB,NSMST),NCOOSC(NOCTPA,NOCTPB,NSMST)
      DIMENSION ICOOSE(NOCTPA,NOCTPB,NSMST),NCOOSE(NOCTPA,NOCTPB,NSMST)
      DIMENSION IASOOS(NOCTPA,NOCTPB,NSMST),IACOOS(NOCTPA,NOCTPB,NSMST)
      DIMENSION I1(*),I2(*),XI1S(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION RHO1S(*)
      DIMENSION NOBPTS(*),IOBPTS(*)
*.
      INTEGER LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      INTEGER RASM(4),RBSM(4),RATP(4),RBTP(4),RSGN(5),RTRP(5)
      REAL * 8 INPROD_MCLR,L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*)
      DIMENSION RHO1P(*),XNATO(*)
*
      PLL=0.0D0
      PLR=0.0D0
*
* ================================
* 1 : Arrays for accessing L and R
* ================================
*.L Compact form
      NTEST=0
      CALL ZOOS(ISSMOS,ISBLTP,NSMST,ISOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,IDC,ISOOSC,NSOOSC,NSCMBC,0)
*.L Full determinant form
      CALL ZOOS(ISSMOS,ISBLTP,NSMST,ISOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,1  ,ISOOSE,NSOOSE,NSCMBE,1)
*.R compact form
      CALL ZOOS(ICSMOS,ICBLTP,NSMST,ICOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,IDC,ICOOSC,NCOOSC,NCCMBC,0)
*.R Full determinant form
      CALL ZOOS(ICSMOS,ICBLTP,NSMST,ICOCOC,NSSOA,NSSOB,
     &          NOCTPA,NOCTPB,1  ,ICOOSE,NCOOSE,NCCMBE,1)
*. Initialize loop over batches of L blocks
      ISENSM = 1
      ISENTA = 1
      ISENTB = 1
      IFRSTS = 1
* Loop over batches over L blocks
10001 CONTINUE
*. Next batch of L blocks
        ISSTSM = ISENSM
        ISSTTA = ISENTA
        ISSTTB = ISENTB
        CALL INCOOS(IDC,ISBLTP,NSOOSE,NOCTPA,NOCTPB,ISSTSM,ISSTTA,
     &              ISSTTB,NSMST,ISENSM,ISENTA,ISENTB,IASOOS,
     &              LS,ISFINI,NSBLK,IFRSTS,ISOCOC)
        IF(NSBLK.EQ.0.AND.ISFINI.NE.0) GOTO 10002
        IFRSTS = 0
*. Obtain L blocks
        IS1SM = ISSTSM
        IS1TA = ISSTTA
        IS1TB = ISSTTB
        ISOFF = 1
        DO 201 ISBLK = 1, NSBLK
          ISBSM = ISSMOS(IS1SM)
          IF(NTEST.GE.20)
     &    WRITE(6,*) ' ISBLK ISOFF ', ISBLK,ISOFF
          IF(ISOCOC(IS1TA,IS1TB).EQ.1)
     &    CALL GSTTBL_MCLR(L,SB(ISOFF),IS1TA,IS1SM,IS1TB,ISBSM,ISOCOC,
     &                NOCTPA,NOCTPB,NSSOA,NSSOB,PSL,ISOOSC,IDC,
     &                PLL,LUL,C2)
*    &                PSL,LUL,C2,NSMST)
          ISOFF = ISOFF + NSOOSE(IS1TA,IS1TB,IS1SM)
          IF(ISBLK.NE.NSBLK) THEN
            CALL NXTBLK_MCLR(IS1TA,IS1TB,IS1SM,NOCTPA,NOCTPB,NSMST,
     &      ISBLTP,IDC,NONEWS,ISOCOC)
          END IF
  201   CONTINUE
*. Initialize loop over blocks over L vector
        ICENSM = 1
        ICENTA = 1
        ICENTB = 1
*. Loop over blocks of R vector
        IFRSTC = 1
 9001   CONTINUE
        IF(NTEST.GE.20)
     &  write(6,*) ' >>> next batch of R blocks '
          ICSTSM = ICENSM
          ICSTTA = ICENTA
          ICSTTB = ICENTB
*
          CALL INCOOS(IDC,ICBLTP,NCOOSE,NOCTPA,NOCTPB,ICSTSM,ICSTTA,
     &         ICSTTB,NSMST,ICENSM,ICENTA,ICENTB,IACOOS,
     &         LC,IFINIC,NCBLK,IFRSTC,ICOCOC)
*. If no more R blocks goto next batch of L blocks
          IF(NCBLK.EQ.0.AND.IFINIC.NE.0) GOTO 10001
          IFRSTC = 0
*. Read L blocks into core
          IC1SM = ICSTSM
          IC1TA = ICSTTA
          IC1TB = ICSTTB
          ICOFF = 1
          DO 200 ICBLK = 1, NCBLK
            ICBSM = ICSMOS(IC1SM)
            IF(ICOCOC(IC1TA,IC1TB).EQ.1)
     &      CALL GSTTBL_MCLR(R,CB(ICOFF),IC1TA,IC1SM,IC1TB,ICBSM,ICOCOC,
     &                  NOCTPA,NOCTPB,NSSOA,NSSOB,PSR,ICOOSC,IDC,
     &                  PLR,LUR,C2)
*    &                  PCL,LUR,C2,NSMST)
            ICOFF = ICOFF + NCOOSE(IC1TA,IC1TB,IC1SM)
            IF(ICBLK.NE.NCBLK) THEN
              CALL NXTBLK_MCLR(IC1TA,IC1TB,IC1SM,NOCTPA,NOCTPB,NSMST,
     &        ICBLTP,IDC,NONEWC,ICOCOC)
            END IF
  200     CONTINUE
*. Loop over L and R blocks in core and obtain  contribution from
* given L and R blocks
          ISOFF = 1
          IASM = ISSTSM
          IATP = ISSTTA
          IBTP = ISSTTB
          DO 10000 ISBLK = 1, NSBLK
            IBSM = ISSMOS(IASM)
            NIA = NSSOA(IASM,IATP)
            NIB = NSSOB(IBSM,IBTP)
*. Possible permutations of L blocks
            CALL PRMBLK(IDC,ISTRFL,IASM,IBSM,IATP,IBTP,PSL,PLR,
     &           LATP,LBTP,LASM,LBSM,LSGN,LTRP,NLPERM)
            DO 9999 ILPERM = 1, NLPERM
              IIASM = LASM(ILPERM)
              IIBSM = LBSM(ILPERM)
              IIATP = LATP(ILPERM)
              IIBTP = LBTP(ILPERM)
              NIIA = NSSOA(IIASM,IIATP)
              NIIB = NSSOB(IIBSM,IIBTP)
*
              IF(LTRP(ILPERM).EQ.1) THEN
                LROW = NSSOA(LASM(ILPERM-1),LATP(ILPERM-1))
                LCOL = NSSOB(LBSM(ILPERM-1),LBTP(ILPERM-1))
                CALL TRPMT3(SB(ISOFF),LROW,LCOL,C2)
                CALL COPVEC(C2,SB(ISOFF),LROW*LCOL)
               END IF
               IF(LSGN(ILPERM).EQ.-1)
     &         CALL SCALVE(SB(ISOFF),-1.0D0,NIA*NIB)

               JASM = ICSTSM
               JATP = ICSTTA
               JBTP = ICSTTB
               ICOFF = 1
               DO 9000 ICBLK = 1, NCBLK
                 JBSM = ICSMOS(JASM)
                 NJA = NSSOA(JASM,JATP)
                 NJB = NSSOB(JBSM,JBTP)
                 XNORM2 =
     &           INPROD_MCLR(CB(ICOFF),CB(ICOFF),NJA*NJB)!,NJA*NJB)
                 IF(NIA*NIB*NJA*NJB.NE.0
     &             .AND.ISOCOC(IATP,IBTP).EQ.1
     &             .AND.ICOCOC(JATP,JBTP).EQ.1
     &             .AND.XNORM2.NE.0.0D0) THEN
*. Possible permutations of this block
                    CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,
     &                   PSR,PLR,RATP,RBTP,RASM,RBSM,RSGN,RTRP,
     &                   NRPERM)
                    DO 8999 IRPERM = 1, NRPERM
                      IF(RTRP(IRPERM).EQ.1) THEN
                        LROW = NSSOA(RASM(IRPERM-1),RATP(IRPERM-1))
                        LCOL = NSSOB(RBSM(IRPERM-1),RBTP(IRPERM-1))
                        CALL TRPMT3(CB(ICOFF),LROW,LCOL,C2)
                        CALL COPVEC(C2,CB(ICOFF),LROW*LCOL)
                      END IF
                      IF(RSGN(IRPERM).EQ.-1)
     &                CALL SCALVE(CB(ICOFF),-1.0D0,NJA*NJB)
                      JJASM = RASM(IRPERM)
                      JJBSM = RBSM(IRPERM)
                      JJATP = RATP(IRPERM)
                      JJBTP = RBTP(IRPERM)
                      NJJA = NSSOA(JJASM,JJATP)
                      NJJB = NSSOB(JJBSM,JJBTP)
                      CALL GSDNBB2(I12,RHO1,RHO2,
     &                     IIASM,IIATP,IIBSM,IIBTP,
     &                     JJASM,JJATP,JJBSM,JJBTP,NGAS,
     &                     NELFSPGPA(1,IOCTPA-1+IIATP),
     &                     NELFSPGPB(1,IOCTPB-1+IIBTP),
     &                     NELFSPGPA(1,IOCTPA-1+JJATP),
     &                     NELFSPGPB(1,IOCTPB-1+JJBTP),
     &                     NAEL,NBEL,IAGRP,IBGRP,
     &                     SB(ISOFF),CB(ICOFF),C2,
     &                     ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                     MXPNGAS,NOBPTS,IOBPTS,MAXI,MAXK,
     &                     SSCR,CSCR,
     &                     I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                     X,NSMOB,NSMST,NSMSX,NSMDX,
     &                     NIIA,NIIB,NJJA,NJJB,MXPOBS,
     &                     IPRNT,NACOB,RHO1S,ieaw,n1,n2)
 8999                   CONTINUE
*. Transpose or scale R block to restore order ??
                  IF(RTRP(NRPERM+1).EQ.1) THEN
                    CALL TRPMT3(CB(ICOFF),NJB,NJA,C2)
                    CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
                  END IF
                  IF(RSGN(NRPERM+1).EQ.-1)
     &            CALL SCALVE(CB(ICOFF),-1.0D0,NJA*NJB)
*
              END IF
              ICOFF = ICOFF + NCOOSE(JATP,JBTP,JASM)
*. NeXt C block
              IF(ICBLK.NE.NCBLK) THEN
                CALL NXTBLK_MCLR(JATP,JBTP,JASM,NOCTPA,NOCTPB,NSMST,
     &          ICBLTP,IDC,NONEWC,ICOCOC)
              END IF
 9000       CONTINUE
*. End of loop over R blocks in Batch
 9999     CONTINUE
*. Transpose or scale L block to restore order ??
          IF(LTRP(NLPERM+1).EQ.1) THEN
            CALL TRPMT3(SB(ISOFF),NIB,NIA,C2)
            CALL COPVEC(C2,SB(ISOFF),NIA*NIB)
          END IF
          IF(LSGN(NLPERM+1).EQ.-1)
     &    CALL SCALVE(SB(ISOFF),-1.0D0,NIA*NIB)
*. NeXt L block
               ISOFF = ISOFF + NSOOSE(IATP,IBTP,IASM)
            IF(ISBLK.NE.NSBLK) THEN
               CALL NXTBLK_MCLR(IATP,IBTP,IASM,NOCTPA,NOCTPB,NSMST,
     &         ISBLTP,IDC,NONEWS,ISOCOC)
             END IF
10000     CONTINUE
*. End of loop over L blocks in batch
*. End of loop over batches of R blocks
        IF(IFINIC.EQ.0) GOTO 9001
      IF(ISFINI.EQ.0) GOTO 10001
*. End of loop over batches of L blocks
10002     CONTINUE
*
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(ISSOA)
        CALL Unused_integer_array(ISSOB)
        CALL Unused_integer_array(SXSTSM)
        CALL Unused_integer_array(ASXAD)
        CALL Unused_real_array(RHO1P)
        CALL Unused_real_array(XNATO)
      END IF
      END
