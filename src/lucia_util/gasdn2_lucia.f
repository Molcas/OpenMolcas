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
* Copyright (C) 1991,1995,1997-1999, Jeppe Olsen                       *
************************************************************************
      SUBROUTINE GASDN2_LUCIA(    I12,   RHO1,   RHO2,  RHO2S,  RHO2A,
     &                              L,      R,     CB,     SB,     C2,
     &                         ICOCOC, ISOCOC, ICSMOS, ISSMOS, ICBLTP,
     &                         ISBLTP,  NACOB,  NSSOA,  ISSOA,  NSSOB,
     &                          ISSOB,   NAEL,  IAGRP,   NBEL,  IBGRP,
*
     &                         IOCTPA, IOCTPB, NOCTPA, NOCTPB,  NSMST,
     &                          NSMOB,  NSMSX,  NSMDX,MXPNGAS, NOBPTS,
     &                         IOBPTS,   MAXK,   MAXI,     LC,     LS,
     &                           CSCR,   SSCR, SXSTSM, STSTSX, STSTDX,
     &                         SXDXSX,  ADSXA,  ASXAD,   NGAS,NELFSPGP,
*
     &                            IDC,     I1,   XI1S,     I2,   XI2S,
     &                             I3,   XI3S,     I4,   XI4S,      X,
     &                         MXPOBS,  IPRNT,  RHO1S,    LUL,    LUR,
     &                            PSL,    PSR,  RHO1P, XNATO ,NBATCHL,
     &                          LBATL, LEBATL, I1BATL,IBLOCKL,NBATCHR,
*
     &                          LBATR, LEBATR, I1BATR,IBLOCKR,ICONSPA,
     &                        ICONSPB,
     &                        SCLFAC_L,SCLFAC_R,S2_TERM1,IUSE_PH,IPHGAS,
     &                        IDOSRHO1, SRHO1, IPACK)
*
* SUBROUTINE GASDN2_LUCIA --> 89
*
*
*
* Jeppe Olsen , Winter of 1991
* GAS modificatios, August 1995
*
* Table driven, June 97
*
* Last revision : Jan. 98 (IUSE_PH,IPHGAS added)
*                 Jan. 99 (IDOSRHO1,SRHO1 added)
*
* =====
* Input
* =====
*
* I12    : = 1 => calculate one-electrondensity matrix
*          = 2 => calculate one-and two- electrondensity matrix
* RHO1   : Initial one-electron density matrix
* RHO2   : Initial two-electron density matrix
* RHO2S  : Initial symmetric two-electron density matrix
* RHO2A  : Initial anti-symmetric two-electron density matrix
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
* IPACK : Logical: If true calculate symmetry packed densities (i.e.
*         RHO2S and RHO2A instead of RHO2
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
#include "para_info.fh"
#include "io_util.fh"
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
      INTEGER NOBPTS(MXPNGAS,NSMOB),IOBPTS(MXPNGAS,NSMOB)
      INTEGER NELFSPGP(MXPNGAS,*)
      LOGICAL IPACK
      DIMENSION IPHGAS(*),SRHO1(*)
*. Info on batches and blocks
      INTEGER  LBATL(NBATCHL),LEBATL(NBATCHL),I1BATL(NBATCHL),
     &         IBLOCKL(8,*)
      INTEGER  LBATR(NBATCHR),LEBATR(NBATCHR),I1BATR(NBATCHR),
     &         IBLOCKR(8,*)
*. Interaction between supergroups
      INTEGER ICONSPA(NOCTPA,NOCTPA),ICONSPB(NOCTPB,NOCTPB)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),XI1S(*),XI2S(*),I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION RHO1S(*)
      DIMENSION SCLFAC_L(*),SCLFAC_R(*)
      DIMENSION ICOOSC(NOCTPA,NOCTPB),ISOOSC(NOCTPA,NOCTPB)
*.
      INTEGER LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      INTEGER RASM(4),RBSM(4),RATP(4),RBTP(4),RSGN(5),RTRP(5)
C-jwk-cleanup      REAL * 8 INPROD,L
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*),RHO2S(*),RHO2A(*)
      DIMENSION RHO1P(*),XNATO(*)

      DIMENSION ISTRFL(1),LBL(1),IDUMMY(1)
      INTEGER SXSTST(1),DXSTST(1)
* Some dummy initializations
      INTERACT = 0 ! jwk-cleanup
*

      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
         WRITE(6,*) ' ================='
         WRITE(6,*) ' GASDN2 speaking :'
         WRITE(6,*) ' ================='
         WRITE(6,*)
         WRITE(6,*) ' NACOB,MAXK,NGAS,IDC,MXPOBS',
     &        NACOB,MAXK,NGAS,IDC,MXPOBS
         WRITE(6,*) ' LUL, LUR ', LUL,LUR
      END IF
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' Initial L vector '
         IF(LUL.EQ.0) THEN
            CALL WRTRS2(       L,  ISSMOS,  ISBLTP,  ISOCOC,  NOCTPA,
     &                    NOCTPB,   NSSOA,   NSSOB,   NSMST)
         ELSE
            CALL WRTVCD(L,LUL,1,-1)
         END IF
         WRITE(6,*) ' Initial R vector '
         IF(LUR.EQ.0) THEN
            CALL WRTRS2(       R,  ICSMOS,  ICBLTP,  ICOCOC,  NOCTPA,
     &                    NOCTPB,   NSSOA,   NSSOB,   NSMST)
         ELSE
            CALL WRTVCD(R,LUR,1,-1)
         END IF
      END IF
* Loop over batches over L blocks
      IF(LUL.NE.0) IDISK(LUL)=0
      DO 10001 IBATCHL = 1, NBATCHL
*. Obtain L blocks
         NBLKL = LBATL(IBATCHL)
        IF(NTEST.GE.200)
     &        WRITE(6,*) ' Left batch, number of blocks',IBATCHL,NBLKL
        DO IIL  = 1,NBLKL
           IL  = I1BATL(IBATCHL)-1+IIL
           IATP = IBLOCKL(1,IL)
           IBTP = IBLOCKL(2,IL)
           IASM = IBLOCKL(3,IL)
           IBSM = IBLOCKL(4,IL)
           IOFF = IBLOCKL(5,IL)
           IF(NTEST.GE.200)
     &          WRITE(6,*) 'IATP IBTP IASM IBSM',IATP,IBTP,IASM,IBSM
           ISCALE = 0
           IF(NTEST.GE.200)
     &          WRITE(6,*) 'IOFF ',IOFF
           CALL GSTTBL(       L,SB(IOFF),    IATP,    IASM,    IBTP,
     &                     IBSM,  ISOCOC,  NOCTPA,  NOCTPB,   NSSOA,
     &                    NSSOB,     PSL,  ISOOSC,     IDC,     PSL,
     &                      LUL,      C2,   NSMST,  ISCALE,SCLFAC_L(IL))
        END DO
*. Loop over batches  of R vector
        IF(LUR.NE.0) IDISK(LUR)=0
        DO 9001 IBATCHR = 1, NBATCHR
*. Read R blocks into core
           NBLKR = LBATR(IBATCHR)
           IF(NTEST.GE.200) WRITE(6,*) ' Right batch, number of blocks',
     &          IBATCHR,NBLKR
           DO IIR  = 1,NBLKR
              IR  = I1BATR(IBATCHR)-1+IIR
              JATP = IBLOCKR(1,IR)
              JBTP = IBLOCKR(2,IR)
              JASM = IBLOCKR(3,IR)
              JBSM = IBLOCKR(4,IR)
              JOFF = IBLOCKR(5,IR)
              IF(NTEST.GE.200) WRITE(6,*) ' JATP JBTP JASM JBSM ',
     &             JATP,JBTP,JASM,JBSM
*. Read R blocks into core
*
*. Only blocks interacting with current batch of L are read in
*. Loop over L  blocks in batch
              DO IIL = 1, NBLKL
                 IL  = I1BATL(IBATCHL)-1+IIL
                 IATP = IBLOCKL(1,IL)
                 IBTP = IBLOCKL(2,IL)
                 IASM = IBLOCKL(3,IL)
                 IBSM = IBLOCKL(4,IL)
*. Well, permutations of L blocks
                 PS=1.0D0
                 PL=1.0D0
                 CALL PRMBLK(    IDC, ISTRFL,   IASM,   IBSM,   IATP,
     &                          IBTP,     PS,     PL,   LATP,   LBTP,
     &                          LASM,   LBSM,   LSGN,   LTRP,  NPERM)
                 DO IPERM = 1, NPERM
                    IIASM = LASM(IPERM)
                    IIBSM = LBSM(IPERM)
                    IIATP = LATP(IPERM)
                    IIBTP = LBTP(IPERM)

                    IAEXC = ICONSPA(IIATP,JATP)
                    IBEXC = ICONSPB(IIBTP,JBTP)
                    IF(IAEXC.EQ.0.AND.IIASM.NE.JASM) IAEXC = 1
                    IF(IBEXC.EQ.0.AND.IIBSM.NE.JBSM) IBEXC = 1
                    IABEXC = IAEXC + IBEXC
                    IF(IABEXC.LE.I12) THEN
                       INTERACT = 1
                    END IF
                 END DO
              END DO
*.            ^ End of checking whether C-block is needed
              ISCALE = 0
              IF(INTERACT.EQ.1) THEN
                 ISCALE = 0
                 CALL GSTTBL(      R,CB(JOFF),  JATP,  JASM,  JBTP,
     &                          JBSM, ICOCOC, NOCTPA, NOCTPB,  NSSOA,
     &                         NSSOB,    PSR, ICOOSC,    IDC,    PCL,
     &                           LUR,   C2,  NSMST, ISCALE,SCLFAC_R(IR))
              ELSE
C             WRITE(6,*) ' TTSS for C block skipped  '
C             CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
                 CALL IDAFILE(LUR,2,LBL,1,IDISK(LUR))
                 CALL IDAFILE(LUR,2,IDUMMY,1,IDISK(LUR))
                 CALL SKPRCD2(LBL(1),-1,LUR)
                 SCLFAC_R(IR) = 0.0D0
              END IF
*
*
              IF(NTEST.GE.100) THEN
                 IF(INTERACT.EQ.1) THEN
                    WRITE(6,*) ' TTSS for C block read in  '
                    CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
                 ELSE
                    WRITE(6,*) ' TTSS for C block skipped  '
                    CALL IWRTMA(IBLOCKR(1,IR),4,1,4,1)
                 END IF
              END IF
           END DO

*. Loop over L and R blocks in batches and obtain  contribution from
* given L and R blocks
           PLR=1.D0
           DO 10000 IIL = 1, NBLKL
              IL  = I1BATL(IBATCHL)-1+IIL
              IF(SCLFAC_L(IL).NE.0.0D0) THEN
                 IATP = IBLOCKL(1,IL)
                 IBTP = IBLOCKL(2,IL)
                 IASM = IBLOCKL(3,IL)
                 IBSM = IBLOCKL(4,IL)
                 IOFF = IBLOCKL(5,IL)
*
                 NIA = NSSOA(IASM,IATP)
                 NIB = NSSOB(IBSM,IBTP)
*. Possible permutations of L blocks
                 CALL PRMBLK(    IDC, ISTRFL,   IASM,   IBSM,   IATP,
     &                          IBTP,    PSL,    PLR,   LATP,   LBTP,
     &                          LASM,   LBSM,   LSGN,   LTRP, NLPERM)
                 DO 9999 ILPERM = 1, NLPERM
C             write(6,*) ' Loop 9999 ILPERM = ', ILPERM
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
                       CALL TRPMT3(SB(IOFF),LROW,LCOL,C2)
                       CALL COPVEC(C2,SB(IOFF),LROW*LCOL)
                    END IF
                    IF(LSGN(ILPERM).EQ.-1)
     &                   CALL SCALVE(SB(IOFF),-1.0D0,NIA*NIB)

                    DO 9000 IIR = 1, NBLKR
                       IR  = I1BATR(IBATCHR)-1+IIR
                       IF(SCLFAC_R(IR).NE.0.0D0) THEN
                          JATP = IBLOCKR(1,IR)
                          JBTP = IBLOCKR(2,IR)
                          JASM = IBLOCKR(3,IR)
                          JBSM = IBLOCKR(4,IR)
                          JOFF = IBLOCKR(5,IR)
*
                          NJA = NSSOA(JASM,JATP)
                          NJB = NSSOB(JBSM,JBTP)
*
                          IAEXC = ICONSPA(JATP,IIATP)
                          IBEXC = ICONSPB(JBTP,IIBTP)
*
                          IF(IAEXC.EQ.0.AND.JASM.NE.IIASM) IAEXC = 1
                          IF(IBEXC.EQ.0.AND.JBSM.NE.IIBSM) IBEXC = 1
                          IABEXC = IAEXC + IBEXC
*
                          IF(IABEXC.LE.I12) THEN
                             INTERACT = 1
                          ELSE
                             INTERACT = 0
                          END IF
*
                          IF(INTERACT.EQ.1) THEN
*. Possible permutations of this block
                             CALL PRMBLK(  IDC,ISTRFL,JASM,JBSM,JATP,
     &                                    JBTP,  PSR,  PLR, RATP, RBTP,
     &                                    RASM, RBSM, RSGN, RTRP,NRPERM)
*. Well, spin permutations are simple to handle
* if there are two terms just calculate and and multiply with
* 1+PSL*PSR
                             IF(NRPERM.EQ.1) THEN
                                FACTOR = 1.0D0
                             ELSE
                                FACTOR = 1.0D0 +PSL*PSR
                             END IF
                             SCLFAC = FACTOR*SCLFAC_L(IL)*SCLFAC_R(IR)
                             IF(INTERACT.EQ.1.AND.SCLFAC.NE.0.0D0) THEN
                                IF(NTEST.GE.20) THEN
                                WRITE(6,*) ' RSDNBB will be called for '
                                WRITE(6,*) ' L block : '
                                WRITE(6,'(A,5I5)')
     &                               ' IIASM IIBSM IIATP IIBTP',
     &                               IIASM,IIBSM,IIATP,IIBTP
                                WRITE(6,*) ' R  block : '
                                WRITE(6,'(A,5I5)')
     &                               ' JASM JBSM JATP JBTP',
     &                               JASM,JBSM,JATP,JBTP
                                WRITE(6,*) ' IOFF,JOFF ', IOFF,JOFF
                                WRITE(6,*) ' SCLFAC = ', SCLFAC
                                END IF



                                CALL GSDNBB2_LUCIA(I12,
     &                                             RHO1,
     &                                             RHO2,
     &                                             RHO2S,
     &                                             RHO2A,
*
     &                                             IIASM,
     &                                             IIATP,
     &                                             IIBSM,
     &                                             IIBTP,
     &                                             JASM,
*
     &                                             JATP,
     &                                             JBSM,
     &                                             JBTP,
     &                                             NGAS,
     &                                       NELFSPGP(1,IOCTPA-1+IIATP),
*
     &                                       NELFSPGP(1,IOCTPB-1+IIBTP),
     &                                        NELFSPGP(1,IOCTPA-1+JATP),
     &                                        NELFSPGP(1,IOCTPB-1+JBTP),
     &                                             NAEL,
     &                                             NBEL,
*
     &                                             IAGRP,
     &                                             IBGRP,
     &                                             SB(IOFF),
     &                                             CB(JOFF),
     &                                              C2,
*
     &                                             ADSXA,
     &                                             SXSTST,
     &                                             STSTSX,
     &                                             DXSTST,
     &                                             STSTDX,
*
     &                                             SXDXSX,
     &                                             MXPNGAS,
     &                                             NOBPTS,
     &                                             IOBPTS,
     &                                             MAXI,
*
     &                                             MAXK,
     &                                             SSCR,CSCR,I1,XI1S,I2,
     &                                             XI2S,I3,XI3S,I4,XI4S,
     &                                               X,
     &                                             NSMOB,
*
     &                                             NSMST,
     &                                             NSMSX,
     &                                             NSMDX,
     &                                             NIIA,
     &                                             NIIB,
*
     &                                             NJA,
     &                                             NJB,
     &                                             MXPOBS,
     &                                             IPRNT,
     &                                             NACOB,
*
     &                                             RHO1S,
     &                                             SCLFAC,
     &                                             S2_TERM1,
     &                                             IUSE_PH,
     &                                             IPHGAS,
*
     &                                             IDOSRHO1,
     &                                             SRHO1,
     &                                             IPACK)
*
* CALL GSDNBB2_LUCIA --> 66
*
                                IF(NTEST.GE.500) THEN
                                   write(6,*) ' Updated rho1 '
                                   call wrtmat(rho1,nacob,nacob,nacob,
     &                                  nacob)
                                   write(6,*) ' Updated srho1 '
                                   call wrtmat(srho1,nacob,nacob,nacob,
     &                                  nacob)
                                END IF
*
                             END IF
                          END IF
                       END IF
 9000               CONTINUE
*. End of loop over R blocks in Batch
 9999            CONTINUE
*. Transpose or scale L block to restore order ??
                 IF(LTRP(NLPERM+1).EQ.1) THEN
                    CALL TRPMT3(SB(IOFF),NIB,NIA,C2)
                    CALL COPVEC(C2,SB(IOFF),NIA*NIB)
                 END IF
                 IF(LSGN(NLPERM+1).EQ.-1)
     &                CALL SCALVE(SB(IOFF),-1.0D0,NIA*NIB)
*
              END IF
10000      CONTINUE
*. End of loop over L blocks in batch
 9001   CONTINUE
*.      ^ End of loop over batches of R blocks
10001 CONTINUE
*.    ^ End of loop over batches of L blocks

      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer_array(ISSOA)
        CALL Unused_integer_array(ISSOB)
        CALL Unused_integer(LC)
        CALL Unused_integer(LS)
        CALL Unused_integer_array(SXSTSM)
        CALL Unused_integer_array(ASXAD)
        CALL Unused_real_array(RHO1P)
        CALL Unused_real_array(XNATO)
        CALL Unused_integer_array(LEBATL)
        CALL Unused_integer_array(LEBATR)
      END IF
      END
