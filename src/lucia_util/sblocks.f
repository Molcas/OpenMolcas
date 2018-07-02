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
* Copyright (C) 1991,1999, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE SBLOCKS( NSBLOCK, ISBLOCK,      CB,      SB,      C2,
     &                     ICOCOC,  ICSMOS,  ICBLTP,      NSSOA,
     &                      NSSOB,    NAEL,   IAGRP,    NBEL,   IBGRP,
     &                     IOCTPA,  IOCTPB,  NOCTPA,  NOCTPB,   NSMST,
     &                      NSMOB,   NSMSX,   NSMDX,  NOBPTS,  IOBPTS,
*
     &                    MXPNGAS,   ITSOB,      MAXK,    MAXI,
     &                                  LC,             XINT,    CSCR,
     &                       SSCR,    STSTSX,  STSTDX,  SXDXSX,
     &                      ADSXA,      NGAS,NELFSPGP,     IDC,
     &                         I1,    XI1S,      I2,    XI2S,   IDOH2,
*
     &                     MXPOBS,  ISTRFL,      PS,   IPRNT,     LUC,
     &                    ICJKAIB,   CJRES,   SIRES,      I3,    XI3S,
     &                         I4,    XI4S,  MXSXST,  MXSXBL,   MOCAA,
     &                         LCBLOCK,LECBLOCK,I1CBLOCK,
     &                    ICBLOCK,IRESTRICT,ICONSPA,ICONSPB, SCLFAC,
*
     &                      IPERTOP,IH0INSPC,  IH0SPC,ICBAT_RES,
     &                   ICBAT_INI,ICBAT_END,IUSE_PH, IPHGAS,I_RES_AB,
     &                          ISIMSYM,   XINT2)
*
* SUBROUTINE SBLOCKS --> 91
*
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen , Winter of 1991
*               Last modification : April 99
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If IRESTRICT.NE. 0 THEN we are after :
* sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)


      IMPLICIT REAL*8(A-H,O-Z)
*. Specific input
      INTEGER ISBLOCK(8,*)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),XINT2(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),XI1S(*),XI2S(*),XI3S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
      DIMENSION ISTRFL(*)
*. Zero order Hamiltonian
      INTEGER IH0SPC(NOCTPA,NOCTPB)
      INTEGER IH0INSPC(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION SCLFAC(*)
#include "bk_approx.fh"
#include "io_util.fh"
*
      COMMON/H_OCC_CONS/IH_OCC_CONS

      DIMENSION C(1),ICOOSC(1),IPHGAS(*)
      INTEGER DXSTST(1)
* IH_OCC_CONS =1 implies that we should employ occupation conserving
* part of Hamiltonian


*.
C-jwk-cleanup      REAL * 8 INPROD
*
      CALL QENTER('SBLOC')
*
C?    WRITE(6,*) ' IPERTOP in SBLOCKS = ', IPERTOP
c      IF(IH_OCC_CONS.EQ.1) THEN
c         WRITE(6,*) ' Occupation conserving part of Hamiltonian '
c      END IF
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
         WRITE(6,*) ' ================='
         WRITE(6,*) ' SBLOCKS speaking :'
         WRITE(6,*) ' ================='
         WRITE(6,*)
         WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &        NSBLOCK
         WRITE(6,*) ' TTSS for each ACTIVE sigma block'
         DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0)
     &           WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
         END DO
         WRITE(6,*) ' IDC PS IPERTOP', IDC,PS,IPERTOP
         WRITE(6,*) ' IDOH2 = ',IDOH2
         WRITE(6,*) ' I_RES_AB=',I_RES_AB
         if ( DoBKAP ) then
          write(6,*) ' I am doing BK-type of approximation '
          write(6,*) ' It is based on orbital splitting '
          write(6,*) ' Min and Max for subspace with exact Hamiltonian '
          write(6,*) ' =============================================== '
          write(6,*) 'NGASBK : ',NGASBK
          write(6,*) '              Min. Occ.      Max. Occ.           '
          Do IGAS = 1, NGASBK
            write(6,'(A,I2,10X,I3,9X,I3)')
     &      '   GAS',IGAS,IOCCPSPC(IGAS,1),IOCCPSPC(IGAS,2)
          End do
         end if
      END IF
*
      IF(NTEST.GE.50) THEN
         WRITE(6,*) ' Initial C vector '
         CALL WRTVCD(CB,LUC,1,-1)
      END IF
* ===========================
* 1 : Arrays for accessing C
* ============================
*. Find batches of C - strings
      CALL PART_CIV2(      IDC,   ICBLTP,    NSSOA,    NSSOB,   NOCTPA,
     &                  NOCTPB,    NSMST,       LC,   ICOCOC,   ICSMOS,
     &                 NCBATCH,  LCBLOCK, LECBLOCK, I1CBLOCK,  ICBLOCK,
     &                       0,  ISIMSYM)
*. Find the active blocks on LUC, store info in SCLFAC
      CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)
*
* Initialize sigma blocks
      DO JSBLOCK = 1, NSBLOCK
         IATP = ISBLOCK(1,JSBLOCK)
         IBTP = ISBLOCK(2,JSBLOCK)
         IASM = ISBLOCK(3,JSBLOCK)
         IBSM = ISBLOCK(4,JSBLOCK)
         IOFF = ISBLOCK(5,JSBLOCK)
         NASTR = NSSOA(IASM,IATP)
         NBSTR = NSSOB(IBSM,IBTP)
         ZERO = 0.0D0
         IF(ISBLOCK(1,JSBLOCK).GT.0)
     &        CALL SETVEC(SB(IOFF),ZERO,NASTR*NBSTR)
      END DO
* Loop over batches over C blocks
      IF(IDOH2.EQ.1) THEN
         MXEXC  = 2
      ELSE
         MXEXC = 1
      END IF
      IDISK(LUC)=0
      IF(ICBAT_RES.EQ.1) THEN
         WRITE(6,*) ' Restricted set of C batches '
         WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
         JCBAT_INI = ICBAT_INI
         JCBAT_END = ICBAT_END
      ELSE
         JCBAT_INI = 1
         JCBAT_END = NCBATCH
      END IF
*
      JOFF = 0 ! jwk-cleanup
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
*
*. Read C blocks into core
*
         ICOFF = 1
         NJBLOCK = LCBLOCK(JCBATCH)
         DO JJCBLOCK = 1, NJBLOCK
            JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
*. Will this block be needed ??
            INTERACT = 0
            IF(SCLFAC(JBLOCK).EQ. 1.0D0) THEN
               JATP = ICBLOCK(1,JBLOCK)
               JBTP = ICBLOCK(2,JBLOCK)
               JASM = ICBLOCK(3,JBLOCK)
               JBSM = ICBLOCK(4,JBLOCK)
               JOFF = ICBLOCK(5,JBLOCK)
               PL=1.D0
               CALL PRMBLK(     IDC,  ISTRFL,    JASM,    JBSM,    JATP,
     &                         JBTP,      PS,      PL,    LATP,    LBTP,
     &                         LASM,    LBSM,    LSGN,    LTRP,   NPERM)
               DO IPERM = 1, NPERM
                  LLASM = LASM(IPERM)
                  LLBSM = LBSM(IPERM)
                  LLATP = LATP(IPERM)
                  LLBTP = LBTP(IPERM)
*.Loop over Sigma blocks in batch
                  DO JSBLOCK = 1, NSBLOCK
                     IDENT = 0
                     IF(ISBLOCK(1,JSBLOCK).GT.0) THEN
                        IATP = ISBLOCK(1,JSBLOCK)
                        IBTP = ISBLOCK(2,JSBLOCK)
                        IASM = ISBLOCK(3,JSBLOCK)
                        IBSM = ISBLOCK(4,JSBLOCK)
*. Are the two blocks connected by allowed excitation
            CALL CON_BLOCKS( IATP, IBTP,LLATP,LLBTP, IASM,
     &                                   IBSM,LLASM,LLBSM,
     &                                  ICONSPA,
     &                                  ICONSPB,
     &                                  NOCTPA,
*
     &                                  NOCTPB,
     &                                  MXEXC,
     &                                  IH_OCC_CONS,
     &                                  INTERACT)
                        IDENT = 0
                        IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &                       IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
                     END IF
                  END DO
               END DO
*.             ^ End of checking whether C-block is needed
            END IF
*           ^ Checking was only done for nonvanishing blocks
*
            ISCALE = 0
            IF(INTERACT.EQ.1) THEN
               CALL GSTTBL(       C,CB(JOFF),    JATP,    JASM,    JBTP,
     &                         JBSM,  ICOCOC,  NOCTPA,  NOCTPB,   NSSOA,
     &                        NSSOB,      PS,  ICOOSC,     IDC,      PL,
     &                          LUC,C2,   NSMST,  ISCALE,SCLFAC(JBLOCK))
*. Note in GSTTBL : ICOOSC only used for CI vectors in core,
            ELSE
*. not relevant
               CALL IDAFILE(LUC,2,LBL,1,IDISK(LUC))
               CALL IDAFILE(LUC,2,iDUMMY,1,IDISK(LUC))
               CALL SKPRCD2(LBL,-1,LUC)
               SCLFAC(JBLOCK) = 0.0D0
            END IF
*
            IF(NTEST.GE.100) THEN
               IF(INTERACT.EQ.1) THEN
                  WRITE(6,*) ' TTSS for C block read in  '
                  CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
               ELSE
                  WRITE(6,*) ' TTSS for C block skipped  '
                  CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
               END IF
            END IF
*
         END DO
*        ^ End of loop over Blocks
*
*. Loop over blocks of sigma and C in core and obtain  contribution from
*  given C block to given S block
*. Loop over C blocks
         DO 9000 ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK
            JATP = ICBLOCK(1,ICBLK)
            JBTP = ICBLOCK(2,ICBLK)
            JASM = ICBLOCK(3,ICBLK)
            JBSM = ICBLOCK(4,ICBLK)
            ICOFF = ICBLOCK(5,ICBLK)
            NJA = NSSOA(JASM,JATP)
            NJB = NSSOB(JBSM,JBTP)
*
            IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
               CALL PRMBLK(     IDC,  ISTRFL,    JASM,    JBSM,    JATP,
     &                         JBTP,      PS,      PL,    LATP,    LBTP,
     &                         LASM,    LBSM,    LSGN,    LTRP,   NPERM)
*. Start with transposed block
               DO 8765 IPERM = NPERM,1, -1
                  LLASM = LASM(IPERM)
                  LLBSM = LBSM(IPERM)
                  LLATP = LATP(IPERM)
                  LLBTP = LBTP(IPERM)
                  NLLA = NSSOA(LLASM,LLATP)
                  NLLB = NSSOB(LLBSM,LLBTP)
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
                  IF(IPERM.EQ.1) THEN
                     IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM) THEN
*. Diagonal blocks, Transposing corresponds to scaling
                        IF(PS.EQ.-1.0D0) THEN
                           CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                        END IF
                     ELSE
*. ofdiagonal blocks, explicit transposing
                        CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
                        CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
                     END IF
                  END IF
*
                  DO 10000 ISBLK = 1, NSBLOCK
                     IF(ISBLOCK(1,ISBLK) .GT. 0 ) THEN
                        IATP = ISBLOCK(1,ISBLK)
                        IBTP = ISBLOCK(2,ISBLK)
                        IASM = ISBLOCK(3,ISBLK)
                        IBSM = ISBLOCK(4,ISBLK)
                        ISOFF = ISBLOCK(5,ISBLK)
                        NIA = NSSOA(IASM,IATP)
                        NIB = NSSOB(IBSM,IBTP)
*
                        IF(NIA*NIB.EQ.0) GOTO 10000
                        IF(IRESTRICT.EQ.1.AND.
     &                       (JASM.GT.IASM.OR.
     &                       JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &                  JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &                       GOTO 10000
*. Are the two blocks connected by allowed excitation
            CALL CON_BLOCKS( IATP, IBTP,LLATP,LLBTP, IASM,IBSM,
     &           LLASM,LLBSM,ICONSPA,ICONSPB,
     &           NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)

*. IF BK approximation is active, check whether block should
* be calculated exactly (1), by diagonal (-1) or is set to zero (0).
             I_DO_EXACT_BLK = 1
             IF(DoBKAP) THEN
               CALL CHECK_BLOCKS_FOR_BK_APPROX(
     &              IATP,IBTP,LLATP,LLBTP,
     &              IASM,IBSM,LLASM,LLBSM,
     &              IOCTPA,IOCTPB,I_DO_EXACT_BLK)
             END IF
C. BK-like approximation stuff
             IF(INTERACT.EQ.0.OR.I_DO_EXACT_BLK.EQ.0) GOTO 10000

             IF(NTEST.GE.100) THEN
               WRITE(6,*) ' Next s block in batch : '
               write(6,*) ' ISBLK IASM IBSM IATP IBTP'
               write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
             END IF
*
             IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &         ((LLBSM.GT.LLASM).OR.
     &         (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &         GOTO 8764
*
             IF(NTEST.GE.60) THEN
               WRITE(6,*) ' RSSBCB will be called for '
               WRITE(6,*) ' Sigma block : '
               WRITE(6,*) ' ISOFF ', ISOFF
               WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
               WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
               WRITE(6,*) ' C     block : '
               WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
               WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
               WRITE(6,*) ' ICOFF ', ICOFF
               WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
             END IF
*
             IF(IRESTRICT.EQ.1.AND.
     &            ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &            IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &            (IDC.EQ.2.AND.
     &            IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &         IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                XFAC = 0.5D0*SCLFAC(ICBLK)
             ELSE
                XFAC = SCLFAC(ICBLK)
             END IF
*. Form of operator in action
C               IF(IPERTOP.NE.0) THEN
*. Not exact Hamiltonian in use
                      IPTSPC = IH0SPC(IATP,IBTP)
                      JPTSPC = IH0SPC(JATP,JBTP)
                      IJOP   = IH0INSPC(IPTSPC)
*
                      IF(IPTSPC.NE.JPTSPC) GOTO 8764
*. BK-like approximation stuff
                IF(I_DO_EXACT_BLK.EQ.1) THEN
                  CALL RSSBCB2(IASM,  IATP,  IBSM,  IBTP,
     &                         LLASM, LLATP, LLBSM, LLBTP,
     &                         NGAS,
     &                         NELFSPGP(1,IATP+IOCTPA-1),
     &                         NELFSPGP(1,IBTP+IOCTPB-1),
     &                         NELFSPGP(1,LLATP+IOCTPA-1),
     &                         NELFSPGP(1,LLBTP+IOCTPB-1),
     &                         NAEL, NBEL, IAGRP, IBGRP,SB(ISOFF),
     &                         CB(ICOFF),IDOH2,ADSXA,STSTSX,
     &                         DXSTST,STSTDX,SXDXSX,NOBPTS,IOBPTS,
     &                         MXPNGAS,ITSOB, MAXI, MAXK, SSCR,
     &                         CSCR,    I1,  XI1S,    I2,  XI2S,
     &                         XINT,    C2, NSMOB, NSMST, NSMSX,
     &                         NSMDX,   NIA,   NIB,  NLLA,  NLLB,
     &                         MXPOBS,   IDC,    CJRES,
     &                         SIRES,    I3,  XI3S,    I4,  XI4S,
     &                         MXSXBL,MXSXST, MOCAA,
     &                         IPRNT,IPERTOP,
     &                         XFAC,IUSE_PH,IPHGAS,
     &                         I_RES_AB,
     &                         XINT2)
                ELSE IF(I_DO_EXACT_BLK.EQ.-1) THEN
*. Giovanni.... transposing sigma and CI vectors:
                  CALL TRPMT3(SB(ISOFF),NIB,NIA,C2)
                  CALL COPVEC(C2,SB(ISOFF),NIA*NIB)
                  CALL TRPMT3(CB(ICOFF),NLLB,NLLA,C2)
                  CALL COPVEC(C2,CB(ICOFF),NLLA*NLLB)
                  FACTOR = 0.0D0
                  CALL ADDDIA_TERM(FACTOR,CB(ICOFF),SB(ISOFF),
     &                             IATP,IBTP,IASM,IBSM)
*. Giovanni.... transposing back sigma and CI vectors:
                  CALL TRPMT3(SB(ISOFF),NIA,NIB,C2)
                  CALL COPVEC(C2,SB(ISOFF),NIA*NIB)
                  CALL TRPMT3(CB(ICOFF),NLLA,NLLB,C2)
                  CALL COPVEC(C2,CB(ICOFF),NLLA*NLLB)
                END IF ! End BK stuff
* CALL RSSBCB2 --> 82
 8764                 CONTINUE
                   END IF
*                  ^ End if S-block should be calculated
10000           CONTINUE
*.              ^  End of loop over sigma blocks
 8765        CONTINUE
          END IF
*         ^ End of C-block is nonvanishing
 9000  CONTINUE
*.     ^ End of loop over C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
*
* Order
      DO  ISBLK = 1 , NSBLOCK
         IF(ISBLOCK(1,ISBLK).GT.0) THEN
            IATP = ISBLOCK(1,ISBLK)
            IBTP = ISBLOCK(2,ISBLK)
            IASM = ISBLOCK(3,ISBLK)
            IBSM = ISBLOCK(4,ISBLK)
            ISOFF  = ISBLOCK(5,ISBLK)
            ISOFFP = ISBLOCK(6,ISBLK)
            NIA = NSSOA(IASM,IATP)
            NIB = NSSOB(IBSM,IBTP)
            IF(ICJKAIB.NE.0) THEN
*. Tranpose sigma block was obtained, transpose to obtain correct block
               CALL TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),
     &              NSSOA(IASM,IATP),C2)
               CALL COPVEC(C2,SB(ISOFF),
     &              NSSOA(IASM,IATP)* NSSOB(IBSM,IBTP))
            END IF
            IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP)
     &           CALL TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))
*
         END IF
      END DO
*
      IF(NTEST.GE.50) THEN
         WRITE(6,*) ' output blocks from SBLOCKS '
         CALL WRTTTS(       SB,  ISBLOCK,  NSBLOCK,    NSMST,
     &                   NSSOA,    NSSOB,        1)
      END IF
*
      CALL QEXIT('SBLOC')
      RETURN
      END
