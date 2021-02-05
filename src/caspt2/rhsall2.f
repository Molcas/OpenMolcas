************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE RHSALL2(IVEC)
      USE CHOVEC_IO
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb
      IMPLICIT REAL*8 (A-H,O-Z)
* ----------------------------------------------------------------
* Code for processing all the cholesky vectors
* in construction of caspt2 right-hand-side array
* Also form the active two-electron integrals 'TUVX'.
* ================================================================
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "WrkSpc.fh"
*
      Integer Active, Inactive, Virtual
      Parameter (Inactive=1, Active=2, Virtual=3)
      Integer nSh(8,3)
#ifdef _DEBUGPRINT_
      INTEGER NUMERR
      SAVE NUMERR
      DATA NUMERR / 0 /
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
      Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
      Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)
*                                                                      *
************************************************************************
*                                                                      *

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'
      END IF

*
* TUVX RHSX        Na^4
* TJVX RHSA        Na^3 Ni
* TJVL RHSB        Na^2 Ni^2
* AJVX RHSD1       Na^2 Ni   Ns
* AJCL RHSH             Ni^2 Ns^2     N^4
* AUVX RHSC        Na^3      Ns
* AUCX RHSF        Na^2      Ns^2
* AUVL RHSD2       Na^2 Ni   Ns
* AUCL RHSG        Na   Ni   Ns^2     N^3
* AJVL RHSE        Na   Ni^2 Ns       N^3
*
*     Initialize RHS array as 'vector' nr IVEC=IRHS on LUSOLV:
*
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate and clear TUVX, two-electron integrals for active
*     orbital indices only: Simple storage, same as for GAMMA2.
*     TUVX is kept allocated until end of subroutine.
      NG1=NASHT**2
      NG2=NG1**2
      NTUVX=NG2
      CALL GETMEM('TUVX','ALLO','REAL',LTUVX,NTUVX)
      CALL DCOPY_(NTUVX,[0.0D0],0,WORK(LTUVX),1)
*                                                                      *
************************************************************************
*                                                                      *
      DO JSYM=1,NSYM
*
       IB1=NBTCHES(JSYM)+1
       IB2=NBTCHES(JSYM)+NBTCH(JSYM)
*
       MXBGRP=IB2-IB1+1
       IF (MXBGRP.LE.0) CYCLE
       CALL GETMEM('BGRP','ALLO','INTE',LBGRP,2*MXBGRP)
       IBGRP=1
       DO IB=IB1,IB2
        IWORK(LBGRP  +2*(IBGRP-1))=IB
        IWORK(LBGRP+1+2*(IBGRP-1))=IB
        IBGRP=IBGRP+1
       END DO
       NBGRP=MXBGRP

       CALL MEMORY_ESTIMATE(JSYM,IWORK(LBGRP),NBGRP,
     &                      NCHOBUF,MXPIQK,NADDBUF)
       IF (IPRGLB.GT.VERBOSE) THEN
         WRITE(6,*)
         WRITE(6,'(A,I12)') '  Number of Cholesky batches: ',IB2-IB1+1
         WRITE(6,'(A,I12)') '  Number of batch groups:     ',NBGRP
         WRITE(6,*)
       END IF
* buffers are kept allocated until the end of JSYM loop.
       CALL GetMem('PIQK','ALLO','REAL',LPIQK,MXPIQK)
       CALL GetMem('BUFF','ALLO','REAL',LBUFF,NADDBUF)
       CALL GetMem('IDXB','ALLO','INTE',LIDXB,NADDBUF)
       CALL GETMEM('BRABUF','ALLO','REAL',LBRA,NCHOBUF)
       CALL GETMEM('KETBUF','ALLO','REAL',LKET,NCHOBUF)
*
*      Loop over groups of batches of Cholesky vectors
*
*      IBSTEP=1
*
*      DO IBSTA=IB1,IB2,IBSTEP
       DO IBGRP=1,NBGRP
        IBSTA=IWORK(LBGRP  +2*(IBGRP-1))
        IBEND=IWORK(LBGRP+1+2*(IBGRP-1))

        NV=0
        DO IB=IBSTA,IBEND
           NV=NV+NVLOC_CHOBATCH(IB)
        END DO

        IF (IPRGLB.GT.VERBOSE) THEN
         WRITE(6,'(A,I12)') '  Cholesky vectors in this group = ', NV
         WRITE(6,*)
        END IF
*                                                                      *
************************************************************************
*                                                                      *
*      Read kets (Cholesky vectors) in the form L(VX), all symmetries:
*
       Call Get_Cholesky_Vectors(Active,Active,JSYM,
     &                           Work(LKET),nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
*      Assemble constributions to TUVX integrals
*      Reuse the ket vectors as L(TU) bra vectors
*
       LBRASM=LKET
       DO ISYI=1,NSYM
        NI=NASH(ISYI)
        iOffi=NAES(iSYI)
        IF(NI.EQ.0) GOTO 115
        ISYP=MUL(ISYI,JSYM)
        NP=NASH(ISYP)
        iOffp=NAES(iSYP)
        IF(NP.EQ.0) GOTO 115
        NPI=NP*NI
        NBRASM=NPI*NV
        LKETSM=LKET

        DO ISYK=1,NSYM
         NK=NASH(ISYK)
         iOffK=NAES(iSYK)
         IF(NK.EQ.0) GOTO 112
         ISYQ=MUL(ISYK,JSYM)
         NQ=NASH(ISYQ)
         iOffQ=NAES(iSYQ)
         IF(NQ.EQ.0) GOTO 112
         NQK=NQ*NK
         NKETSM=NQK*NV
*
         IF (NPI*NQK.GT.mxPIQK) THEN
           WRITE(6,*) 'NPIQK larger than mxPIQK in TUVX, bug?'
           Call AbEnd()
         END IF
*        CALL GETMEM('PIQK','ALLO','REAL',LPIQK,NPI*NQK)
         CALL DGEMM_('N','T',NPI,NQK,NV,1.0D0,WORK(LBRASM),NPI,
     &        WORK(LKETSM),NQK,0.0D0,WORK(LPIQK),NPI)
*
         Call ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,
     &                WORK(LTUVX),nTUVX,Work(LPIQK),NPI*NQK,
     &                NUMERR)
*        CALL GETMEM('PIQK','FREE','REAL',LPIQK,NPI*NQK)
*
         LKETSM=LKETSM+NKETSM
 112     CONTINUE
        END DO
        LBRASM=LBRASM+NBRASM
 115    CONTINUE
       END DO
*                                                                      *
************************************************************************
*                                                                      *
*      Read bra (Cholesky vectors) in the form L(TJ): All symmetries
*
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                           Work(LBRA),nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
*      Assemble contributions to TJVX
*      Loop over the bras and kets, form <A|0>
*
      Call Process_RHS_Block(Inactive,Active,Active,Active,
     &                       'A ',
     &                       Work(LBRA),nBra,Work(LKET),nKet,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
*      TJVL RHSB
*      TJVL: Use TJ buffer as if it was VL, form <B|0>
*
      Call Process_RHS_Block(Inactive,Active,Inactive,Active,
     &                       'B ',
     &                       Work(LBRA),nBra,Work(LBRA),nBra,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
* We still have L(VX) vectors in core, at WORK(LKETS).
*
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                           Work(LBRA),nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVX RHSD1
* Loop over the bra and ket vectors.
*
      Call Process_RHS_Block(Inactive,Virtual,Active,Active,
     &                       'D1',
     &                       Work(LBRA),nBra,Work(LKET),nKet,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* AJCL RHSH
* AJCL: Use AJ buffer still in core as if it was CL, form <H|0>
*
      Call Process_RHS_Block(Inactive,Virtual,Inactive,Virtual,
     &                       'H ',
     &                       Work(LBRA),nBra,Work(LBRA),nBra,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read Bra (Cholesky vectors)= L(AU)
*
       Call Get_Cholesky_Vectors(Active,Virtual,JSYM,
     &                           Work(LBRA),nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUVX RHSC
* AUVX: Loop over the bras and kets
*
      Call Process_RHS_Block(Active,Virtual,Active,Active,
     &                       'C ',
     &                       Work(LBRA),nBra,Work(LKET),nKet,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* AUCX RHSF
* AUCX: Use AU buffer still in core as if it was CX, form <F|0>
*
      Call Process_RHS_Block(Active,Virtual,Active,Virtual,
     &                       'F ',
     &                       Work(LBRA),nBra,Work(LBRA),nBra,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(VL), all symmetries:
*
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                           Work(LKET),nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUVL RHSD2
* Loop over bras and kets, form <D2|0>.
*
      Call Process_RHS_Block(Active,Virtual,Inactive,Active,
     &                       'D2',
     &                       Work(LBRA),nBra,Work(LKET),nKet,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(CL), all symmetries:
*
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                           Work(LKET),nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUCL RHSG
* Loop over bras and kets, form  <G|0>
*
      Call Process_RHS_Block(Active,Virtual,Inactive,Virtual,
     &                       'G ',
     &                       Work(LBRA),nBra,Work(LKET),nKet,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read bra vectors AJ
*
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                           Work(LBRA),nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets in the form L(VL)
*
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                           Work(LKET),nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVL RHSE
* AJVL: Loop over bras and kets. Form <E|0>
*
      Call Process_RHS_Block(Inactive,Virtual,Inactive,Active,
     &                       'E ',
     &                       Work(LBRA),nBra,Work(LKET),nKet,
     &                       Work(LPIQK),mxPIQK,
     &                       Work(LBUFF),iWork(LidxB),nAddBuf,
     &                       nSh,JSYM,
     &                       IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* End of loop over batches, IB
      END DO
*                                                                      *
************************************************************************
*                                                                      *
*SVC  Call GetMem('ADDRHS','Free','Real',ipAdd,nAdd)
      CALL GETMEM('BRABUF','FREE','REAL',LBRA,NCHOBUF)
      CALL GETMEM('KETBUF','FREE','REAL',LKET,NCHOBUF)
      CALL GetMem('PIQK','FREE','REAL',LPIQK,MXPIQK)
      CALL GetMem('BUFF','FREE','REAL',LBUFF,NADDBUF)
      CALL GetMem('IDXB','FREE','INTE',LIDXB,NADDBUF)
      CALL GETMEM('BGRP','FREE','INTE',LBGRP,2*MXBGRP)
*                                                                      *
************************************************************************
*                                                                      *
* End of loop over JSYM
      END DO
*                                                                      *
************************************************************************
*                                                                      *
* Synchronized add RHS partial arrays from all nodes into each node.

C-SVC: read the DRA's from disk and copy them all to LUSOLV to continue
C      in serial mode.  FIXME: this call has to be removed when we reach
C      full parallel capabilities
*     CALL SYNRHS(IVEC)
C-SVC: at this point, the RHS elements are on disk, both in LUSOLV and
C      as DRAs with the name RHS_XX_XX_XX with XX a number representing
C      the case, symmetry, and rhs vector respectively.

* The RHS elements of Cases A, C, D1  need a correction:
      CALL MODRHS(IVEC,WORK(LFIMO))

#ifdef _DEBUGPRINT_
* compute and print RHS fingerprints
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') 'Case','Sym','Fingerprint'
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') '====','===','==========='
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NAS*NIS.NE.0) THEN
            CALL RHS_ALLO (NAS,NIS,lg_W)
            CALL RHS_READ (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
            DNRM2 = RHS_DDOT(NAS,NIS,lg_W,lg_W)
            WRITE(6,'(1X,I4,1X,I3,1X,F18.11)') ICASE,ISYM,DNRM2
          END IF
        END DO
      END DO
#endif

* Synchronized add tuvx partial arrays from all nodes into each node.
      CALL CHO_GADGOP(WORK(LTUVX),NTUVX,'+')
* Put TUVX on disk for possible later use:
      CALL PT2_PUT(NTUVX,'TUVX',WORK(LTUVX))
      CALL GETMEM('TUVX','FREE','REAL',LTUVX,NTUVX)

*                                                                      *
************************************************************************
*                                                                      *

      RETURN
      END
      Subroutine Get_Cholesky_Vectors(ITK,ITQ,JSYM,
     &                                Array,nArray,
     &                                IBSTA,IBEND)
      USE CHOVEC_IO
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      Real*8  Array(*)

      ! ugly hack to convert separate k/q orbital types into a specific
      ! case
      ICASE=ITK*ITQ
      IF (ICASE.EQ.3) THEN
        ICASE=4
      ELSE
        ICASE=ICASE/2
      END IF

      LKETSM=1
      DO ISYK=1,NSYM
        NQK=NPQ_CHOTYPE(ICASE,ISYK,JSYM)
        IF(NQK.EQ.0) CYCLE
        DO IB=IBSTA,IBEND
          NV=NVLOC_CHOBATCH(IB)
          NKETSM=NQK*NV
          IDISK=IDLOC_CHOGROUP(ICASE,ISYK,JSYM,IB)
          CALL DDAFILE(LUDRA,2,Array(LKETSM),NKETSM,IDISK)
          LKETSM=LKETSM+NKETSM
        END DO
      END DO
      nArray=LKETSM-1
*
      Return
      End
      Subroutine Process_RHS_Block(ITI,ITP,ITK,ITQ,
     &                             Case,
     &                             Cho_Bra,nBra,Cho_Ket,nKet,
     &                             PIQK,mxPIQK,
     &                             BUFF,idxBuff,nBUFF,
     &                             nSh,JSYM,
     &                             IVEC,NV)
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      DIMENSION Cho_Bra(nBra), Cho_Ket(nKet)
      DIMENSION BUFF(nBuff),idxBuff(nBuff),PIQK(mxPIQK)
      Integer nSh(8,3)
      Character Case*2
*
*
      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'Processing RHS block '//Case
      END IF
      LBRASM=1
      DO ISYI=1,NSYM
         NI=NSH(ISYI,ITI)
         IF(NI.EQ.0) GOTO 125
         ISYP=MUL(ISYI,JSYM)
         NP=NSH(ISYP,ITP)
         IF(NP.EQ.0) GOTO 125
         NPI=NP*NI
         NBRASM=NPI*NV
*
         LKETSM=1
         DO ISYK=1,NSYM
            NK=NSH(ISYK,ITK)
            IF(NK.EQ.0) GOTO 122
            ISYQ=MUL(ISYK,JSYM)
            NQ=NSH(ISYQ,ITQ)
            IF(NQ.EQ.0) GOTO 122
            NQK=NQ*NK
            NKETSM=NQK*NV
*
C SVC: we need an NPI*NQK to store the 2-electron integrals, and 2
C buffers (values+indices) for sorting them.  Later, we can try to get
C rid of the buffer that stores the values and only use an index buffer
C and the two-electron integrals for the scatter operation.  For the
C buffer, any size can be taken, but assuming there is enough memory
C available, it's set to the size of the two-electron integrals unless
C larger than some predefined maximum buffer size.
            NPIQK=NPI*NQK
            IF (NPIQK.GT.MXPIQK) THEN
              IF (Case.eq.'H') THEN
                KPI=MXPIQK/NQK
                NPIQK=KPI*NQK
              ELSE IF (Case.eq.'G') THEN
                KQK=MXPIQK/NPI
                NPIQK=NPI*KQK
              ELSE
                WRITE(6,*) ' NPIQK > MXPIQK and case != G or H'
                WRITE(6,'(A,A2)')  ' CASE =   ', Case
                WRITE(6,'(A,I12)') ' NPIQK =  ', NPIQK
                WRITE(6,'(A,I12)') ' MXPIQK = ', MXPIQK
                WRITE(6,*) ' This should not happen, please report.'
                CALL AbEnd()
              END IF
            END IF

C-SVC: sanity check
            IF (NPIQK.LE.0) THEN
              WRITE(6,'(1X,A)') ' ADDRHS: zero-sized NPIQK'
              CALL AbEnd()
            END IF
*
            If (Case.eq.'A ') Then
               CALL ADDRHSA(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxBuff,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'B ') Then
               CALL ADDRHSB(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxBuff,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'D1') Then
               CALL ADDRHSD1(IVEC,JSYM,ISYI,ISYK,
     &                       NP,NI,NQ,NK,PIQK,
     &                       nBuff,Buff,idxBuff,
     &                       Cho_Bra(LBRASM),
     &                       Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'H ') Then
               CALL ADDRHSH(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,NPIQK,
     &                      nBuff,Buff,idxBuff,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'C ') Then
               CALL ADDRHSC(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxBuff,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'F ') Then
               CALL ADDRHSF(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxBuff,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'D2') Then
               CALL ADDRHSD2(IVEC,JSYM,ISYI,ISYK,
     &                       NP,NI,NQ,NK,PIQK,
     &                       nBuff,Buff,idxBuff,
     &                       Cho_Bra(LBRASM),
     &                       Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'G ') Then
               CALL ADDRHSG(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,NPIQK,
     &                      nBuff,Buff,idxBuff,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'E ') Then
               CALL ADDRHSE(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxBuff,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else
               Call Abend()
            End If
*
         LKETSM=LKETSM+NKETSM
 122     CONTINUE
        END DO
        LBRASM=LBRASM+NBRASM
 125    CONTINUE
      END DO
*
      Return
      End
      Subroutine ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,
     &                   TUVX,nTUVX,PIQK,nPIQK,
     &                   NUMERR)
      Implicit Real*8 (A-H,O-Z)
      Real*8 TUVX(nTUVX), PIQK(nPIQK)
*
* Add into correct positions in TUVX:
*
      DO iX=0,NK-1
         iX1=NASHT*(iX+iOffK)
         iX2=NQ   * iX
         DO iV=0,NQ-1
            iVX1=NASHT*(iV+iOffQ+iX1)
            iVX2=NI   *(iV      +iX2)
            DO iU=0,NI-1
               iUVX1=NASHT*(iU+iOffI+iVX1)
               iUVX2=NP   *(iU+      iVX2)
#ifdef __INTEL_COMPILER
*  This to avoid Intel over optimization
               Call DaXpY_(nP,1.0D0,PIQK(1+      iUVX2),1,
     &                             TUVX(1+iOffP+iUVX1),1)
#else
               DO iT=0,NP-1
                  iTUVX=iT+iOffP+iUVX1
                  iPIQK=iT      +iUVX2
#ifdef _DEBUGPRINT_
* Temporary test statements -- remove after debug!
                  IF(ITUVX.LT.0 .or. ITUVX.gt.NTUVX) THEN
                     ITUVX=NTUVX
                     NUMERR=NUMERR+1
                     IF (NUMERR.GT.100) THEN
                        WRITE(6,*)' THIS IS TOO MUCH -- STOP.'
                        CALL QUIT(_RC_INTERNAL_ERROR_)
                     END IF
                  END IF
* End of temporary test statements
#else
* Avoid unused argument warnings
      IF (.FALSE.) Call Unused_integer(NUMERR)
#endif
                  TUVX(1+iTUVX)=TUVX(1+iTUVX)+PIQK(1+iPIQK)
               END DO
#endif
            END DO
         END DO
      END DO
*
      Return
      End
      SUBROUTINE MEMORY_ESTIMATE(JSYM,LBGRP,NBGRP,
     &                           NCHOBUF,NPIQK,NADDBUF)
      USE CHOVEC_IO
      use output, only:iPrGlb,verbose
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      DIMENSION LBGRP(2,NBGRP)
      Integer Active, Inactive, Virtual
      Parameter (Inactive=1, Active=2, Virtual=3)
      Integer nSh(8,3)
      DIMENSION ITYPE(4,9)
      DATA ITYPE /
     &  Inactive,Active,Active,Active,
     &  Inactive,Active,Inactive,Active,
     &  Inactive,Virtual,Active,Active,
     &  Inactive,Virtual,Inactive,Virtual,
     &  Active,Virtual,Active,Active,
     &  Active,Virtual,Active,Virtual,
     &  Active,Virtual,Inactive,Active,
     &  Active,Virtual,Inactive,Virtual,
     &  Inactive,Virtual,Inactive,Active /

      Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
      Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
      Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)

CSVC: compute the maximum RHS size that will be allocated per node, this
C     will be later allocated when addrhs routines are called. If global
C     arrays are used, this may actually be done outside of GETMEM, but
C     let's just reserve it anyway.
      ISYM=JSYM
      MXRHS=2*NTU(ISYM)*NISUP(ISYM,5)
      DO ISYI=1,NSYM
        ISYM=ISYI
C case A
        MXRHS=Max(MXRHS,NTUV(ISYM)*NISH(ISYM))
C case GP,GM
        MXRHS=Max(MXRHS,NASH(ISYM)*NISUP(ISYM,10)
     &                 +NASH(ISYM)*NISUP(ISYM,11))

        ISYM=MUL(JSYM,ISYI)
C case C
        MXRHS=Max(MXRHS,NTUV(ISYM)*NSSH(ISYM))
C case EP,EM
        MXRHS=Max(MXRHS,NASH(ISYM)*NISUP(ISYM,6)
     &                 +NASH(ISYM)*NISUP(ISYM,7))
        DO ISYK=1,NSYM
          ISYM=MUL(ISYI,ISYK)
C case BP,BM
          MXRHS=Max(MXRHS,NTGEU(ISYM)*NIGEJ(ISYM))
          MXRHS=Max(MXRHS,NTGTU(ISYM)*NIGTJ(ISYM))
C case HP,HM
          MXRHS=Max(MXRHS,NAGEB(ISYM)*NIGEJ(ISYM))
          MXRHS=Max(MXRHS,NAGTB(ISYM)*NIGTJ(ISYM))
C case F
          MXRHS=Max(MXRHS,NTGEU(ISYM)*NAGEB(ISYM))
          MXRHS=Max(MXRHS,NTGTU(ISYM)*NAGTB(ISYM))

          ISYM=MUL(ISYI,MUL(JSYM,ISYK))
C case D1,D2
          MXRHS=Max(MXRHS,2*NTU(ISYM)*NISUP(ISYM,5))
        End Do
      End Do
      MXRHS = iPARDIV(MXRHS,0)

CSVC: determine maximum pair index size per symmetry and in total.
*     MXBFSZ=0
      MXNPITOT=0
      DO ISYI=1,NSYM
        ISYP=MUL(ISYI,JSYM)
        NI=MAX(NISH(ISYI),NASH(ISYI))
        NP=MAX(NASH(ISYP),NSSH(ISYP))
        MXNPITOT=MXNPITOT+NP*NI
*       MXBFSZ=MXBFSZ+NP*NI*NJSCT
      END DO

CSVC: determine maximum and minimum size to hold integral matrix.
C     cases here correspond to A, B, D1, H, C, F, D2, G, E.
C     with number icase =      1, 2,  3, 4, 5, 6,  7, 8, 9.
C     both case G and H can use blocking, currently this is done over
C     the number of inactive orbitals in one of the pair indices.
C     for case G, the minimum needed is NAU*NC, because of blocking NL.
C     for case H, the minimum needed is NA*NCL, because of blocking NJ.
C     to start, reserve space for TUVX integrals (NASHT**4)
      MAXPIQK=NASHT**4
      MINPIQK=NASHT**4
      DO ICASE=1,9
        DO ISYI=1,NSYM
          NI=NSH(ISYI,ITYPE(1,ICASE))
          ISYP=MUL(ISYI,JSYM)
          NP=NSH(ISYP,ITYPE(2,ICASE))
          NPI=NP*NI
          DO ISYK=1,NSYM
            NK=NSH(ISYK,ITYPE(3,ICASE))
            ISYQ=MUL(ISYK,JSYM)
            NQ=NSH(ISYQ,ITYPE(4,ICASE))
            NQK=NQ*NK
            MAXPIQK=MAX(MAXPIQK,NPI*NQK)
            IF (ICASE.EQ.4) THEN
              MINPIQK=MAX(MINPIQK,NP*NQK)
            ELSE IF (ICASE.EQ.8) THEN
              MINPIQK=MAX(MINPIQK,NPI*NQ)
            ELSE
              MINPIQK=MAX(MINPIQK,NPI*NQK)
            END IF
          END DO
        END DO
      END DO
      MAXBUFF=NINT(SQRT(DBLE(MAXPIQK)))
      MINBUFF=NINT(SQRT(DBLE(MINPIQK)))

CSVC: total number of cholesky vectors
      MXBATCH=0
      NVECTOT=0
      IB1=NBTCHES(JSYM)+1
      IB2=NBTCHES(JSYM)+NBTCH(JSYM)
      DO IB=IB1,IB2
        NV=NVLOC_CHOBATCH(IB)
        MXBATCH=MAX(MXBATCH,NV)
        NVECTOT=NVECTOT+NV
      END DO
      MAXCHOL=MXNPITOT*NVECTOT
      MINCHOL=MXNPITOT*MXBATCH

CSVC: can we fit this all in memory?
      CALL GetMem('MAXSIZE','MAX','Real',iDum,MXAVAIL)

      MINNICE=MXRHS+MAXPIQK+2*MAXBUFF+2*MAXCHOL
      MINGOOD=MXRHS+MINPIQK+2*MINBUFF+2*MAXCHOL
      MINSLOW=MXRHS+MINPIQK+2*MINBUFF+2*MINCHOL

      IF (IPRGLB.GT.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(A,I1)') '  Memory estimates in RHSALL, SYM ', JSYM
        WRITE(6,'(A,2X,I16)') '   allocatable:    ', MXAVAIL
        WRITE(6,'(A,2X,I16)') '   recommended:    ', MINNICE
        WRITE(6,'(A,2X,I16)') '   convenient:     ', MINGOOD
        WRITE(6,'(A,2X,I16)') '   minimum:        ', MINSLOW
        WRITE(6,*)
        IF (MXAVAIL.GE.MINNICE) THEN
          WRITE(6,*) ' I can use all cholesky vectors at once'
          WRITE(6,*) ' as well as the whole integral matrixi.'
        ELSE IF (MXAVAIL.GE.MINGOOD) THEN
          WRITE(6,*) ' I will group batches of cholesky vectors'
          WRITE(6,*) ' and then maximize use of the integral matrix.'
        ELSE IF (MXAVAIL.GE.MINSLOW) THEN
          WRITE(6,*) ' I will at least try to group batches.'
        ELSE
          WRITE(6,*) ' Do you see my problem?'
        END IF
      END IF

      IF (MXAVAIL.GE.MINNICE) THEN
C group all batches, take the maximum needed for integrals and try
C to max out the buffer size, check it is larger than the minimum.
        NCHOBUF=MAXCHOL
        NBGRP=1
        LBGRP(1,1)=IB1
        LBGRP(2,1)=IB2
        NADDBUF=MAXBUFF
        NPIQK=MAXPIQK
      ELSE IF (MXAVAIL.GE.MINGOOD) THEN
C group all batches, take smaller buffer size, and try to max out
C integrals, and check they are lager than minimum needed
        NCHOBUF=MAXCHOL
        NBGRP=1
        LBGRP(1,1)=IB1
        LBGRP(2,1)=IB2
        NADDBUF=MINBUFF
        NCHUNK=(MXAVAIL-MXRHS-2*NADDBUF-2*NCHOBUF)/MINPIQK
        NPIQK=MINPIQK*NCHUNK
      ELSE IF (MXAVAIL.GE.MINSLOW) THEN
        NADDBUF=MINBUFF
        NPIQK=MINPIQK
        NCHOBUF=(MXAVAIL-MXRHS-NPIQK-2*NADDBUF)/2
C create batch groups that have at most MXCHOVEC cholesky vectors
        MXCHOVEC=MAX(NCHOBUF/MXNPITOT,1)
        NCHOVEC=0
        IBGRP=1
        LBGRP(1,IBGRP)=IB1
        DO IB=IB1,IB2
          NV=NVLOC_CHOBATCH(IB)
          NCHOVEC=NCHOVEC+NV
          IF (NCHOVEC.GT.MXCHOVEC) THEN
            LBGRP(2,IBGRP)=IB-1
            IBGRP=IBGRP+1
            LBGRP(1,IBGRP)=IB
            NCHOVEC=NV
          END IF
        END DO
        LBGRP(2,IBGRP)=IB2
        NBGRP=IBGRP
      ELSE
        WRITE(6,*)
        WRITE(6,*) '  Not enough memory in RHSLL2...'
        WRITE(6,'(A,I16)') '   allocatable:    ', MXAVAIL
        WRITE(6,'(A,I16)') '   minimum:        ', MINSLOW
        Call AbEnd()
      END IF

CSVC: sanity check, should not happen.
      IF (MXRHS.GT.MXAVAIL-2*NCHOBUF-NPIQK-2*NADDBUF) THEN
        WRITE(6,*)
        WRITE(6,*) 'RHSALL2: RHS allocation will starve.'
        WRITE(6,*) 'Possible bug in memory estimate.'
        WRITE(6,*) 'This should not happen, please report.'
        WRITE(6,*)
        WRITE(6,'(2X,A8,2X,I14)') 'MXAVAIL ', MXAVAIL
        WRITE(6,'(2X,A8,2X,I14)') 'MXRHS   ', MXRHS
        WRITE(6,'(2X,A8,2X,I14)') 'NCHOBUF ', NCHOBUF
        WRITE(6,'(2X,A8,2X,I14)') 'NPIQK   ', NPIQK
        WRITE(6,'(2X,A8,2X,I14)') 'NADDBUF ', NADDBUF
        CALL AbEnd()
      END IF

      IF (IPRGLB.GT.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(A16,2A16)') '  Buffer sizes:',
     &                        '           used',
     &                        '          ideal'
        WRITE(6,'(A16,2I16)') '  ChoVecs:  ', 2*NCHOBUF, 2*MAXCHOL
        WRITE(6,'(A16,2I16)') '  Integral: ',   NPIQK  ,   MAXPIQK
        WRITE(6,'(A16,2I16)') '  Scatter:  ', 2*NADDBUF, 2*MAXBUFF
        WRITE(6,*)
      END IF

      END
