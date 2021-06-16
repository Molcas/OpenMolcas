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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
C     based on rhsall2.f
C
C     In principle, For E = T_{ij}^{ab}*(ia|jb).
C     With p and q for general orbitals,
C     L_{pq} = (pa|jb)*T_{qj}^{ab} + (ip|jb)*T_{ij}^{qb}
C            + (ia|pb)*T_{iq}^{ab} + (ia|jp)*T_{ij}^{aq}
C     For the first term, with P and Q for auxiliary orbitals
C     L_{pq}(1) = (pa|jb)*T_{qj}^{ab}
C               = (pa|P)*(P|jb) * T_{qj}^{ab} (2c-2e is omitted?)
C               = (pa|P) * tilde{T}_{qa}^P
C               = C_{mu p}*C_{nu a}*(mu nu|P) * tilde{T}_{qa}^P
C               = C_{mu p} * (mu nu|P) * V_{nu q}^P
C     where
C     tilde{T}_{ia}^P = T_{ij}^{ab} * (P|jb)
C     V_{mu p}^P      = tilde{T}_{pa}^P * C_{mu a}
C
C     Dimension of tilde{T} will be the same as that of (ia|P) in MO,
C     where (i,a) = (inact,active), (inact,virtual), (active,virtual)
C
C     tilde{T} is constructed in this file.
C     tilde{T} -> V_{mu p}^P transformations, contraction with 3c-2e,
C     and construction of the orbital Lagrangian (L_{pq}) is elsewhere
C     ... maybe in OLagVVVO
C
C-----------------------------------------------------------------------
C
C     For the ERI derivative calculation,
C     d(mu nu|rho sigma)/da
C     = d/da (mu nu|P) (P|Q)^-1 (Q|rho sigma)
C     = d(mu nu|P)/da (P|Q)^-1 (Q|rho sigma)
C       + (mu nu|P) d(P|Q)^-1/da (Q|rho sigma)
C       + (mu nu|P) (P|Q)^-1 d(Q|rho sigma)/da
C     = d(mu nu|P)/da (P|Q)^-1 (Q|rho sigma)
C       - (mu nu|P) (P|R)^-1 d(R|S)/da (S|Q)^-1 (Q|rho sigma)
C       + (mu nu|P) (P|Q)^-1 d(Q|rho sigma)/da
C     = d(mu nu|P)/da (tP|rho sigma)
C       - (mu nu|tP) d(P|Q)/da (tQ|rho sigma)
C       + (mu nu|tP) d(P|rho sigma)/da
C     where (mu nu|tP) = (mu nu|Q)*(Q|P)^-1
C
C     D_{mu nu rho sigma}*d(mu nu|rho sigma)/da
C     = d(mu nu|P)/da (tP|rho sigma) * D_{mu nu rho sigma}
C       - (mu nu|tP) d(P|Q)/da (tQ|rho sigma) * D_{mu nu rho sigma}
C       + (mu nu|tP) d(P|rho sigma)/da * D_{mu nu rho sigma}
C     = d(mu nu|P)/da tD_{mu nu}^tP
C       - (mu nu|tP) d(P|Q)/da tD_{mu nu}^tQ
C       + tD_{rho sigma}^tP d(P|rho sigma)/da
C     where tD_{mu nu}^tP = D_{mu nu rho sigma} * (rho sigma|tP)
C     In practice, tD_{pq}^tP is constructed and saved in disk.
C     This will be read when 3c-2e ERI derivatives are concerned,
C     and MO->AO transformations will be done on-the-fly.
C     Note that MO coefficients of CASPT2 have to be used.
C
C     For 2c-2e ERI derivatives,
C     D(tP,tQ) = tD_{pq} * C_{mu p} C_{nu q} * (mu nu|tP)
C     then saved.
C
C     Subroutine OLagNS_RI(iSym0,WRK1,WRK2,DPT2C,BRAD,A_PT2,nChoVec)
      Subroutine OLagNS_RI(iSym0,WRK1,WRK2,DPT2C,
     *                     BraAI,BraSI,BraAA,BraSA,A_PT2,nChoVec)
C
      Use CHOVEC_IO
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Integer Active, Inactive, Virtual
      Parameter (Inactive=1, Active=2, Virtual=3)
      Integer nSh(8,3)
#ifdef _DEBUG_
      INTEGER NUMERR
      SAVE NUMERR
      DATA NUMERR / 0 /
#endif
C
      Dimension WRK1(*),WRK2(*),DPT2C(*),A_PT2(nChoVec,nChoVec)
      Dimension BraAI(nAsh(iSym0),nIsh(iSym0),nChoVec),
     *          BraSI(nSsh(iSym0),nIsh(iSym0),nChoVec),
     *          BraAA(nAsh(iSym0),nAsh(iSym0),nChoVec),
     *          BraSA(nSsh(iSym0),nAsh(iSym0),nChoVec)
C     Real*8, Target :: BraD(*)
C     Real*8, Pointer :: BraAI(:,:,:),BraSI(:,:,:),
C    *                   BraAA(:,:,:),BraSA(:,:,:)
C
      Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
      Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
      Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)
C
      Call QEnter('OLagNS_RI')

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'
      END IF
C
      iVec = iVecX
      iSym = iSym0
      !! Set pointers
      !! active-inactive
C     n = nAsh(iSym)*nIsh(iSym)*nChoVec
C     i = 1
C     j = n + i-1
C     BraAI(1:nAsh(iSym),1:nIsh(iSym),1:nChoVec) => BraD(i:j)
C     !! secondary-inactive
C     n = nSsh(iSym)*nIsh(iSym)*nChoVec
C     i = j+1
C     j = n + i-1
C     BraSI(1:nSsh(iSym),1:nIsh(iSym),1:nChoVec) => BraD(i:j)
C     !! active-active
C     n = nAsh(iSym)*nAsh(iSym)*nChoVec
C     i = j+1
C     j = n + i-1
C     BraAA(1:nAsh(iSym),1:nAsh(iSym),1:nChoVec) => BraD(i:j)
C     !! secondary-active
C     n = nSsh(iSym)*nAsh(iSym)*nChoVec
C     i = j+1
C     j = n + i-1
C     BraSA(1:nSsh(iSym),1:nAsh(iSym),1:nChoVec) => BraD(i:j)
*
      SCLNEL = 1.0D+00/DBLE(MAX(1,NACTEL))
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
     &                     NCHOBUF,MXPIQK,NADDBUF)
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
      CALL GETMEM('WRK','ALLO','REAL',ipWRK,nOrb(jSym)*nChoVec)
      CALL GETMEM('WRK2','ALLO','REAL',ipWRK2,
     *            max(norb(jsym)*nchovec,nAsh(jSym)**2))
C
C     Loop over groups of batches of Cholesky vectors
C
      DO IBGRP=1,NBGRP
C
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
*     Read kets (Cholesky vectors) in the form L(VX), all symmetries:
*
      Call Get_Cholesky_Vectors(Active,Active,JSYM,
     &                          Work(LKET),nKet,
     &                          IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
*     Assemble constributions to TUVX integrals
*     Reuse the ket vectors as L(TU) bra vectors
*
      LBRASM=LKET
      DO ISYI=1,NSYM
        NI=NASH(ISYI)
        iOffi=NAES(iSYI)
        IF(NI.EQ.0) Cycle
        ISYP=MUL(ISYI,JSYM)
        NP=NASH(ISYP)
        iOffp=NAES(iSYP)
        IF(NP.EQ.0) Cycle
        NPI=NP*NI
        NBRASM=NPI*NV
        LKETSM=LKET

        DO ISYK=1,NSYM
          NK=NASH(ISYK)
          iOffK=NAES(iSYK)
          IF(NK.EQ.0) Cycle
          ISYQ=MUL(ISYK,JSYM)
          NQ=NASH(ISYQ)
          iOffQ=NAES(iSYQ)
          IF(NQ.EQ.0) Cycle
          NQK=NQ*NK
          NKETSM=NQK*NV
*
          IF (NPI*NQK.GT.mxPIQK) THEN
            WRITE(6,*) 'NPIQK larger than mxPIQK in TUVX, bug?'
            Call AbEnd()
          END IF
*         CALL GETMEM('PIQK','ALLO','REAL',LPIQK,NPI*NQK)
          CALL DGEMM_('N','T',NPI,NQK,NV,1.0D0,WORK(LBRASM),NPI,
     &         WORK(LKETSM),NQK,0.0D0,WORK(LPIQK),NPI)
*
          Call ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,
     &                 WORK(LTUVX),nTUVX,Work(LPIQK),NPI*NQK,
     &                 NUMERR)
*         CALL GETMEM('PIQK','FREE','REAL',LPIQK,NPI*NQK)
*
          LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
      END DO
*                                                                      *
************************************************************************
*                                                                      *
*       Read bra (Cholesky vectors) in the form L(TJ): All symmetries
*
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
*      Assemble contributions to TJVX
*      Loop over the bras and kets, form <A|0>
*
      Call OLagNS_RI2(Inactive,Active,Active,Active,
     &                'A ',Work(LBRA),Work(LKET))
*                                                                      *
************************************************************************
*                                                                      *
*      TJVL RHSB
*      TJVL: Use TJ buffer as if it was VL, form <B|0>
*
      nKet = nBra
      Call OLagNS_RI2(Inactive,Active,Inactive,Active,
     &                'B ',Work(LBRA),Work(LBRA))
*                                                                      *
************************************************************************
*                                                                      *
* Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
* We still have L(VX) vectors in core, at WORK(LKETS).
*
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVX RHSD1
* Loop over the bra and ket vectors.
*
      Call OLagNS_RI2(Inactive,Virtual,Active,Active,
     &                'D1',Work(LBRA),Work(LKET))
*                                                                      *
************************************************************************
*                                                                      *
* AJCL RHSH
* AJCL: Use AJ buffer still in core as if it was CL, form <H|0>
*
      nKet = nBra
      Call OLagNS_RI2(Inactive,Virtual,Inactive,Virtual,
     &                'H ',Work(LBRA),Work(LBRA))
*                                                                      *
************************************************************************
*                                                                      *
* Read Bra (Cholesky vectors)= L(AU)
*
      Call Get_Cholesky_Vectors(Active,Virtual,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
C                                                                      *
************************************************************************
*                                                                      *
* AUVX RHSC
* AUVX: Loop over the bras and kets
*
      Call OLagNS_RI2(Active,Virtual,Active,Active,
     &                'C ',Work(LBRA),Work(LKET))
*                                                                      *
************************************************************************
*                                                                      *
* AUCX RHSF
* AUCX: Use AU buffer still in core as if it was CX, form <F|0>
*
      nKet = nBra
      Call OLagNS_RI2(Active,Virtual,Active,Virtual,
     &                'F ',Work(LBRA),Work(LBRA))
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(VL), all symmetries:
*
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                          Work(LKET),nKet,
     &                          IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUVL RHSD2
* Loop over bras and kets, form <D2|0>.
*
      Call OLagNS_RI2(Active,Virtual,Inactive,Active,
     &                'D2',Work(LBRA),Work(LKET))
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(CL), all symmetries:
*
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                          Work(LKET),nKet,
     &                          IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUCL RHSG
* Loop over bras and kets, form  <G|0>
*
      Call OLagNS_RI2(Active,Virtual,Inactive,Virtual,
     &                'G ',Work(LBRA),Work(LKET))
*                                                                      *
************************************************************************
*                                                                      *
* Read bra vectors AJ
*
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets in the form L(VL)
*
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                          Work(LKET),nKet,
     &                          IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVL RHSE
* AJVL: Loop over bras and kets. Form <E|0>
*
      Call OLagNS_RI2(Inactive,Virtual,Inactive,Active,
     &                'E ',Work(LBRA),Work(LKET))
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
      CALL GETMEM('WRK','FREE','REAL',ipWRK,nOrb(jSym)*nChoVec)
      CALL GETMEM('WRK2','FREE','REAL',ipWRK2,
     *            max(norb(jsym)*nchovec,nAsh(jSym)**2))
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
C     CALL MODRHS(IVEC,WORK(LFIMO))


* Synchronized add tuvx partial arrays from all nodes into each node.
C     CALL CHO_GADGOP(WORK(LTUVX),NTUVX,'+')
* Put TUVX on disk for possible later use:
C     CALL PT2_PUT(NTUVX,'TUVX',WORK(LTUVX))
      CALL GETMEM('TUVX','FREE','REAL',LTUVX,NTUVX)
C
      !! Symmetrize A_PT2...?
      Do iChoVec = 1, nChoVec
        Do jChoVec = 1, iChoVec-1
          Tmp = (A_PT2(iChoVec,jChoVec)+A_PT2(jChoVec,iChoVec))*0.5D+00
          A_PT2(iChoVec,jChoVec) = Tmp
          A_PT2(jChoVec,iChoVec) = Tmp
        End Do
      End Do
C
      Call DScal_(nBasSq,SCLNEL,DPT2C,1)
C
      Call QExit('OLagNS_RI')
C
      Return
C
      Contains
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_RI2(ITI,ITP,ITK,ITQ,Case,Cho_Bra,Cho_Ket)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
C     DIMENSION Cho_Bra(nBra), Cho_Ket(nKet)
      DIMENSION Cho_Bra(*), Cho_Ket(*)
      Character Case*2
C
      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'Processing RHS block '//Case
      END IF
C
      LBRASM=1
      CALL CWTime(TotCPU0,TotWall0)
      DO ISYI=1,NSYM
        NI=NSH(ISYI,ITI)
        IF(NI.EQ.0) CYCLE
        ISYP=MUL(ISYI,JSYM)
        NP=NSH(ISYP,ITP)
        IF(NP.EQ.0) CYCLE
        NPI=NP*NI
        NBRASM=NPI*NV
C
        LKETSM=1
        DO ISYK=1,NSYM
          NK=NSH(ISYK,ITK)
          IF(NK.EQ.0) CYCLE
          ISYQ=MUL(ISYK,JSYM)
          NQ=NSH(ISYQ,ITQ)
          IF(NQ.EQ.0) CYCLE
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
C
          IF (NPIQK.LE.0) THEN
            WRITE(6,'(1X,A)') ' ADDRHS: zero-sized NPIQK'
            CALL AbEnd()
          END IF
*
          !! NBUFF(=nAddBuf) is removed
          If (Case.eq.'A ') Then
             CALL OLagNS_RI_A(NP,NI,NQ,NK,
     &                        Cho_Bra(LBRASM),
     &                        Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'B ') Then
             CALL OLagNS_RI_B(NP,NI,NQ,NK,
     &                        Cho_Bra(LBRASM),
     &                        Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'D1') Then
             CALL OLagNS_RI_D1(NP,NI,NQ,NK,
     &                         Cho_Bra(LBRASM),
     &                         Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'H ') Then
             CALL OLagNS_RI_H(NP,NI,NQ,NK,NPIQK,
     &                        Cho_Bra(LBRASM),
     &                        Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'C ') Then
             CALL OLagNS_RI_C(NP,NI,NQ,NK,
     &                        Cho_Bra(LBRASM),
     &                        Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'F ') Then
             CALL OLagNS_RI_F(NP,NI,NQ,NK,
     &                        Cho_Bra(LBRASM),
     &                        Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'D2') Then
             CALL OLagNS_RI_D2(NP,NI,NQ,NK,
     &                         Cho_Bra(LBRASM),
     &                         Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'G ') Then
             CALL OLagNS_RI_G(NP,NI,NQ,NK,NPIQK,
     &                        Cho_Bra(LBRASM),
     &                        Cho_Ket(LKETSM),NV)
          Else If (Case.eq.'E ') Then
             CALL OLagNS_RI_E(NP,NI,NQ,NK,
     &                        Cho_Bra(LBRASM),
     &                        Cho_Ket(LKETSM),NV)
          Else
             Call Abend()
          End If
C
          LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
      END DO
      CALL CWTime(TotCPU1,TotWall1)
      write(6,'("CPU/Wall Time (Case ",A2,"):",2f10.2)')
     *  Case,totcpu1-totcpu0,totwall1-totwall0
C
      Return
C
      End Subroutine OLagNS_RI2
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_A(NT,NJ,NV,NX,
     &                       Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NX,NCHO)
C
      ISYJ = ISYI
      ISYX = ISYK
C
      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYX)
      ISYM=ISYJ
      IF(NINDEP(ISYM,1).EQ.0) RETURN
      NAS=NTUV(ISYM)
      NIS=NISH(ISYM)
      NWA=NAS*NIS
      IF(NWA.EQ.0) RETURN
C
C     ---- A
C
      !! Read the T-amplitude
      ICASE=1
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN.ne.0) Then
        nIS = nISup(iSym,iCase)
        nVec = nIN*nIS
        If (nVec.ne.0) Then
          Call RHS_ALLO(nAS,nIS,ipT)
          CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
        End If
      End If
C
      nOrbT = nFro(iSyT)+nIsh(iSyT)+nAsh(iSyT)+nSsh(iSyT)
      DO IT=1,NT
        ITABS=IT+NAES(ISYT)
        iTtot = iT + nFro(iSyT) + nIsh(iSyT)
        DO IJ=1,NJ
          IJABS=IJ+NIES(ISYJ)
          iJtot = iJ + nFro(iSyJ)
C
          Call DCopy_(NV*NX,[0.0D+00],0,WRK1,1)
          Call DCopy_(nAsh(iSym)*nAsh(iSym),[0.0D+00],0,WRK2,1)
          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
            IF (ISYV.EQ.ISYX) THEN !! not sure
              !! ONEADD contributions
              IW1=KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
              IW2=IJ
              IW=IW1+NAS*(IW2-1)
C
              ValAF = Work(ipT+IW-1)*2.0D+00
              DPT2C(iTtot+nOrbT*(iJtot-1))
     *          = DPT2C(iTtot+nOrbT*(iJtot-1)) + ValAF
            END IF
            DO IX=1,NX
              IXABS=IX+NAES(ISYX)
              IW1=KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
              IW2=IJ
              IW=IW1+NAS*(IW2-1)
C
              ValA = Work(ipT+IW-1)*2.0D+00
              WRK1(iX+NX*(iV-1)) = ValA
              WRK2(iXabs+nAsh(iSym)*(iVabs-1)) = ValA
C             Do iChoVec = 1, nCho
C               BraAA(iVabs,iXabs,iChoVec)
C    *            = BraAA(iVabs,iXabs,iChoVec)
C    *            + ValA*Cho_Bra(iT,iJ,iChoVec)
C               BraAI(iTabs,iJabs,iChoVec)
C    *            = BraAI(iTabs,iJabs,iChoVec)
C    *            + ValA*Cho_Ket(iV,iX,iChoVec)
C               Do jChoVec = 1, nCho
C                 A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *             + ValA*Cho_Bra(iT,iJ,iChoVec)
C    *                   *Cho_Ket(iV,iX,jChoVec)
C    *             + ValA*Cho_Bra(iT,iJ,jChoVec)
C    *                   *Cho_Ket(iV,iX,iChoVec)
C               End Do
C             End Do
            END DO
          END DO
          Call DGEMM_('N','N',nAsh(iSym)*nAsh(iSym),nCho,1,
     *                1.0D+00,WRK2,nAsh(iSym)*nAsh(iSym),
     *                        Cho_Bra(iT,iJ,1),NT*NJ,
     *                1.0D+00,BraAA(1,1,1),nAsh(iSym)*nAsh(iSym))
          Call DGEMV_('T',NV*NX,nCho,
     *                1.0D+00,Cho_Ket,NV*NX,WRK1,1,
     *                0.0D+00,WRK2,1)
          Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                BraAI(iTabs,iJabs,1),nAsh(iSym)*nIsh(iSym))
          Call DGEMM_('T','N',nCho,nCho,1,
     *                2.0D+00,Cho_Bra(iT,iJ,1),NT*NJ,WRK2,1,
     *                1.0D+00,A_PT2,nCho)
        END DO
      END DO
C
      CALL RHS_FREE(nAS,nIS,ipT)
C
      RETURN
C
      End Subroutine OLagNS_RI_A
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_B(NT,NJ,NV,NL,
     &                       Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NL,NCHO)
C
      ISYJ = ISYI
      ISYL = ISYK
C
      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYL)
      IF(ISYT.LT.ISYV) RETURN
      SQ2=SQRT(2.0D0)
      ISYM=MUL(ISYJ,ISYL) !!
C
      IF(NINDEP(ISYM,2).GT.0) THEN
* The plus combination:
       ICASE=2
       NASP=NTGEU(ISYM)
       NISP=NIGEJ(ISYM)
       NWBP=NASP*NISP
      ELSE
       NWBP=0
      ENDIF
      IF(NINDEP(ISYM,3).GT.0) THEN
* The minus combination:
       ICASE=3
       NASM=NTGTU(ISYM)
       NISM=NIGTJ(ISYM)
       NWBM=NASM*NISM
      ELSE
       NWBM=0
      ENDIF
      If (Max(NWBP,NWBM).le.0) RETURN
C
      IF(NWBP.GT.0.AND.NINDEP(ISYM,2).GT.0) THEN
        !! Read the T-amplitude
        ICASE=2
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP.ne.0) Then
          nISP = nISup(iSym,iCase)
          nVec = nINP*nISP
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASP,nISP,ipTP)
            CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
          End If
        End If
C
        DO IT=1,NT
          ITABS=IT+NAES(ISYT)
          IVMAX=NV
          IF(ISYV.EQ.ISYT) IVMAX=IT
          DO IV=1,IVMAX
            IVABS=IV+NAES(ISYV)
            SCL1=0.5D0
            IW1=KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
            IF(ITABS.EQ.IVABS) SCL1=0.25D0
C
            Call DCopy_(NJ*NL,[0.0D+00],0,WRK1,1)
            Call DCopy_(nIsh(iSym)*nIsh(iSym),[0.0D+00],0,WRK2,1)
            DO IJ=1,NJ
              IJABS=IJ+NIES(ISYJ)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                SCL=SCL1
                IF(IJABS.GE.ILABS) THEN
                 IW2=KIGEJ(IJABS,ILABS)-NIGEJES(ISYM)
                 IF(IJABS.EQ.ILABS) SCL=SQ2*SCL1
                ELSE
                 IW2=KIGEJ(ILABS,IJABS)-NIGEJES(ISYM)
                END IF
                IW=IW1+NASP*(IW2-1)
C
                ValBP = SCL*Work(ipTP+IW-1)*2.0d+00
                WRK1(iL+NL*(iJabs-1)) = ValBP
                WRK2(iLabs+nIsh(iSym)*(iJ-1)) = ValBP
C               Do iChoVec = 1, nCho
C                 BraAI(iVabs,iLabs,iChoVec)
C    *              = BraAI(iVabs,iLabs,iChoVec)
C    *              + ValBP*Cho_Bra(iT,iJ,iChoVec)
C                 BraAI(iTabs,iJabs,iChoVec)
C    *              = BraAI(iTabs,iJabs,iChoVec)
C    *              + ValBP*Cho_Ket(iV,iL,iChoVec)
C                 Do jChoVec = 1, nCho
C                   A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *               + ValBP*Cho_Bra(iT,iJ,iChoVec)
C    *                      *Cho_Ket(iV,iL,jChoVec)
C    *               + ValBP*Cho_Bra(iT,iJ,jChoVec)
C    *                      *Cho_Ket(iV,iL,iChoVec)
C                 End Do
C               End Do
              END DO
            END DO
            Call DCopy_(NL*NCHO,Cho_Ket(iV,1,1),NV,Work(ipWRK),1)
            Call DGEMM_('T','N',nIsh(iSym),NCHO,NL,
     *                  1.0D+00,WRK1,NL,Work(ipWRK),NL,
     *                  0.0D+00,Work(ipWRK2),nIsh(iSym))
            Do iChoVec = 1, nCho
              Call DaXpY_(nIsh(iSym),1.0D+00,
     *                    Work(ipWRK2+nIsh(iSym)*(iChoVec-1)),1,
     *                    BraAI(iTabs,1,iChoVec),nAsh(iSym))
            End Do
            !
            Call DCopy_(NJ*NCHO,Cho_Bra(iT,1,1),NT,Work(ipWRK),1)
            !! Check for NL.ne.nIsh etc
            Call DGEMM_('T','N',nCho,nCho,NL,
     *                  2.0D+00,Work(ipWRK2),NL,Work(ipWRK),NJ,
     *                  1.0D+00,A_PT2,nCho)
            Call DGEMM_('N','N',nIsh(iSym),NCHO,NJ,
     *                  1.0D+00,WRK2,NL,Work(ipWRK),NJ,
     *                  0.0D+00,Work(ipWRK2),nIsh(iSym))
            Do iChoVec = 1, nCho
              Call DaXpY_(nIsh(iSym),1.0D+00,
     *                    Work(ipWRK2+nIsh(iSym)*(iChoVec-1)),1,
     *                    BraAI(iVabs,1,iChoVec),nAsh(iSym))
            End Do
          END DO
        END DO
C
        CALL RHS_FREE(nASP,nISP,ipTP)
      END IF
C
      IF(NINDEP(ISYM,3).GT.0) THEN
        !! Read the T-amplitude
        ICASE=3
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM.ne.0) Then
          nISM = nISup(iSym,iCase)
          nVec = nINM*nISM
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASM,nISM,ipTM)
            CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
          End If
        End If
C
        DO IT=1,NT
          ITABS=IT+NAES(ISYT)
          IVMAX=NV
          IF(ISYV.EQ.ISYT) IVMAX=IT-1
          DO IV=1,IVMAX
            IVABS=IV+NAES(ISYV)
            IW1=KTGTU(ITABS,IVABS)-NTGTUES(ISYM)
C
            Call DCopy_(NJ*NL,[0.0D+00],0,WRK1,1)
            Call DCopy_(nIsh(iSym)*nIsh(iSym),[0.0D+00],0,WRK2,1)
            DO IJ=1,NJ
              IJABS=IJ+NIES(ISYJ)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                IF(IJABS.GT.ILABS) THEN
                  IW2=KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
                  SCL =  0.5D+00
                ELSE IF (IJABS.LT.ILABS) THEN
                  IW2=KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
                  SCL = -0.5D+00
                ELSE
                  CYCLE
                END IF
                IW=IW1+NASM*(IW2-1)
C
                ValBM = SCL*Work(ipTM+IW-1)*2.0D+00
                WRK1(iL+NL*(iJabs-1)) = ValBM
                WRK2(iLabs+nIsh(iSym)*(iJ-1)) = ValBM
C               Do iChoVec = 1, nCho
C                 BraAI(iVabs,iLabs,iChoVec)
C    *              = BraAI(iVabs,iLabs,iChoVec)
C    *              + ValBM*Cho_Bra(iT,iJ,iChoVec)
C                 BraAI(iTabs,iJabs,iChoVec)
C    *              = BraAI(iTabs,iJabs,iChoVec)
C    *              + ValBM*Cho_Ket(iV,iL,iChoVec)
C                 Do jChoVec = 1, nCho
C                   A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *               + ValBM*Cho_Bra(iT,iJ,iChoVec)
C    *                      *Cho_Ket(iV,iL,jChoVec)
C    *               + ValBM*Cho_Bra(iT,iJ,jChoVec)
C    *                      *Cho_Ket(iV,iL,iChoVec)
C                 End Do
C               End Do
              END DO
            END DO
            Call DCopy_(NL*NCHO,Cho_Ket(iV,1,1),NV,Work(ipWRK),1)
            Call DGEMM_('T','N',nIsh(iSym),NCHO,NL,
     *                  1.0D+00,WRK1,NL,Work(ipWRK),NL,
     *                  0.0D+00,Work(ipWRK2),nIsh(iSym))
            Do iChoVec = 1, nCho
              Call DaXpY_(nIsh(iSym),1.0D+00,
     *                    Work(ipWRK2+nIsh(iSym)*(iChoVec-1)),1,
     *                    BraAI(iTabs,1,iChoVec),nAsh(iSym))
            End Do
            !
            Call DCopy_(NJ*NCHO,Cho_Bra(iT,1,1),NT,Work(ipWRK),1)
            !! Check for NL.ne.nIsh etc
            Call DGEMM_('T','N',nCho,nCho,NL,
     *                  2.0D+00,Work(ipWRK2),NL,Work(ipWRK),NJ,
     *                  1.0D+00,A_PT2,nCho)
            Call DGEMM_('N','N',nIsh(iSym),NCHO,NJ,
     *                  1.0D+00,WRK2,NL,Work(ipWRK),NJ,
     *                  0.0D+00,Work(ipWRK2),nIsh(iSym))
            Do iChoVec = 1, nCho
              Call DaXpY_(nIsh(iSym),1.0D+00,
     *                    Work(ipWRK2+nIsh(iSym)*(iChoVec-1)),1,
     *                    BraAI(iVabs,1,iChoVec),nAsh(iSym))
            End Do
          END DO
        END DO
C
        CALL RHS_FREE(nASM,nISM,ipTM)
      END IF
C
      RETURN
C
      End Subroutine OLagNS_RI_B
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_C(NA,NU,NV,NX,
     &                       Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NX,NCHO)
C
      ISYU = ISYI
      ISYX = ISYK
C
      ISYA=MUL(JSYM,ISYU)
      ISYV=MUL(JSYM,ISYX)
      ISYM=ISYA !!
      IF(NINDEP(ISYM,4).EQ.0) RETURN
      NAS=NTUV(ISYM)
      NIS=NSSH(ISYM)
      NWC=NAS*NIS
      IF(NWC.EQ.0) RETURN
C
C     ---- C
C
      !! Read the T-amplitude
      ICASE=4
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN.ne.0) Then
        nIS = nISup(iSym,iCase)
        nVec = nIN*nIS
        If (nVec.ne.0) Then
          Call RHS_ALLO(nAS,nIS,ipT)
          CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
        End If
      End If
C
      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      DO IA=1,NA
        IAABS=IA+NSES(ISYA)
        iAtot = iA + nFro(iSyA) + nIsh(iSyA) + nAsh(iSyA)
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          iUtot = iU + nFro(iSyU) + nIsh(iSyU)
C
          Call DCopy_(NV*NX,[0.0D+00],0,WRK1,1)
          Call DCopy_(nAsh(iSym)*nAsh(iSym),[0.0D+00],0,WRK2,1)
          Call DCopy_(NU*NX,[0.0D+00],0,Work(ipWRK),1)
          Call DCopy_(nAsh(iSym)*nAsh(iSym),[0.0D+00],0,Work(ipWRK2),1)
          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
            iVtot = iV + nFro(iSyV) + nIsh(iSyV)
            ValCF = 0.0D+00
            IF (ISYV.EQ.ISYX) THEN !! not sure
              !! ONEADD contributions
              IW1=KTUV(IUABS,IVABS,IVABS)-NTUVES(ISYM)
              IW2=IA
              IW=IW1+NAS*(IW2-1)
C
              ValCF = Work(ipT+IW-1)*2.0D+00
              DPT2C(iAtot+nOrbA*(iUtot-1))
     *          = DPT2C(iAtot+nOrbA*(iUtot-1)) + ValCF
              ValCF = ValCF*SCLNEL
            END IF
            DO IX=1,NX
              IXABS=IX+NAES(ISYX)
              iXtot = iX + nFro(iSyX) + nIsh(iSyX)
              IW1=KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
              IW2=IA
              IW=IW1+NAS*(IW2-1)
C
              ValC = Work(ipT+IW-1)*2.0D+00
              WRK1(iV+NV*(iX-1)) = ValC
              WRK2(iVabs+nAsh(iSym)*(iXabs-1)) = ValC
              Work(ipWRK +iU-1+NU*(iXabs-1))
     *          = Work(ipWRK +iU-1+NU*(iXabs-1)) - ValCF
              Work(ipWRK2+iU-1+NU*(iX-1))
     *          = Work(ipWRK2+iU-1+NU*(iX-1)) - ValCF
C             Do iChoVec = 1, nCho
C               BraAA(iVabs,iXabs,iChoVec)
C    *            = BraAA(iVabs,iXabs,iChoVec)
C    *            + ValC*Cho_Bra(iA,iU,iChoVec)
C               BraSA(iAabs,iUabs,iChoVec)
C    *            = BraSA(iAabs,iUabs,iChoVec)
C    *            + ValC*Cho_Ket(iV,iX,iChoVec)
C               BraAA(iUabs,iXabs,iChoVec)
C    *            = BraAA(iUabs,iXabs,iChoVec)
C    *            - ValCF*Cho_Bra(iA,iX,iChoVec)
C               BraSA(iAabs,iXabs,iChoVec)
C    *            = BraSA(iAabs,iXabs,iChoVec)
C    *            - ValCF*Cho_Ket(iU,iX,iChoVec)
C               Do jChoVec = 1, nCho
C                 A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *             + ValC*Cho_Bra(iA,iU,iChoVec)
C    *                   *Cho_Ket(iV,iX,jChoVec)
C    *             + ValC*Cho_Bra(iA,iU,jChoVec)
C    *                   *Cho_Ket(iV,iX,iChoVec)
C    *             - ValCF*Cho_Bra(iA,iX,iChoVec)
C    *                    *Cho_Ket(iU,iX,jChoVec)
C    *             - ValCF*Cho_Bra(iA,iX,jChoVec)
C    *                    *Cho_Ket(iU,iX,iChoVec)
C               End Do
C             End Do
            END DO
          END DO
          Call DGEMM_('N','N',nAsh(iSym)*nAsh(iSym),nCho,1,
     *                1.0D+00,WRK2,nAsh(iSym)*nAsh(iSym),
     *                        Cho_Bra(iA,iU,1),NA*NU,
     *                1.0D+00,BraAA(1,1,1),nAsh(iSym)*nAsh(iSym))
          Call DGEMV_('T',NV*NX,nCho,
     *                1.0D+00,Cho_Ket,NV*NX,WRK1,1,
     *                0.0D+00,WRK2,1)
          Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                BraSA(iAabs,iUabs,1),nSsh(iSym)*nAsh(iSym))
          Call DGEMM_('T','N',nCho,nCho,1,
     *                2.0D+00,Cho_Bra(iA,iU,1),NA*NU,WRK2,1,
     *                1.0D+00,A_PT2,nCho)
C
          Call DGEMM_('T','N',1,NX*nCho,NU,
     *                1.0D+00,Work(ipWRK),nAsh(iSym),
     *                        Cho_Ket(1,1,1),NU,
     *                1.0D+00,BraSA(iAabs,1,1),nSsh(iSym))
          Call DCopy_(NX*NCHO,[0.0D+00],0,Work(ipWRK),1)
          Do iX = 1, NX
            Call DaXpY_(NCHO,Work(ipWRK2+iU-1+NU*(iX-1)),
     *                  Cho_Bra(iA,iX,1),NA*NX,
     *                  Work(ipWRK+NCHO*(iX-1)),1)
            IXABS=IX+NAES(ISYX)
            Call DaXpY_(NCHO,1.0D+00,Work(ipWRK+NCHO*(iX-1)),1,
     *                  BraAA(iUabs,iXabs,1),nAsh(iSym)**2)
            Call DGEMM_('T','N',nCho,nCho,1,
     *                  2.0D+00,Cho_Ket(iU,iX,1),NU*NX,
     *                          Work(ipWRK+NCHO*(iX-1)),1,
     *                  1.0D+00,A_PT2,nCho)
          End Do
        END DO
      END DO
C
      CALL RHS_FREE(nAS,nIS,ipT)
C
      RETURN
C
      End Subroutine OLagNS_RI_C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_D1(NA,NJ,NV,NX,
     &                        Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION Cho_Bra(NA,NJ,NCHO), Cho_Ket(NV,NX,NCHO)
*      Logical Incore
      DIMENSION IOFFD(8,8)
C
      ISYJ = ISYI
      ISYX = ISYK
C
      DO ISW=1,NSYM
       IO=0
       DO ISA=1,NSYM
        IOFFD(ISA,ISW)=IO
        ISI=MUL(ISA,ISW)
        IO=IO+NSSH(ISA)*NISH(ISI)
       END DO
      END DO

      ISYA=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYX)
      ISYM=JSYM !!!
      IF(NINDEP(ISYM,5).EQ.0) RETURN
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD.EQ.0) RETURN
C
C     ---- D1
C
      !! Read the T-amplitude
      ICASE=5
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN.ne.0) Then
        nIS = nISup(iSym,iCase)
        nVec = nIN*nIS
        If (nVec.ne.0) Then
          Call RHS_ALLO(nAS,nIS,ipT)
          CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
        End If
      End If
C
      NBXSZA=NSECBX
      NBXSZJ=NINABX
C
      DO IASTA=1,NA,NBXSZA
        IAEND=MIN(IASTA-1+NBXSZA,NA)
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
C
      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        iJtot = iJ + nFro(iSyJ)
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
          iAtot = iA + nFro(iSyA) + nIsh(iSyA) + nAsh(iSyA)
C
          Call DCopy_(NV*NX,[0.0D+00],0,WRK1,1)
          Call DCopy_(nAsh(iSym)*nAsh(iSym),[0.0D+00],0,WRK2,1)
          DO IX=1,NX
            IXABS=IX+NAES(ISYX)
            iXtot = iX + nFro(iSyX) + nIsh(iSyX)
            DO IV=1,NV
              IVABS=IV+NAES(ISYV)
              iVtot = iV + nFro(iSyV) + nIsh(iSyV)
              IW1=KTU(IVABS,IXABS)-NTUES(ISYM)
              IW2=IOFFD(ISYA,ISYM)+IJ+NJ*(IA-1)
              IW=IW1+NAS*(IW2-1)
C
              ValD = Work(ipT+IW-1)*2.0D+00
              If (iVtot.eq.iXtot) Then
                DPT2C(iAtot+nOrbA*(iJtot-1))
     *            = DPT2C(iAtot+nOrbA*(iJtot-1)) + ValD
              End If
              WRK1(iV+NV*(iX-1)) = ValD
              WRK2(iVabs+nAsh(iSym)*(iXabs-1)) = ValD
C             Do iChoVec = 1, nCho
C               BraAA(iVabs,iXabs,iChoVec)
C    *            = BraAA(iVabs,iXabs,iChoVec)
C    *            + ValD*Cho_Bra(iA,iJ,iChoVec)
C               BraSI(iAabs,iJabs,iChoVec)
C    *            = BraSI(iAabs,iJabs,iChoVec)
C    *            + ValD*Cho_Ket(iV,iX,iChoVec)
C               Do jChoVec = 1, nCho
C                 A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *             + ValD*Cho_Bra(iA,iJ,iChoVec)
C    *                   *Cho_Ket(iV,iX,jChoVec)
C    *             + ValD*Cho_Bra(iA,iJ,jChoVec)
C    *                   *Cho_Ket(iV,iX,iChoVec)
C               End Do
C             End Do
            END DO
          END DO
          Call DGEMM_('N','N',nAsh(iSym)*nAsh(iSym),nCho,1,
     *                1.0D+00,WRK2,nAsh(iSym)*nAsh(iSym),
     *                        Cho_Bra(iA,iJ,1),NA*NJ,
     *                1.0D+00,BraAA(1,1,1),nAsh(iSym)*nAsh(iSym))
          Call DGEMV_('T',NV*NX,nCho,
     *                1.0D+00,Cho_Ket,NV*NX,WRK1,1,
     *                0.0D+00,WRK2,1)
          Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                BraSI(iAabs,iJabs,1),nSsh(iSym)*nIsh(iSym))
          Call DGEMM_('T','N',nCho,nCho,1,
     *                2.0D+00,Cho_Bra(iA,iJ,1),NA*NJ,WRK2,1,
     *                1.0D+00,A_PT2,nCho)
        END DO
      END DO
C
        ENDDO
      ENDDO
C
      CALL RHS_FREE(nAS,nIS,ipT)
C
      RETURN
C
      End Subroutine OLagNS_RI_D1
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_D2(NA,NU,NV,NL,
     &                        Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
      DIMENSION IOFFD(8,8)
C
      ISYU = ISYI
      ISYL = ISYK
C
      DO ISYW=1,NSYM
       IO=0
       DO ISYA=1,NSYM
        IOFFD(ISYA,ISYW)=IO
        ISYII=MUL(ISYA,ISYW)
        IO=IO+NSSH(ISYA)*NISH(ISYII)
       END DO
      END DO

      ISYA=MUL(JSYM,ISYU)
      ISYV=MUL(JSYM,ISYL)
      ISYM=MUL(ISYU,ISYV)
      IF(NINDEP(ISYM,5).EQ.0) RETURN
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD.EQ.0) RETURN
C
C     ---- D2
C
      !! Read the T-amplitude
      ICASE=5
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN.ne.0) Then
        nIS = nISup(iSym,iCase)
        nVec = nIN*nIS
        If (nVec.ne.0) Then
          Call RHS_ALLO(nAS,nIS,ipT)
          CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
        End If
      End If

      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      DO IA=1,NA
        IAABS=IA+NSES(ISYA)
        iAtot = iA + nFro(iSyA) + nIsh(iSyA) + nAsh(iSyA)
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          iUtot = iU + nFro(iSyU) + nIsh(iSyU)
C
          Call DCopy_(NV*NL,[0.0D+00],0,WRK1,1)
          Call DCopy_(nAsh(iSym)*nIsh(iSym),[0.0D+00],0,WRK2,1)
          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
            iVtot = iV + nFro(iSyV) + nIsh(iSyV)
            DO IL=1,NL
              ILABS=IL+NIES(ISYL)
              iLtot = iL + nFro(iSyL)
              IW1=NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
              IW2=IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
              IW=IW1+NAS*(IW2-1)
C
              ValD = Work(ipT+IW-1)*2.0D+00
              WRK1(iV+NV*(iL-1)) = ValD
              WRK2(iVabs+nAsh(iSym)*(iLabs-1)) = ValD
C             Do iChoVec = 1, nCho
C               BraAI(iVabs,iLabs,iChoVec)
C    *            = BraAI(iVabs,iLabs,iChoVec)
C    *            + ValD*Cho_Bra(iA,iU,iChoVec)
C               BraSA(iAabs,iUabs,iChoVec)
C    *            = BraSA(iAabs,iUabs,iChoVec)
C    *            + ValD*Cho_Ket(iV,iL,iChoVec)
C               Do jChoVec = 1, nCho
C                 A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *             + ValD*Cho_Bra(iA,iU,iChoVec)
C    *                   *Cho_Ket(iV,iL,jChoVec)
C    *             + ValD*Cho_Bra(iA,iU,jChoVec)
C    *                   *Cho_Ket(iV,iL,iChoVec)
C               End Do
C             End Do
            END DO
          END DO
          Call DGEMM_('N','N',nAsh(iSym)*nIsh(iSym),nCho,1,
     *                1.0D+00,WRK2,nAsh(iSym)*nIsh(iSym),
     *                        Cho_Bra(iA,iU,1),NA*NU,
     *                1.0D+00,BraAI(1,1,1),nAsh(iSym)*nIsh(iSym))
          Call DGEMV_('T',NV*NL,nCho,
     *                1.0D+00,Cho_Ket,NV*NL,WRK1,1,
     *                0.0D+00,WRK2,1)
          Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                BraSA(iAabs,iUabs,1),nSsh(iSym)*nAsh(iSym))
          Call DGEMM_('T','N',nCho,nCho,1,
     *                2.0D+00,Cho_Bra(iA,iU,1),NA*NU,WRK2,1,
     *                1.0D+00,A_PT2,nCho)
        END DO
      END DO
C
      CALL RHS_FREE(nAS,nIS,ipT)
C
      RETURN
C
      End Subroutine OLagNS_RI_D2
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_E(NA,NJ,NV,NL,
     &                       Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION Cho_Bra(NA,NJ,NCHO), Cho_Ket(NV,NL,NCHO)
*      Logical Incore
      DIMENSION IOFF1(8),IOFF2(8)
C
      ISYJ = ISYI
      ISYL = ISYK
C
      SQ32=SQRT(1.5D0)
      ISYA=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYL)
      ISYM=ISYV
      ISYJL=MUL(ISYJ,ISYL)

C Set up offset table:
      IO1=0
      IO2=0
      DO ISA=1,NSYM
        IOFF1(ISA)=IO1
        IOFF2(ISA)=IO2
        ISIJ=MUL(ISA,ISYM)
        IO1=IO1+NSSH(ISA)*NIGEJ(ISIJ)
        IO2=IO2+NSSH(ISA)*NIGTJ(ISIJ)
      END DO

      NAS=NASH(ISYM)
      NISP=NISUP(ISYM,6)
      NISM=NISUP(ISYM,7)
      NIS=NISP+NISM
      NWP=NAS*NISP
      NWM=NAS*NISM
      NW=NWP+NWM
      If (NW.eq.0) RETURN
C
C     ---- EP
C
      IF (NWP.GT.0) THEN
        !! Read the T-amplitude
        ICASE=6
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP.ne.0) Then
          nISP = nISup(iSym,iCase)
          nVec = nINP*nISP
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASP,nISP,ipTP)
            CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
          End If
        End If

        NBXSZA=NSECBX
        NBXSZJ=NINABX

        DO IASTA=1,NA,NBXSZA
          IAEND=MIN(IASTA-1+NBXSZA,NA)
          NASZ=IAEND-IASTA+1
          DO IJSTA=1,NJ,NBXSZJ
            IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
            NJSZ=IJEND-IJSTA+1
C
        DO IJ=IJSTA,IJEND
          IJABS=IJ+NIES(ISYJ)
          DO IA=IASTA,IAEND
            IAABS=IA+NSES(ISYA)
C
            Call DCopy_(NV*NL,[0.0D+00],0,WRK1,1)
            Call DCopy_(nAsh(iSym)*nIsh(iSym),[0.0D+00],0,WRK2,1)
            DO IV=1,NV
              IVABS=IV+NAES(ISYV)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                SCL=SQRT(0.5D0)
                IF(IJABS.GE.ILABS) THEN
                  JGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                  IF(IJABS.EQ.ILABS) SCL=1.0D0
                ELSE
                  JGEL=KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
                END IF
                IW1=IV
                IW2=IA+NA*(JGEL-1)+IOFF1(ISYA)
                IW=IW1+NAS*(IW2-1)
C
                ValEP = SCL*Work(ipTP+IW-1)*2.0d+00
                WRK1(iV+NV*(iL-1)) = ValEP
                WRK2(iVabs+nAsh(iSym)*(iLabs-1)) = ValEP
C               Do iChoVec = 1, nCho
C                 BraAI(iVabs,iLabs,iChoVec)
C    *              = BraAI(iVabs,iLabs,iChoVec)
C    *              + ValEP*Cho_Bra(iA,iJ,iChoVec)
C                 BraSI(iAabs,iJabs,iChoVec)
C    *              = BraSI(iAabs,iJabs,iChoVec)
C    *              + ValEP*Cho_Ket(iV,iL,iChoVec)
C                 Do jChoVec = 1, nCho
C                   A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *               + ValEP*Cho_Bra(iA,iJ,iChoVec)
C    *                      *Cho_Ket(iV,iL,jChoVec)
C    *               + ValEP*Cho_Bra(iA,iJ,jChoVec)
C    *                      *Cho_Ket(iV,iL,iChoVec)
C                 End Do
C               End Do
              END DO
            END DO
            Call DGEMM_('N','N',nAsh(iSym)*nIsh(iSym),nCho,1,
     *                  1.0D+00,WRK2,nAsh(iSym)*nIsh(iSym),
     *                          Cho_Bra(iA,iJ,1),NA*NJ,
     *                  1.0D+00,BraAI(1,1,1),nAsh(iSym)*nIsh(iSym))
            Call DGEMV_('T',NV*NL,nCho,
     *                  1.0D+00,Cho_Ket,NV*NL,WRK1,1,
     *                  0.0D+00,WRK2,1)
            Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                  BraSI(iAabs,iJabs,1),nSsh(iSym)*nIsh(iSym))
            Call DGEMM_('T','N',nCho,nCho,1,
     *                  2.0D+00,Cho_Bra(iA,iJ,1),NA*NJ,WRK2,1,
     *                  1.0D+00,A_PT2,nCho)
          END DO
        END DO
C
          ENDDO
        ENDDO
C
        CALL RHS_FREE(nASP,nISP,ipTP)
      END IF
C
C     ---- EM
C
      IF (NWM.GT.0) THEN
        !! Read the T-amplitude
        ICASE=7
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM.ne.0) Then
          nISM = nISup(iSym,iCase)
          nVec = nINM*nISM
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASM,nISM,ipTM)
            CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
          End If
        End If
C
        NBXSZA=NSECBX
        NBXSZJ=NINABX
C
        DO IASTA=1,NA,NBXSZA
          IAEND=MIN(IASTA-1+NBXSZA,NA)
          NASZ=IAEND-IASTA+1
          DO IJSTA=1,NJ,NBXSZJ
            IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
            NJSZ=IJEND-IJSTA+1
C
        DO IJ=IJSTA,IJEND
          IJABS=IJ+NIES(ISYJ)
          DO IA=IASTA,IAEND
            IAABS=IA+NSES(ISYA)
C
            Call DCopy_(NV*NL,[0.0D+00],0,WRK1,1)
            Call DCopy_(nAsh(iSym)*nIsh(iSym),[0.0D+00],0,WRK2,1)
            DO IV=1,NV
              IVABS=IV+NAES(ISYV)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                IF(IJABS.NE.ILABS) THEN
                  IF(IJABS.GT.ILABS) THEN
                    SCL=SQ32
                    JGTL=KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                  ELSE
                    SCL=-SQ32
                    JGTL=KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
                  END IF
                  IW1=IV
                  IW2=IA+NA*(JGTL-1)+IOFF2(ISYA)
                  IW=IW1+NAS*(IW2-1)
C
                 ValEM = SCL*Work(ipTM+IW-1)*2.0d+00
                 WRK1(iV+NV*(iL-1)) = ValEM
                 WRK2(iVabs+nAsh(iSym)*(iLabs-1)) = ValEM
C                Do iChoVec = 1, nCho
C                  BraAI(iVabs,iLabs,iChoVec)
C    *               = BraAI(iVabs,iLabs,iChoVec)
C    *               + ValEM*Cho_Bra(iA,iJ,iChoVec)
C                  BraSI(iAabs,iJabs,iChoVec)
C    *               = BraSI(iAabs,iJabs,iChoVec)
C    *               + ValEM*Cho_Ket(iV,iL,iChoVec)
C                  Do jChoVec = 1, nCho
C                    A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *                + ValEM*Cho_Bra(iA,iJ,iChoVec)
C    *                       *Cho_Ket(iV,iL,jChoVec)
C    *                + ValEM*Cho_Bra(iA,iJ,jChoVec)
C    *                       *Cho_Ket(iV,iL,iChoVec)
C                  End Do
C                End Do
                END IF
              END DO
            END DO
            Call DGEMM_('N','N',nAsh(iSym)*nIsh(iSym),nCho,1,
     *                  1.0D+00,WRK2,nAsh(iSym)*nIsh(iSym),
     *                          Cho_Bra(iA,iJ,1),NA*NJ,
     *                  1.0D+00,BraAI(1,1,1),nAsh(iSym)*nIsh(iSym))
            Call DGEMV_('T',NV*NL,nCho,
     *                  1.0D+00,Cho_Ket,NV*NL,WRK1,1,
     *                  0.0D+00,WRK2,1)
            Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                  BraSI(iAabs,iJabs,1),nSsh(iSym)*nIsh(iSym))
            Call DGEMM_('T','N',nCho,nCho,1,
     *                  2.0D+00,Cho_Bra(iA,iJ,1),NA*NJ,WRK2,1,
     *                  1.0D+00,A_PT2,nCho)
          END DO
        END DO
C
          ENDDO
        ENDDO
C
        CALL RHS_FREE(nASM,nISM,ipTM)
      END IF
C
      RETURN
C
      End Subroutine OLagNS_RI_E
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_F(NA,NU,NC,NX,
     &                       Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NX,NCHO)
C
      ISYU = ISYI
      ISYX = ISYK
C
      IF(ISYU.LT.ISYX) RETURN
C
      ISYA=MUL(JSYM,ISYU)
      ISYC=MUL(JSYM,ISYX)
      ISYM=MUL(ISYU,ISYX) !!
C
      IF(NINDEP(ISYM,8).GT.0) THEN
* The plus combination:
       NASP=NTGEU(ISYM)
       NISP=NAGEB(ISYM)
       NWFP=NASP*NISP
      ELSE
       NWFP=0
      ENDIF
      IF(NINDEP(ISYM,9).GT.0) THEN
       ICASE=9
* The minus combination:
       NASM=NTGTU(ISYM)
       NISM=NAGTB(ISYM)
       NWFM=NASM*NISM
      ELSE
       NWFM=0
      ENDIF
      If (NWFP+NWFM.le.0) RETURN
C
C     ---- FP
C
      IF (NWFP.GT.0.AND.NINDEP(ISYM,8).GT.0) THEN
        !! Read the T-amplitude
        ICASE=8
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP.ne.0) Then
          nISP = nISup(iSym,iCase)
          nVec = nINP*nISP
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASP,nISP,ipTP)
            CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
          End If
        End If
C
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          IXMAX=NX
          IF(ISYU.EQ.ISYX) IXMAX=IU
          DO IX=1,IXMAX
            IXABS=IX+NAES(ISYX)
            SCL1=0.5D0
            IF(IUABS.EQ.IXABS) SCL1=0.25D0
            IW1=KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
C
            Call DCopy_(NA*NC,[0.0D+00],0,WRK1,1)
            Call DCopy_(nSsh(iSym)*nSsh(iSym),[0.0D+00],0,WRK2,1)
            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              DO IC=1,NC
                ICABS=IC+NSES(ISYC)
                SCL=SCL1
                IF(IAABS.GE.ICABS) THEN
                  IW2=KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
                  IF(IAABS.EQ.ICABS) SCL=SQRT(2.0D0)*SCL1
                ELSE
                  IW2=KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
                END IF
                IW=IW1+NASP*(IW2-1)
C
                ValFP = SCL*Work(ipTP+IW-1)*2.0d+00
                WRK1(iC+NC*(iA-1)) = ValFP
                WRK2(iCabs+nSsh(iSym)*(iA-1)) = ValFP
C               Do iChoVec = 1, nCho
C                 BraSA(iCabs,iXabs,iChoVec)
C    *              = BraSA(iCabs,iXabs,iChoVec)
C    *              + ValFP*Cho_Bra(iA,iU,iChoVec)
C                 BraSA(iAabs,iUabs,iChoVec)
C    *              = BraSA(iAabs,iUabs,iChoVec)
C    *              + ValFP*Cho_Ket(iC,iX,iChoVec)
C                 Do jChoVec = 1, nCho
C                   A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *               + ValFP*Cho_Bra(iA,iU,iChoVec)
C    *                      *Cho_Ket(iC,iX,jChoVec)
C    *               + ValFP*Cho_Bra(iA,iU,jChoVec)
C    *                      *Cho_Ket(iC,iX,iChoVec)
C                 End Do
C               End Do
              END DO
            END DO
            Call DGEMM_('N','N',nSsh(iSym),nCho,NA,
     *                  1.0D+00,WRK2,nSsh(iSym),
     *                          Cho_Ket(1,iX,1),NC*NX,
     *                  1.0D+00,BraSA(1,iUabs,1),nSsh(iSym)*nAsh(iSym))
            Call DGEMM_('N','N',nSsh(iSym),nCho,NA,
     *                  1.0D+00,WRK2,nSsh(iSym),Cho_Bra(1,iU,1),NA*NU,
     *                  0.0D+00,Work(ipWRK),nSsh(iSym))
            Do iChoVec = 1, nCho
              Call DaXpY_(nSsh(iSym),1.0D+00,
     *                    Work(ipWRK+nSsh(iSym)*(iChoVec-1)),1,
     *                    BraSA(1,iXabs,iChoVec),1)
            End Do
C
            Call DGEMM_('N','N',NC,nCho,NA,
     *                  1.0D+00,WRK1,NC,Cho_Bra(1,iU,1),NA*NU,
     *                  0.0D+00,Work(ipWRK),NC)
            Call DGEMM_('T','N',nCho,nCho,NC,
     *                  2.0D+00,Work(ipWRK),NC,
     *                          Cho_Ket(1,iX,1),NC*NX,
     *                  1.0D+00,A_PT2,nCho)
          END DO
        END DO
C
        CALL RHS_FREE(nASP,nISP,ipTP)
      END IF
C
C     ---- FM
C
      IF (NWFM.GT.0.AND.NINDEP(ISYM,9).GT.0) THEN
        !! Read the T-amplitude
        ICASE=9
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM.ne.0) Then
          nISM = nISup(iSym,iCase)
          nVec = nINM*nISM
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASM,nISM,ipTM)
            CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
          End If
        End If
C
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          IXMAX=NX
          IF(ISYU.EQ.ISYX) IXMAX=IU-1
          DO IX=1,IXMAX
            IXABS=IX+NAES(ISYX)
            IW1=KTGTU(IUABS,IXABS)-NTGTUES(ISYM)
C
            Call DCopy_(NA*NC,[0.0D+00],0,WRK1,1)
            Call DCopy_(nSsh(iSym)*nSsh(iSym),[0.0D+00],0,WRK2,1)
            Call DCopy_(nSsh(iSym)*nSsh(iSym),[0.0D+00],0,Work(ipWRK),1)
            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              DO IC=1,NC
                ICABS=IC+NSES(ISYC)
                IF(IAABS.GT.ICABS) THEN
                  IW2=KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
                  SCL = -0.5D+00
                ELSE IF(IAABS.LT.ICABS) THEN
                  IW2=KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
                  SCL =  0.5D+00
                ELSE
                  CYCLE
                END IF
                IW=IW1+NASM*(IW2-1)
C
                ValFM = SCL*Work(ipTM+IW-1)*2.0D+00
                WRK1(iC+NC*(iA-1)) = ValFM
                WRK2(iCabs+nSsh(iSym)*(iA-1)) = ValFM
                Work(ipWRK+iC-1+NC*(iAabs-1)) = ValFM
C               Do iChoVec = 1, nCho
C                 BraSA(iCabs,iXabs,iChoVec)
C    *              = BraSA(iCabs,iXabs,iChoVec)
C    *              + ValFM*Cho_Bra(iA,iU,iChoVec)
C                 BraSA(iAabs,iUabs,iChoVec)
C    *              = BraSA(iAabs,iUabs,iChoVec)
C    *              + ValFM*Cho_Ket(iC,iX,iChoVec)
C                 Do jChoVec = 1, nCho
C                   A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *               + ValFM*Cho_Bra(iA,iU,iChoVec)
C    *                      *Cho_Ket(iC,iX,jChoVec)
C    *               + ValFM*Cho_Bra(iA,iU,jChoVec)
C    *                      *Cho_Ket(iC,iX,iChoVec)
C                 End Do
C               End Do
              END DO
            END DO
C
            Call DGEMM_('T','N',nSsh(iSym),nCho,NC,
     *                  1.0D+00,Work(ipWRK),NC,
     *                          Cho_Ket(1,iX,1),NC*NX,
     *                  1.0D+00,BraSA(1,iUabs,1),nSsh(iSym)*nAsh(iSym))
            Call DGEMM_('N','N',nSsh(iSym),nCho,NA,
     *                  1.0D+00,WRK2,nSsh(iSym),Cho_Bra(1,iU,1),NA*NU,
     *                  0.0D+00,Work(ipWRK),nSsh(iSym))
            Do iChoVec = 1, nCho
              Call DaXpY_(nSsh(iSym),1.0D+00,
     *                    Work(ipWRK+nSsh(iSym)*(iChoVec-1)),1,
     *                    BraSA(1,iXabs,iChoVec),1)
            End Do
C
            Call DGEMM_('N','N',NC,nCho,NA,
     *                  1.0D+00,WRK1,NC,Cho_Bra(1,iU,1),NA*NU,
     *                  0.0D+00,Work(ipWRK),NC)
            Call DGEMM_('T','N',nCho,nCho,NC,
     *                  2.0D+00,Work(ipWRK),NC,
     *                          Cho_Ket(1,iX,1),NC*NX,
     *                  1.0D+00,A_PT2,nCho)
          END DO
        END DO
C
        CALL RHS_FREE(nASM,nISM,ipTM)
      END IF
C
      RETURN
C
      End Subroutine OLagNS_RI_F
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_G(NA,NU,NC,NL,NAUCL,
     &                       Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     DIMENSION Buff(nBuff)
C     DIMENSION idxBuf(nBuff)
C     DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC*NL,NCHO)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NL,NCHO)
*      Logical Incore
      DIMENSION IOFF1(8),IOFF2(8)
C
      ISYU = ISYI
      ISYL = ISYK
C
      ISYA=MUL(JSYM,ISYU)
      ISYC=MUL(JSYM,ISYL)
      ISYM=ISYU
      ISYAC=MUL(ISYA,ISYC)
C Set up offset table:
      IO1=0
      IO2=0
      DO ISI=1,NSYM
        IOFF1(ISI)=IO1
        IOFF2(ISI)=IO2
        ISAB=MUL(ISI,ISYM)
        IO1=IO1+NISH(ISI)*NAGEB(ISAB)
        IO2=IO2+NISH(ISI)*NAGTB(ISAB)
      END DO

C   Allocate W with parts WP,WM
      NAS=NASH(ISYM)
      NISP=NISUP(ISYM,10)
      NISM=NISUP(ISYM,11)
      NIS=NISP+NISM
      NWGP=NAS*NISP
      NWGM=NAS*NISM
      NWG=NWGP+NWGM
C
      LDGP=NAS
      LDGM=NAS
C
C     ---- GP
C
      IF (NWGP.GT.0) THEN
        !! Read the T-amplitude
        ICASE=10
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP.ne.0) Then
          nISP = nISup(iSym,iCase)
          nVec = nINP*nISP
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASP,nISP,ipTP)
            CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
          End If
        End If
C
        NBXSZC=NSECBX
        KCL=NAUCL/(NA*NU)
        NBXSZL=KCL/NC
        IF (NBXSZL.LE.0) THEN
          Write (6,*) 'Not enough memory in ADDRHSG, I give up'
          CALL Abend()
        ENDIF
C
        DO ICSTA=1,NC,NBXSZC
          ICEND=MIN(ICSTA-1+NBXSZC,NC)
C         NCSZ=ICEND-ICSTA+1
          DO ILSTA=1,NL,NBXSZL
            ILEND=MIN(ILSTA-1+NBXSZL,NL)
C           NLSZ=ILEND-ILSTA+1
C
        DO IL=ILSTA,ILEND
          ILABS=IL+NIES(ISYL)
          DO IC=ICSTA,ICEND
            ICABS=IC+NSES(ISYC)
C
            Call DCopy_(NA*NU,[0.0D+00],0,WRK1,1)
            Call DCopy_(nSsh(iSym)*nAsh(iSym),[0.0D+00],0,WRK2,1)
            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              SCL=SQRT(0.5D0)
              IF(IAABS.GE.ICABS) THEN
               IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
               IF(IAABS.EQ.ICABS) SCL=1.0D0
              ELSE
               IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
              END IF
              DO IU=1,NU
                IUABS=IU+NAES(ISYU)
                IW1=IU
                IW2=IL+NL*(IAGEC-1)+IOFF1(ISYL)
                IW=IW1+NAS*(IW2-1)
C
                ValGP = SCL*Work(ipTP+IW-1)*2.0d+00
                WRK1(iA+NA*(iU-1)) = ValGP
                WRK2(iAabs+nSsh(iSym)*(iUabs-1)) = ValGP
C               Do iChoVec = 1, nCho
C                 BraSA(iAabs,iUabs,iChoVec)
C    *              = BraSA(iAabs,iUabs,iChoVec)
C    *              + ValGP*Cho_Ket(iC,iL,iChoVec)
C                 BraSI(iCabs,iLabs,iChoVec)
C    *              = BraSI(iCabs,iLabs,iChoVec)
C    *              + ValGP*Cho_Bra(iA,iU,iChoVec)
C                 Do jChoVec = 1, nCho
C                   A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *               + ValGP*Cho_Bra(iA,iU,iChoVec)
C    *                      *Cho_Ket(iC,iL,jChoVec)
C    *               + ValGP*Cho_Bra(iA,iU,jChoVec)
C    *                      *Cho_Ket(iC,iL,iChoVec)
C                 End Do
C               End Do
              END DO
            END DO
            Call DGEMM_('N','N',nSsh(iSym)*nAsh(iSym),nCho,1,
     *                  1.0D+00,WRK2,nSsh(iSym)*nAsh(iSym),
     *                          Cho_Ket(iC,iL,1),nSsh(iSym)*nIsh(iSym),
     *                  1.0D+00,BraSA(1,1,1),nSsh(iSym)*nAsh(iSym))
            Call DGEMV_('T',NA*NU,nCho,
     *                  1.0D+00,Cho_Bra,NA*NU,WRK1,1,
     *                  0.0D+00,WRK2,1)
            Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                  BraSI(iCabs,iLabs,1),nSsh(iSym)*nIsh(iSym))
            Call DGEMM_('T','N',nCho,nCho,1,
     *                  2.0D+00,Cho_Ket(iC,iL,1),NC*NL,WRK2,1,
     *                  1.0D+00,A_PT2,nCho)
          END DO
        END DO
C
          ENDDO
        ENDDO
C
        CALL RHS_FREE(nASP,nISP,ipTP)
      END IF
C
C     ---- GM
C
      IF (NWGM.GT.0) THEN
        !! Read the T-amplitude
        ICASE=11
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM.ne.0) Then
          nISM = nISup(iSym,iCase)
          nVec = nINM*nISM
          If (nVec.ne.0) Then
            Call RHS_ALLO(nASM,nISM,ipTM)
            CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
          End If
        End If
C
        NBXSZC=NSECBX
        KCL=NAUCL/(NA*NU)
        NBXSZL=KCL/NC
        IF (NBXSZL.LE.0) THEN
          Write (6,*) 'Not enough memory in ADDRHSG, I give up'
          CALL Abend()
        ENDIF


        DO ICSTA=1,NC,NBXSZC
          ICEND=MIN(ICSTA-1+NBXSZC,NC)
          NCSZ=ICEND-ICSTA+1
          DO ILSTA=1,NL,NBXSZL
            ILEND=MIN(ILSTA-1+NBXSZL,NL)
            NLSZ=ILEND-ILSTA+1
C
        ICL=0
        IBUF=0
        DO IL=ILSTA,ILEND
          ILABS=IL+NIES(ISYL)
          DO IC=ICSTA,ICEND
            ICABS=IC+NSES(ISYC)
            ICL=ICL+1
C
            Call DCopy_(NA*NU,[0.0D+00],0,WRK1,1)
            Call DCopy_(nSsh(iSym)*nAsh(iSym),[0.0D+00],0,WRK2,1)
            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              IF(IAABS.GT.ICABS) THEN
                IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                SCL=SQRT(1.5D0)
              ELSE IF (IAABS.LT.ICABS) Then
                IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL=-SQRT(1.5D0)
              ELSE
                CYCLE
              End If
              DO IU=1,NU
                IUABS=IU+NAES(ISYU)
                IW1=IU
                IW2=IL+NL*(IAGTC-1)+IOFF2(ISYL)
                IW=IW1+NAS*(IW2-1)
C
                ValGM = SCL*Work(ipTM+IW-1)*2.0d+00
                WRK1(iA+NA*(iU-1)) = ValGM
                WRK2(iAabs+nSsh(iSym)*(iUabs-1)) = ValGM
C               Do iChoVec = 1, nCho
C                 BraSA(iAabs,iUabs,iChoVec)
C    *              = BraSA(iAabs,iUabs,iChoVec)
C    *              + ValGM*Cho_Ket(iC,iL,iChoVec)
C                 BraSI(iCabs,iLabs,iChoVec)
C    *              = BraSI(iCabs,iLabs,iChoVec)
C    *              + ValGM*Cho_Bra(iA,iU,iChoVec)
C                 Do jChoVec = 1, nCho
C                  A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *               + ValGM*Cho_Bra(iA,iU,iChoVec)
C    *                      *Cho_Ket(iC,iL,jChoVec)
C    *               + ValGM*Cho_Bra(iA,iU,jChoVec)
C    *                      *Cho_Ket(iC,iL,iChoVec)
C                 End Do
C               End Do
              END DO
            END DO
            Call DGEMM_('N','N',nSsh(iSym)*nAsh(iSym),nCho,1,
     *                  1.0D+00,WRK2,nSsh(iSym)*nAsh(iSym),
     *                          Cho_Ket(iC,iL,1),nSsh(iSym)*nIsh(iSym),
     *                  1.0D+00,BraSA(1,1,1),nSsh(iSym)*nAsh(iSym))
            Call DGEMV_('T',NA*NU,nCho,
     *                  1.0D+00,Cho_Bra,NA*NU,WRK1,1,
     *                  0.0D+00,WRK2,1)
            Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                  BraSI(iCabs,iLabs,1),nSsh(iSym)*nIsh(iSym))
            Call DGEMM_('T','N',nCho,nCho,1,
     *                  2.0D+00,Cho_Ket(iC,iL,1),NC*NL,WRK2,1,
     *                  1.0D+00,A_PT2,nCho)
          END DO
        END DO

          ENDDO
        ENDDO

        CALL RHS_FREE(nASM,nISM,ipTM)
      END IF
C
      RETURN
C
      End Subroutine OLagNS_RI_G
C
C-----------------------------------------------------------------------
C
      !! ADDRHSH
      Subroutine OLagNS_RI_H(NA,NJ,NC,NL,NAJCL,
     &                   Cho_Bra,Cho_Ket,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     DIMENSION AJCL(NC*NL,*)
      DIMENSION Cho_Bra(NA,NJ,NCHO), Cho_Ket(NC,NL,NCHO)
C
      ISYJ = ISYI
      ISYL = ISYK
C
      IF(ISYJ.LT.ISYL) Return
      ISYA=MUL(JSYM,ISYJ)
      ISYC=MUL(JSYM,ISYL)
      ISYM=MUL(ISYA,ISYC)
      ISYAC=ISYM
      ISYJL=ISYM
C
      NASP=NAGEB(ISYM)
      NISP=NIGEJ(ISYM)
      NWHP=NASP*NISP
      IF(NWHP.EQ.0) Return
      if (nwhp.eq.0) write (6,*) cho_bra(1,1,1) !! avoid unused tenta
      NASM=NAGTB(ISYM)
      NISM=NIGTJ(ISYM)
C     NWHM=NASM*NISM
C
C     LDHP=NASP
C     LDHM=NASM
C
C     ---- HP
C
      !! Read the T-amplitude
      ICASE=12
      nINP = nINDEP(iSym,iCase)
      nASP = nASup(iSym,iCase)
      If (nINP.ne.0) Then
        nISP = nISup(iSym,iCase)
        nVec = nINP*nISP
        If (nVec.ne.0) Then
          Call RHS_ALLO(nINP,nISP,ipTP)
          Call RHS_READ_SR(ipTP,iCase,iSym,iVec)
        End If
      End If
C
      NBXSZA=NSECBX
      KAJ=NAJCL/(NC*NL)
      NBXSZJ=KAJ/NA
      IF (NBXSZJ.LE.0) THEN
        Write (6,*) 'Not enough memory in ADDRHSH, I give up'
        CALL Abend()
      ENDIF
C
      NBXSZC=NSECBX
      NBXSZL=NINABX
C
      DO IASTA=1,NA,NBXSZA
        IAEND=MIN(IASTA-1+NBXSZA,NA)
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          DO ICSTA=1,NC,NBXSZC
            ICEND=MIN(ICSTA-1+NBXSZC,NC)
            DO ILSTA=1,NL,NBXSZL
              ILEND=MIN(ILSTA-1+NBXSZL,NL)
C
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        ILMAX=NL
        IF(ISYJ.EQ.ISYL) ILMAX=IJ
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
C
          Call DCopy_(NC*NL,[0.0D+00],0,WRK1,1)
          Call DCopy_(nSsh(iSym)*nIsh(iSym),[0.0D+00],0,WRK2,1)
          DO IL=ILSTA,MIN(ILEND,ILMAX)
            ILABS=IL+NIES(ISYL)
            SCL1=1.0D0
            IJGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
            IF(IJABS.EQ.ILABS) SCL1=SQRT(0.5D0)
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
C
              SCL=SCL1
              IF(IAABS.GE.ICABS) THEN
                IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                IF(IAABS.EQ.ICABS) SCL=SQRT(2.0D0)*SCL1
              ELSE
                IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
              END IF
              IW=IAGEC+NAGEB(ISYM)*(IJGEL-1)
C
              ValHP = SCL*Work(ipTP+IW-1)*2.0d+00
              WRK1(iC+NC*(iL-1)) = ValHP
              WRK2(iCabs+nSsh(iSym)*(iLabs-1)) = ValHP
C             Do iChoVec = 1, nCho
C               BraSI(iAabs,iJabs,iChoVec) = BraSI(iAabs,iJabs,iChoVec)
C    *            + ValHP*Cho_Ket(iC,iL,iChoVec)
C               BraSI(iCabs,iLabs,iChoVec) = BraSI(iCabs,iLabs,iChoVec)
C    *            + ValHP*Cho_Ket(iA,iJ,iChoVec)
C               Do jChoVec = 1, nCho
C                 A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *             + ValHP*Cho_Ket(iA,iJ,iChoVec)*Cho_Bra(iC,iL,jChoVec)
C    *             + ValHP*Cho_Ket(iA,iJ,jChoVec)*Cho_Bra(iC,iL,iChoVec)
C               End Do
C             End Do
            END DO
          END DO
          Call DGEMM_('N','N',nSsh(iSym)*nIsh(iSym),nCho,1,
     *                1.0D+00,WRK2,nSsh(iSym)*nIsh(iSym),
     *                        Cho_Ket(iA,iJ,1),nSsh(iSym)*nIsh(iSym),
     *                1.0D+00,BraSI(1,1,1),nSsh(iSym)*nIsh(iSym))
          Call DGEMV_('T',NC*NL,nCho,
     *                1.0D+00,Cho_Ket,NC*NL,WRK1,1,
     *                0.0D+00,WRK2,1)
          Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                BraSI(iAabs,iJabs,1),nSsh(iSym)*nIsh(iSym))
          Call DGEMM_('T','N',nCho,nCho,1,
     *                2.0D+00,Cho_Ket(iA,iJ,1),NA*NJ,WRK2,1,
     *                1.0D+00,A_PT2,nCho)
        END DO
      END DO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      CALL RHS_FREE(nINP,nISP,ipTP)
C
C     ---- HM
C
      !! Read the T-amplitude
      ICASE=13
      nINM = nINDEP(iSym,iCase)
      nASM = nASup(iSym,iCase)
      If (nINM.ne.0) Then
        nISM = nISup(iSym,iCase)
        nVec = nINM*nISM
        If (nVec.ne.0) Then
          Call RHS_ALLO(nINM,nISM,ipTM)
          Call RHS_READ_SR(ipTM,iCase,iSym,iVec)
        End If
      End If
C
      NBXSZA=NSECBX
      KAJ=NAJCL/(NC*NL)
      NBXSZJ=KAJ/NA
      IF (NBXSZJ.LE.0) THEN
        Write (6,*) 'Not enough memory in ADDRHSH, I give up'
        CALL Abend()
      ENDIF
      NBXSZC=NSECBX
      NBXSZL=NINABX
C
      DO IASTA=1,NA,NBXSZA
        IAEND=MIN(IASTA-1+NBXSZA,NA)
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          DO ICSTA=1,NC,NBXSZC
            ICEND=MIN(ICSTA-1+NBXSZC,NC)
            DO ILSTA=1,NL,NBXSZL
              ILEND=MIN(ILSTA-1+NBXSZL,NL)
C
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        ILMAX=NL
        IF(ISYJ.EQ.ISYL) ILMAX=IJ-1
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
C
          Call DCopy_(NC*NL,[0.0D+00],0,WRK1,1)
          Call DCopy_(nSsh(iSym)*nIsh(iSym),[0.0D+00],0,WRK2,1)
          DO IL=ILSTA,MIN(ILMAX,ILEND)
            ILABS=IL+NIES(ISYL)
            IJGTL=KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
C
              IF (IAABS.GT.ICABS) THEN
                IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                SCL= SQRT(3.0D0)
              ELSE IF(IAABS.LT.ICABS) THEN
                IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL=-SQRT(3.0D0)
              ELSE
                Cycle
              ENDIF
              IW=IAGTC+NAGTB(ISYM)*(IJGTL-1)
C
              ValHM = SCL*Work(ipTM+IW-1)*2.0d+00
              WRK1(iC+NC*(iL-1)) = ValHM
              WRK2(iCabs+nSsh(iSym)*(iLabs-1)) = ValHM
C             Do iChoVec = 1, nCho
C               BraSI(iAabs,iJabs,iChoVec) = BraSI(iAabs,iJabs,iChoVec)
C    *            + ValHM*Cho_Ket(iC,iL,iChoVec)
C               BraSI(iCabs,iLabs,iChoVec) = BraSI(iCabs,iLabs,iChoVec)
C    *            + ValHM*Cho_Ket(iA,iJ,iChoVec)
C               Do jChoVec = 1, nCho
C                 A_PT2(iChoVec,jChoVec) = A_PT2(iChoVec,jChoVec)
C    *             + ValHM*Cho_Ket(iA,iJ,iChoVec)*Cho_Bra(iC,iL,jChoVec)
C    *             + ValHM*Cho_Ket(iA,iJ,jChoVec)*Cho_Bra(iC,iL,iChoVec)
C               End Do
C             End Do
            END DO
          END DO
          Call DGEMM_('N','N',nSsh(iSym)*nIsh(iSym),nCho,1,
     *                1.0D+00,WRK2,nSsh(iSym)*nIsh(iSym),
     *                        Cho_Ket(iA,iJ,1),nSsh(iSym)*nIsh(iSym),
     *                1.0D+00,BraSI(1,1,1),nSsh(iSym)*nIsh(iSym))
          Call DGEMV_('T',NC*NL,nCho,
     *                1.0D+00,Cho_Ket,NC*NL,WRK1,1,
     *                0.0D+00,WRK2,1)
          Call DaXpY_(nCho,1.0D+00,WRK2,1,
     *                BraSI(iAabs,iJabs,1),nSsh(iSym)*nIsh(iSym))
          Call DGEMM_('T','N',nCho,nCho,1,
     *                2.0D+00,Cho_Ket(iA,iJ,1),NA*NJ,WRK2,1,
     *                1.0D+00,A_PT2,nCho)
        END DO
      END DO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      CALL RHS_FREE(nINM,nISM,ipTM)
C
      Return
C
      End Subroutine OLagNS_RI_H
C
      End Subroutine OLagNS_RI
