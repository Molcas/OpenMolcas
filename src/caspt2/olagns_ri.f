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
      Subroutine OLagNS_RI(iSym0,DPT2C,DPT2Canti,A_PT2,nChoVec)
C
      Use CHOVEC_IO
      use caspt2_output, only: iPrGlb, verbose
      use caspt2_gradient, only: do_csf
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "warnings.h"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Integer Active, Inactive, Virtual
      Parameter (Inactive=1, Active=2, Virtual=3)
      Integer nSh(8,3)
#ifdef _DEBUGPRINT_
      INTEGER NUMERR
      SAVE NUMERR
      DATA NUMERR / 0 /
#endif
C
      Dimension DPT2C(*),DPT2Canti(*),A_PT2(nChoVec,nChoVec)
C
      Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
      Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
      Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)
C

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'
      END IF
C
      ! iVec = iVecX
      iSym = iSym0
      SCLNEL = 1.0D+00/DBLE(MAX(1,NACTEL))
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
      CALL GETMEM('BRAD','ALLO','REAL',LBRAD,NCHOBUF)
      CALL GETMEM('KETD','ALLO','REAL',LKETD,NCHOBUF)
C
C     Loop over groups of batches of Cholesky vectors
C
      IOFFCV = 1
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
      Call DCopy_(nKet,[0.0D+00],0,Work(LKETD),1)
*                                                                      *
************************************************************************
*                                                                      *
*       Read bra (Cholesky vectors) in the form L(TJ): All symmetries
*
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
      Call DCopy_(nBra,[0.0D+00],0,Work(LBRAD),1)
*                                                                      *
************************************************************************
*                                                                      *
*      Assemble contributions to TJVX
*      Loop over the bras and kets, form <A|0>
*
      Call OLagNS_RI2(Inactive,Active,Active,Active,
     &                'A ',Work(LBRA),Work(LKET),
     &                Work(LBRAD),Work(LKETD))
*                                                                      *
************************************************************************
*                                                                      *
*      TJVL RHSB
*      TJVL: Use TJ buffer as if it was VL, form <B|0>
*
      nKet = nBra
      Call OLagNS_RI2(Inactive,Active,Inactive,Active,
     &                'B ',Work(LBRA),Work(LBRA),
     &                Work(LBRAD),Work(LBRAD))
*                                                                      *
************************************************************************
*                                                                      *
* Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
* We still have L(VX) vectors in core, at WORK(LKETS).
*
      Call Cholesky_Vectors(1,Inactive,Active,JSYM,Work(LBRAD),nBra,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
      Call DCopy_(nBra,[0.0D+00],0,Work(LBRAD),1)
*                                                                      *
************************************************************************
*                                                                      *
* AJVX RHSD1
* Loop over the bra and ket vectors.
*
      Call OLagNS_RI2(Inactive,Virtual,Active,Active,
     &                'D1',Work(LBRA),Work(LKET),
     &                Work(LBRAD),Work(LKETD))
*                                                                      *
************************************************************************
*                                                                      *
* AJCL RHSH
* AJCL: Use AJ buffer still in core as if it was CL, form <H|0>
*
      nKet = nBra
      Call OLagNS_RI2(Inactive,Virtual,Inactive,Virtual,
     &                'H ',Work(LBRA),Work(LBRA),
     &                Work(LBRAD),Work(LBRAD))
*                                                                      *
************************************************************************
*                                                                      *
* Read Bra (Cholesky vectors)= L(AU)
*
      Call Cholesky_Vectors(1,Inactive,Virtual,JSYM,Work(LBRAD),nBra,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Active,Virtual,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
      Call DCopy_(nBra,[0.0D+00],0,Work(LBRAD),1)
C                                                                      *
************************************************************************
*                                                                      *
* AUVX RHSC
* AUVX: Loop over the bras and kets
*
      Call OLagNS_RI2(Active,Virtual,Active,Active,
     &                'C ',Work(LBRA),Work(LKET),
     &                Work(LBRAD),Work(LKETD))
*                                                                      *
************************************************************************
*                                                                      *
* AUCX RHSF
* AUCX: Use AU buffer still in core as if it was CX, form <F|0>
*
      nKet = nBra
      Call OLagNS_RI2(Active,Virtual,Active,Virtual,
     &                'F ',Work(LBRA),Work(LBRA),
     &                Work(LBRAD),Work(LBRAD))
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(VL), all symmetries:
*
      Call Cholesky_Vectors(1,Active,Active,JSYM,Work(LKETD),nKet,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                          Work(LKET),nKet,
     &                          IBSTA,IBEND)
      Call Cholesky_Vectors(2,Inactive,Active,JSYM,Work(LKETD),nKet,
     &                      IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUVL RHSD2
* Loop over bras and kets, form <D2|0>.
*
      Call OLagNS_RI2(Active,Virtual,Inactive,Active,
     &                'D2',Work(LBRA),Work(LKET),
     &                Work(LBRAD),Work(LKETD))
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(CL), all symmetries:
*
      Call Cholesky_Vectors(1,Inactive,Active,JSYM,Work(LKETD),nKet,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                          Work(LKET),nKet,
     &                          IBSTA,IBEND)
      Call Cholesky_Vectors(2,Inactive,Virtual,JSYM,Work(LKETD),nKet,
     &                      IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUCL RHSG
* Loop over bras and kets, form  <G|0>
*
      Call OLagNS_RI2(Active,Virtual,Inactive,Virtual,
     &                'G ',Work(LBRA),Work(LKET),
     &                Work(LBRAD),Work(LKETD))
*                                                                      *
************************************************************************
*                                                                      *
* Read bra vectors AJ
*
      Call Cholesky_Vectors(1,Active,Virtual,JSYM,Work(LBRAD),nBra,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                          Work(LBRA),nBra,
     &                          IBSTA,IBEND)
      Call DCopy_(nBra,Work(LKETD),1,Work(LBRAD),1)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets in the form L(VL)
*
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                          Work(LKET),nKet,
     &                          IBSTA,IBEND)
      Call Cholesky_Vectors(2,Inactive,Active,JSYM,Work(LKETD),nKet,
     &                      IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVL RHSE
* AJVL: Loop over bras and kets. Form <E|0>
*
      Call OLagNS_RI2(Inactive,Virtual,Inactive,Active,
     &                'E ',Work(LBRA),Work(LKET),
     &                Work(LBRAD),Work(LKETD))
*                                                                      *
************************************************************************
*                                                                      *
* End of loop over batches, IB
      Call Cholesky_Vectors(1,Inactive,Virtual,JSYM,Work(LBRAD),nBra,
     &                      IBSTA,IBEND)
      Call Cholesky_Vectors(1,Inactive,Active,JSYM,Work(LKETD),nKet,
     &                      IBSTA,IBEND)
      !! Construct A_PT2
      NVI = NV
      JOFFCV = 1
      DO JBGRP=1,NBGRP
C
        JBSTA=IWORK(LBGRP  +2*(JBGRP-1))
        JBEND=IWORK(LBGRP+1+2*(JBGRP-1))
C
        NVJ=0
        DO JB=JBSTA,JBEND
          NVJ=NVJ+NVLOC_CHOBATCH(JB)
        END DO
C
        !! BraAI
        If (nIsh(iSym0)*nAsh(iSym0).ne.0) Then
        Call Cholesky_Vectors(2,Inactive,Active,JSYM,Work(LBRA),nBra,
     &                        IBSTA,IBEND)
        Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                            Work(LKET),nKet,
     &                            JBSTA,JBEND)
        Call DGEMM_('T','N',NVI,NVJ,nIsh(iSym0)*nAsh(iSym0),
     &              1.0D+00,Work(LBRA),nIsh(iSym0)*nAsh(iSym0),
     &                      Work(LKET),nIsh(iSym0)*nAsh(iSym0),
     &              1.0D+00,A_PT2(IOFFCV,JOFFCV),nChoVec)
        End If
C
        !! BraSI
        If (nIsh(iSym0)*nSsh(iSym0).ne.0) Then
        Call Cholesky_Vectors(2,Inactive,Virtual,JSYM,Work(LBRA),nBra,
     &                        IBSTA,IBEND)
        Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                            Work(LKET),nKet,
     &                            JBSTA,JBEND)
        Call DGEMM_('T','N',NVI,NVJ,nIsh(iSym0)*nSsh(iSym0),
     &              1.0D+00,Work(LBRA),nIsh(iSym0)*nSsh(iSym0),
     &                      Work(LKET),nIsh(iSym0)*nSsh(iSym0),
     &              1.0D+00,A_PT2(IOFFCV,JOFFCV),nChoVec)
        End If
C
        !! BraSA
        If (nAsh(iSym0)*nSsh(iSym0).ne.0) Then
        Call Cholesky_Vectors(2,Active,Virtual,JSYM,Work(LBRA),nBra,
     &                        IBSTA,IBEND)
        Call Get_Cholesky_Vectors(Active,Virtual,JSYM,
     &                            Work(LKET),nKet,
     &                            JBSTA,JBEND)
        Call DGEMM_('T','N',NVI,NVJ,nAsh(iSym0)*nSsh(iSym0),
     &              1.0D+00,Work(LBRA),nAsh(iSym0)*nSsh(iSym0),
     &                      Work(LKET),nAsh(iSym0)*nSsh(iSym0),
     &              1.0D+00,A_PT2(IOFFCV,JOFFCV),nChoVec)
        End If
C
        !! BraAA
        If (nAsh(iSym0)*nAsh(iSym0).ne.0) Then
        Call Cholesky_Vectors(2,Active,Active,JSYM,Work(LBRA),nBra,
     &                        IBSTA,IBEND)
        Call Get_Cholesky_Vectors(Active,Active,JSYM,
     &                            Work(LKET),nKet,
     &                            JBSTA,JBEND)
        Call DGEMM_('T','N',NVI,NVJ,nAsh(iSym0)*nAsh(iSym0),
     &              1.0D+00,Work(LBRA),nAsh(iSym0)*nAsh(iSym0),
     &                      Work(LKET),nAsh(iSym0)*nAsh(iSym0),
     &              1.0D+00,A_PT2(IOFFCV,JOFFCV),nChoVec)
        End If
        JOFFCV = JOFFCV + NVJ
      END DO !! end of JBGRP loop
      IOFFCV = IOFFCV + NVI
      END DO !! end of IBGRP loop
*                                                                      *
************************************************************************
*                                                                      *
*SVC  Call GetMem('ADDRHS','Free','Real',ipAdd,nAdd)
      CALL GETMEM('BRABUF','FREE','REAL',LBRA,NCHOBUF)
      CALL GETMEM('KETBUF','FREE','REAL',LKET,NCHOBUF)
      CALL GETMEM('BRAD','FREE','REAL',LBRAD,NCHOBUF)
      CALL GETMEM('KETD','FREE','REAL',LKETD,NCHOBUF)
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
C
C
      Call DScal_(nChoVec**2,2.0D+00,A_PT2,1)
C
      If (NBGRP.ne.0) SCLNEL = SCLNEL/DBLE(NBGRP)
      Call DScal_(nBasSq,SCLNEL,DPT2C,1)
      If (do_csf) Call DScal_(nBasSq,SCLNEL,DPT2Canti,1)
C
C
      Return
C
      Contains
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_RI2(ITI,ITP,ITK,ITQ,Case,Cho_Bra,Cho_Ket,
     &                      Cho_BraD,Cho_KetD)
C
      use caspt2_output, only:iPrGlb,debug
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
C     DIMENSION Cho_Bra(nBra), Cho_Ket(nKet)
      DIMENSION Cho_Bra(*), Cho_Ket(*)
      DIMENSION Cho_BraD(*), Cho_KetD(*)
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
             CALL OLagNS_RI_A(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'B ') Then
             CALL OLagNS_RI_B(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'D1') Then
             CALL OLagNS_RI_D1(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'H ') Then
             CALL OLagNS_RI_H(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'C ') Then
             CALL OLagNS_RI_C(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'F ') Then
             CALL OLagNS_RI_F(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'D2') Then
             CALL OLagNS_RI_D2(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'G ') Then
             CALL OLagNS_RI_G(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case.eq.'E ') Then
             CALL OLagNS_RI_E(ISYI,ISYK,NP,NI,NQ,NK,Work(LPIQK),NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else
             Call Abend()
          End If
C
          LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
      END DO
      CALL CWTime(TotCPU1,TotWall1)
      IF (IPRGLB.GE.VERBOSE) THEN
        write(6,'(" CPU/Wall Time (Case ",A2,"):",2f10.2)')
     *    Case,totcpu1-totcpu0,totwall1-totwall0
      END IF
C
      Return
C
      End Subroutine OLagNS_RI2
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_A(ISYI,ISYK,NT,NJ,NV,NX,TJVX,NTJVX,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
      use caspt2_gradient, only: do_csf
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION TJVX(NT,NJ,NV,NX)
      DIMENSION Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NX,NCHO),
     *          Cho_BraD(NT,NJ,NCHO), Cho_KetD(NV,NX,NCHO)
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
          If (do_csf) Then
            Call RHS_ALLO(nAS,nIS,ipTanti)
            CALL RHS_READ_C(ipTanti,iCase,iSym,7)
          End If
        End If
      End If
C
      Call DCopy_(NTJVX,[0.0D+00],0,TJVX,1)
C
      nOrbT = nFro(iSyT)+nIsh(iSyT)+nAsh(iSyT)+nSsh(iSyT)
      DO IT=1,NT
        ITABS=IT+NAES(ISYT)
        iTtot = iT + nFro(iSyT) + nIsh(iSyT)
        DO IJ=1,NJ
          ! IJABS=IJ+NIES(ISYJ)
          iJtot = iJ + nFro(iSyJ)
C
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
              If (do_csf) Then
                ValAFanti = Work(ipTanti+IW-1)*2.0D+00
                DPT2Canti(iTtot+nOrbT*(iJtot-1))
     *            = DPT2Canti(iTtot+nOrbT*(iJtot-1)) + ValAFanti
              End If
            END IF
            DO IX=1,NX
              IXABS=IX+NAES(ISYX)
              IW1=KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
              IW2=IJ
              IW=IW1+NAS*(IW2-1)
              TJVX(IT,IJ,IV,IX) = Work(ipT+IW-1)
            END DO
          END DO
        END DO
      END DO
C
      Call DGEMM_('T','N',NV*NX,NCHO,NT*NJ,
     *            1.0D+00,TJVX(1,1,1,1),NT*NJ,Cho_Bra(1,1,1),NT*NJ,
     *            1.0D+00,Cho_KetD(1,1,1),NV*NX)
      Call DGEMM_('N','N',NT*NJ,NCHO,NV*NX,
     *            1.0D+00,TJVX(1,1,1,1),NT*NJ,Cho_Ket(1,1,1),NV*NX,
     *            1.0D+00,Cho_BraD(1,1,1),NT*NJ)
C
      CALL RHS_FREE(nAS,nIS,ipT)
      If (do_csf) CALL RHS_FREE(nAS,nIS,ipTanti)
C
      RETURN
C
      End Subroutine OLagNS_RI_A
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_B(ISYI,ISYK,NT,NJ,NV,NL,TJVL,NTJVL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION TJVL(NT,NJ,NV,NL)
      DIMENSION Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NL,NCHO),
     *          Cho_BraD(NT,NJ,NCHO), Cho_KetD(NV,NL,NCHO)
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
      Call DCopy_(NTJVL,[0.0D+00],0,TJVL,1)
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
                TJVL(IT,IJ,IV,IL) = SCL*Work(ipTP+IW-1)
              END DO
            END DO
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
                TJVL(IT,IJ,IV,IL) = TJVL(IT,IJ,IV,IL)
     *            + SCL*Work(ipTM+IW-1)
              END DO
            END DO
          END DO
        END DO
C
        CALL RHS_FREE(nASM,nISM,ipTM)
      END IF
C
      Call DGEMM_('T','N',NV*NL,NCHO,NT*NJ,
     *            1.0D+00,TJVL(1,1,1,1),NT*NJ,Cho_Bra(1,1,1),NT*NJ,
     *            1.0D+00,Cho_KetD(1,1,1),NV*NL)
      Call DGEMM_('N','N',NT*NJ,NCHO,NV*NL,
     *            1.0D+00,TJVL(1,1,1,1),NT*NJ,Cho_Ket(1,1,1),NV*NL,
     *            1.0D+00,Cho_BraD(1,1,1),NT*NJ)
C
      RETURN
C
      End Subroutine OLagNS_RI_B
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_C(ISYI,ISYK,NA,NU,NV,NX,AUVX,NAUVX,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
      use caspt2_gradient, only: do_csf
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AUVX(NA,NU,NV,NX)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NX,NCHO),
     *          Cho_BraD(NA,NU,NCHO), Cho_KetD(NV,NX,NCHO)
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
          If (do_csf) Then
            Call RHS_ALLO(nAS,nIS,ipTanti)
            CALL RHS_READ_C(ipTanti,iCase,iSym,7)
          End If
        End If
      End If
C
      Call DCopy_(NAUVX,[0.0D+00],0,AUVX,1)
C
      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      DO IA=1,NA
        iAtot = iA + nFro(iSyA) + nIsh(iSyA) + nAsh(iSyA)
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          iUtot = iU + nFro(iSyU) + nIsh(iSyU)
          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
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
              If (do_csf) Then
                ValCFanti = Work(ipTanti+IW-1)*2.0D+00
                DPT2Canti(iAtot+nOrbA*(iUtot-1))
     *            = DPT2Canti(iAtot+nOrbA*(iUtot-1)) + ValCFanti
              End If
              ValCF = ValCF*SCLNEL
            END IF
            DO IX=1,NX
              IXABS=IX+NAES(ISYX)
              IW1=KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
              IW2=IA
              IW=IW1+NAS*(IW2-1)
C
              ValC = Work(ipT+IW-1)
              AUVX(IA,IU,IV,IX) = AUVX(IA,IU,IV,IX) + ValC
              AUVX(IA,IX,IU,IX) = AUVX(IA,IX,IU,IX) - ValCF*0.5D+00
            END DO
          END DO
        END DO
      END DO
C
      Call DGEMM_('T','N',NV*NX,NCHO,NA*NU,
     *            1.0D+00,AUVX(1,1,1,1),NA*NU,Cho_Bra(1,1,1),NA*NU,
     *            1.0D+00,Cho_KetD(1,1,1),NV*NX)
      Call DGEMM_('N','N',NA*NU,NCHO,NV*NX,
     *            1.0D+00,AUVX(1,1,1,1),NA*NU,Cho_Ket(1,1,1),NV*NX,
     *            1.0D+00,Cho_BraD(1,1,1),NA*NU)
C
      CALL RHS_FREE(nAS,nIS,ipT)
      If (do_csf) CALL RHS_FREE(nAS,nIS,ipTanti)
C
      RETURN
C
      End Subroutine OLagNS_RI_C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_D1(ISYI,ISYK,NA,NJ,NV,NX,AJVX,NAJVX,
     &                        Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
      use caspt2_gradient, only: do_csf
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AJVX(NV,NX,*)
      DIMENSION Cho_Bra(NA,NJ,NCHO), Cho_Ket(NV,NX,NCHO),
     *          Cho_BraD(NA,NJ,NCHO), Cho_KetD(NV,NX,NCHO)
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
          If (do_csf) Then
            Call RHS_ALLO(nAS,nIS,ipTanti)
            CALL RHS_READ_C(ipTanti,iCase,iSym,7)
          End If
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
      IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
      Call DCopy_(NAJVX,[0.0D+00],0,AJVX,1)
C
      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      IAJ=0
      DO IJ=IJSTA,IJEND
        iJtot = iJ + nFro(iSyJ)
        DO IA=IASTA,IAEND
          iAtot = iA + nFro(iSyA) + nIsh(iSyA) + nAsh(iSyA)
          IAJ = IAJ + 1
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
              ValD = Work(ipT+IW-1)
              If (iVtot.eq.iXtot) Then
                DPT2C(iAtot+nOrbA*(iJtot-1))
     *            = DPT2C(iAtot+nOrbA*(iJtot-1)) + ValD*2.0d+00
                If (do_csf) Then
                  ValDanti = Work(ipTanti+IW-1)*2.0D+00
                  DPT2Canti(iAtot+nOrbA*(iJtot-1))
     *              = DPT2Canti(iAtot+nOrbA*(iJtot-1)) + ValDanti
                End If
              End If
              AJVX(IV,IX,IAJ) = Work(ipT+IW-1)
            END DO
          END DO
        END DO
      END DO
C
      Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NX,
     *            1.0D+00,AJVX(1,1,1),NV*NX,Cho_Ket(1,1,1),NV*NX,
     *            1.0D+00,Cho_BraD(IAJSTA,1,1),NA*NJ)
      Call DGEMM_('N','N',NV*NX,NCHO,NASZ*NJSZ,
     *            1.0D+00,AJVX(1,1,1),NV*NX,Cho_Bra(IAJSTA,1,1),NA*NJ,
     *            1.0D+00,Cho_KetD(1,1,1),NV*NX)
C
        ENDDO
      ENDDO
C
      CALL RHS_FREE(nAS,nIS,ipT)
      If (do_csf) CALL RHS_FREE(nAS,nIS,ipTanti)
C
      RETURN
C
      End Subroutine OLagNS_RI_D1
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_D2(ISYI,ISYK,NA,NU,NV,NL,AUVL,NAUVL,
     &                        Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AUVL(NA,NU,NV,NL)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NL,NCHO),
     *          Cho_BraD(NA,NU,NCHO), Cho_KetD(NV,NL,NCHO)
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
C
      Call DCopy_(NAUVL,[0.0D+00],0,AUVL,1)
C
      DO IA=1,NA
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
            DO IL=1,NL
              IW1=NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
              IW2=IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
              IW=IW1+NAS*(IW2-1)
              AUVL(IA,IU,IV,IL) = Work(ipT+IW-1)
            END DO
          END DO
        END DO
      END DO
C
      Call DGEMM_('T','N',NV*NL,NCHO,NA*NU,
     *            1.0D+00,AUVL,NA*NU,Cho_Bra,NA*NU,
     *            1.0D+00,Cho_KetD,NV*NL)
      Call DGEMM_('N','N',NA*NU,NCHO,NV*NL,
     *            1.0D+00,AUVL,NA*NU,Cho_Ket,NV*NL,
     *            1.0D+00,Cho_BraD,NA*NU)
C
      CALL RHS_FREE(nAS,nIS,ipT)
C
      RETURN
C
      End Subroutine OLagNS_RI_D2
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_E(ISYI,ISYK,NA,NJ,NV,NL,AJVL,NAJVL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AJVL(NV,NL,*)
      DIMENSION Cho_Bra(NA,NJ,NCHO), Cho_Ket(NV,NL,NCHO),
     *          Cho_BraD(NA,NJ,NCHO), Cho_KetD(NV,NL,NCHO)
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
      ! NIS=NISP+NISM
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
        IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        Call DCopy_(NAJVL,[0.0D+00],0,AJVL,1)
C
        IAJ=0
        DO IJ=IJSTA,IJEND
          IJABS=IJ+NIES(ISYJ)
          DO IA=IASTA,IAEND
            ! IAABS=IA+NSES(ISYA)
            IAJ=IAJ+1
C
            DO IV=1,NV
              ! IVABS=IV+NAES(ISYV)
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
                AJVL(IV,IL,IAJ) = SCL*Work(ipTP+IW-1)
              END DO
            END DO
          END DO
        END DO
C
        Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NL,
     *              1.0D+00,AJVL(1,1,1),NV*NL,Cho_Ket(1,1,1),NV*NL,
     *              1.0D+00,Cho_BraD(IAJSTA,1,1),NA*NJ)
        Call DGEMM_('N','N',NV*NL,NCHO,NASZ*NJSZ,
     *              1.0D+00,AJVL(1,1,1),NV*NL,Cho_Bra(IAJSTA,1,1),NA*NJ,
     *              1.0D+00,Cho_KetD(1,1,1),NV*NL)
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
        IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        Call DCopy_(NAJVL,[0.0D+00],0,AJVL,1)
C
        IAJ=0
        DO IJ=IJSTA,IJEND
          IJABS=IJ+NIES(ISYJ)
          DO IA=IASTA,IAEND
            ! IAABS=IA+NSES(ISYA)
            IAJ=IAJ+1
C
            DO IV=1,NV
              ! IVABS=IV+NAES(ISYV)
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
                  AJVL(IV,IL,IAJ) = SCL*Work(ipTM+IW-1)
                END IF
              END DO
            END DO
          END DO
        END DO
C
        Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NL,
     *              1.0D+00,AJVL(1,1,1),NV*NL,Cho_Ket(1,1,1),NV*NL,
     *              1.0D+00,Cho_BraD(IAJSTA,1,1),NA*NJ)
        Call DGEMM_('N','N',NV*NL,NCHO,NASZ*NJSZ,
     *              1.0D+00,AJVL(1,1,1),NV*NL,Cho_Bra(IAJSTA,1,1),NA*NJ,
     *              1.0D+00,Cho_KetD(1,1,1),NV*NL)
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
      SUBROUTINE OLagNS_RI_F(ISYI,ISYK,NA,NU,NC,NX,AUCX,NAUCX,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AUCX(NA,NU,NC,NX)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NX,NCHO),
     *          Cho_BraD(NA,NU,NCHO), Cho_KetD(NC,NX,NCHO)
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
      Call DCopy_(NAUCX,[0.0D+00],0,AUCX,1)
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
                AUCX(IA,IU,IC,IX) = SCL*Work(ipTP+IW-1)
              END DO
            END DO
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
                AUCX(IA,IU,IC,IX) = AUCX(IA,IU,IC,IX)
     *            + SCL*Work(ipTM+IW-1)
              END DO
            END DO
          END DO
        END DO
C
        CALL RHS_FREE(nASM,nISM,ipTM)
      END IF
C
      Call DGEMM_('T','N',NC*NX,NCHO,NA*NU,
     *            1.0D+00,AUCX,NA*NU,Cho_Bra,NA*NU,
     *            1.0D+00,Cho_KetD,NC*NX)
      Call DGEMM_('N','N',NA*NU,NCHO,NC*NX,
     *            1.0D+00,AUCX,NA*NU,Cho_Ket,NC*NX,
     *            1.0D+00,Cho_BraD,NA*NU)
C
      RETURN
C
      End Subroutine OLagNS_RI_F
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_RI_G(ISYI,ISYK,NA,NU,NC,NL,AUCL,NAUCL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AUCL(NA,NU,*)
C     DIMENSION Buff(nBuff)
C     DIMENSION idxBuf(nBuff)
C     DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC*NL,NCHO)
      DIMENSION Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NL,NCHO),
     *          Cho_BraD(NA,NU,NCHO), Cho_KetD(NC,NL,NCHO)
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
      ! NIS=NISP+NISM
      NWGP=NAS*NISP
      NWGM=NAS*NISM
      ! NWG=NWGP+NWGM
C
C     LDGP=NAS
C     LDGM=NAS
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
          NCSZ=ICEND-ICSTA+1
          DO ILSTA=1,NL,NBXSZL
            ILEND=MIN(ILSTA-1+NBXSZL,NL)
            ! NLSZ=ILEND-ILSTA+1
C
            ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
            Call DCopy_(NAUCL,[0.0D+00],0,AUCL,1)
C
        ICL=0
        ! IBUF=0
        DO IL=ILSTA,ILEND
          ! ILABS=IL+NIES(ISYL)
          DO IC=ICSTA,ICEND
            ICABS=IC+NSES(ISYC)
            ICL=ICL+1
C
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
                ! IUABS=IU+NAES(ISYU)
                IW1=IU
                IW2=IL+NL*(IAGEC-1)+IOFF1(ISYL)
                IW=IW1+NAS*(IW2-1)
C
                ValGP = SCL*Work(ipTP+IW-1)
                AUCL(IA,IU,ICL) = ValGP
              END DO
            END DO
          END DO
        END DO
C
        Call DGEMM_('N','N',NA*NU,NCHO,ICL,
     *              1.0D+00,AUCL,NA*NU,
     *                      Cho_Ket(ICLSTA,1,1),NC*NL,
     *              1.0D+00,Cho_BraD(1,1,1),NA*NU)
        Call DGEMM_('T','N',ICL,NCHO,NA*NU,
     *              1.0D+00,AUCL,NA*NU,
     *                      Cho_Bra(1,1,1),NA*NU,
     *              1.0D+00,Cho_KetD(ICLSTA,1,1),NC*NL)
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
            ! NLSZ=ILEND-ILSTA+1
C
            ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
            Call DCopy_(NAUCL,[0.0D+00],0,AUCL,1)
C
        ICL=0
        ! IBUF=0
        DO IL=ILSTA,ILEND
          ! ILABS=IL+NIES(ISYL)
          DO IC=ICSTA,ICEND
            ICABS=IC+NSES(ISYC)
            ICL=ICL+1
C
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
                ! IUABS=IU+NAES(ISYU)
                IW1=IU
                IW2=IL+NL*(IAGTC-1)+IOFF2(ISYL)
                IW=IW1+NAS*(IW2-1)
C
                ValGM = SCL*Work(ipTM+IW-1)
                AUCL(IA,IU,ICL) = ValGM
              END DO
            END DO
          END DO
        END DO
C
        Call DGEMM_('N','N',NA*NU,NCHO,ICL,
     *              1.0D+00,AUCL,NA*NU,
     *                      Cho_Ket(ICLSTA,1,1),NC*NL,
     *              1.0D+00,Cho_BraD(1,1,1),NA*NU)
        Call DGEMM_('T','N',ICL,NCHO,NA*NU,
     *              1.0D+00,AUCL,NA*NU,
     *                      Cho_Bra(1,1,1),NA*NU,
     *              1.0D+00,Cho_KetD(ICLSTA,1,1),NC*NL)
C
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
      Subroutine OLagNS_RI_H(ISYI,ISYK,NA,NJ,NC,NL,AJCL,NAJCL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)
C
      USE SUPERINDEX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AJCL(NC*NL,*)
      DIMENSION Cho_Bra(NA,NJ,NCHO), Cho_Ket(NC,NL,NCHO),
     *          Cho_BraD(NA,NJ,NCHO), Cho_KetD(NC,NL,NCHO)
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
      if (nwhp.eq.0) write (6,*) cho_ketd(1,1,1) !! avoid unused tenta
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
      nVec = 0
      If (nINP.ne.0) Then
        nISP = nISup(iSym,iCase)
        nVec = nINP*nISP
        If (nVec.ne.0) Then
          Call RHS_ALLO(nASP,nISP,ipTP)
          CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
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
        NASZ=IAEND-IASTA+1
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          NJSZ=IJEND-IJSTA+1
C
           IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
           Call DCopy_(NAJCL,[0.0D+00],0,AJCL,1)
C
          DO ICSTA=1,NC,NBXSZC
            ICEND=MIN(ICSTA-1+NBXSZC,NC)
            NCSZ=ICEND-ICSTA+1
            DO ILSTA=1,NL,NBXSZL
              ILEND=MIN(ILSTA-1+NBXSZL,NL)
              ! NLSZ=ILEND-ILSTA+1
C
              ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
C
      IAJ=0
      ! IBUF=0
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        ILMAX=NL
        IF(ISYJ.EQ.ISYL) ILMAX=IJ
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
          IAJ=IAJ+1
C
          ICL=0
          DO IL=ILSTA,MIN(ILEND,ILMAX)
            ILABS=IL+NIES(ISYL)
            SCL1=1.0D0
            IJGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
            IF(IJABS.EQ.ILABS) SCL1=SQRT(0.5D0)
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
              ICL=ICL+1
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
              ValHP = SCL*Work(ipTP+IW-1)
              AJCL(ICLSTA+ICL-1,IAJ) = ValHP
            END DO
          END DO
        END DO
      END DO
C
            ENDDO
          ENDDO
          Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NC*NL,
     *                1.0D+00,AJCL(1,1),NC*NL,Cho_Ket(1,1,1),NC*NL,
     *                1.0D+00,Cho_BraD(IAJSTA,1,1),NA*NJ)
          Call DGEMM_('N','N',NC*NL,NCHO,NASZ*NJSZ,
     *                1.0D+00,AJCL(1,1),NC*NL,Cho_Ket(IAJSTA,1,1),NC*NL,
     *                1.0D+00,Cho_BraD(1,1,1),NA*NJ)
C
        ENDDO
      ENDDO
C
      IF (nVec.ne.0) CALL RHS_FREE(nINP,nISP,ipTP)
C
C     ---- HM
C
      !! Read the T-amplitude
      ICASE=13
      nINM = nINDEP(iSym,iCase)
      nASM = nASup(iSym,iCase)
      nVec = 0
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
        NASZ=IAEND-IASTA+1
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          NJSZ=IJEND-IJSTA+1
C
          IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
          Call DCopy_(NAJCL,[0.0D+00],0,AJCL,1)
C
          DO ICSTA=1,NC,NBXSZC
            ICEND=MIN(ICSTA-1+NBXSZC,NC)
            NCSZ=ICEND-ICSTA+1
            DO ILSTA=1,NL,NBXSZL
              ILEND=MIN(ILSTA-1+NBXSZL,NL)
              ! NLSZ=ILEND-ILSTA+1
C
              ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
C
      IAJ=0
      ! IBUF=0
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        ILMAX=NL
        IF(ISYJ.EQ.ISYL) ILMAX=IJ-1
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
          IAJ=IAJ+1
C
          ICL=0
          DO IL=ILSTA,MIN(ILMAX,ILEND)
            ILABS=IL+NIES(ISYL)
            IJGTL=KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
              ICL=ICL+1
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
              ValHM = SCL*Work(ipTM+IW-1)
              AJCL(ICLSTA+ICL-1,IAJ) = ValHM
            END DO
          END DO
        END DO
      END DO
C
            ENDDO
          ENDDO
          Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NC*NL,
     *                1.0D+00,AJCL(1,1),NC*NL,Cho_Ket(1,1,1),NC*NL,
     *                1.0D+00,Cho_BraD(IAJSTA,1,1),NA*NJ)
          Call DGEMM_('N','N',NC*NL,NCHO,NASZ*NJSZ,
     *                1.0D+00,AJCL(1,1),NC*NL,Cho_Ket(IAJSTA,1,1),NC*NL,
     *                1.0D+00,Cho_BraD(1,1,1),NA*NJ)
        ENDDO
      ENDDO
C
      If (nVec.ne.0) CALL RHS_FREE(nINM,nISM,ipTM)
C
      Return
C
      End Subroutine OLagNS_RI_H
C
      End Subroutine OLagNS_RI
C
C-----------------------------------------------------------------------
C
      Subroutine Cholesky_Vectors(MODE,ITK,ITQ,JSYM,Array,nArray,
     *                            IBSTA,IBEND)
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
      LUCDER = 63 ! tentative
      DO ISYK=1,NSYM
        NQK=NPQ_CHOTYPE(ICASE,ISYK,JSYM)
        IF(NQK.EQ.0) CYCLE
        DO IB=IBSTA,IBEND
          NV=NVLOC_CHOBATCH(IB)
          NKETSM=NQK*NV
          IDISK=IDLOC_CHOGROUP(ICASE,ISYK,JSYM,IB)
          !! MODE = 1: write
          !! MODE = 2: read
          CALL DDAFILE(LUCDER,MODE,Array(LKETSM),NKETSM,IDISK)
          LKETSM=LKETSM+NKETSM
        END DO
      END DO
      nArray=LKETSM-1
*
      Return
      End
