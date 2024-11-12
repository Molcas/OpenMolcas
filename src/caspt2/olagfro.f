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
      Subroutine OLagFro0(DPT2_ori,DPT2)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension DPT2_ori(*),DPT2(*)
C
      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI1.gt.0) Then
          nFroI = nFro(iSym)
          ! nIshI = nIsh(iSym)
          ! nAshI = nAsh(iSym)
          ! nSshI = nSsh(iSym)
C
          !! Do for all orbitals
          Do iOrb = 1, nOrbI1
            iOrb1 = iOrb
            iOrb2 = iOrb+nFroI
            Do jOrb = 1, nOrbI1
              jOrb1 = jOrb
              jOrb2 = jOrb+nFroI
              DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *          = DPT2_ori(iMO1+iOrb1-1+nOrbI1*(jOrb1-1))
              DPT2(iMO2+jOrb2-1+nOrbI2*(iOrb2-1))
     *          = DPT2_ori(iMO1+jOrb1-1+nOrbI1*(iOrb1-1))
            End Do
          End Do
C
          !! Inactive orbitals
        ! Do iOrb = 1, nIshI
        !   iOrb1 = iOrb
        !   iOrb2 = iOrb+nFroI
        !   Do jOrb = 1, nIshI
        !     jOrb1 = jOrb
        !     jOrb2 = jOrb+nFroI
        !     DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *  !       = DPT2_ori(iMO1+iOrb1-1+nOrbI1*(jOrb1-1))
        !     DPT2(iMO2+jOrb2-1+nOrbI2*(iOrb2-1))
     *  !       = DPT2_ori(iMO1+jOrb1-1+nOrbI1*(iOrb1-1))
        !   End Do
        ! End Do
C
        ! !! External orbitals
        ! Do iOrb = 1, nSshI
        !   iOrb1 = iOrb+nIshI+nAshI
        !   iOrb2 = iOrb+nFroI+nIshI+nAshI
        !   Do jOrb = 1, nSshI
        !     jOrb1 = jOrb+nIshI+nAshI
        !     jOrb2 = jOrb+nFroI+nIshI+nAshI
        !     DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *  !       = DPT2_ori(iMO1+iOrb1-1+nOrbI1*(jOrb1-1))
        !     DPT2(iMO2+jOrb2-1+nOrbI2*(iOrb2-1))
     *  !       = DPT2_ori(iMO1+jOrb1-1+nOrbI1*(iOrb1-1))
        !   End Do
        ! End Do
        End If
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do
C
      End Subroutine OLagFro0
C
C-----------------------------------------------------------------------
C
      Subroutine OLagFroD(DIA,DI,RDMSA,Trf)
C
      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension DIA(*),DI(*)
      Dimension RDMSA(*),Trf(*)
      real(kind=wp),allocatable :: WRK1(:),WRK2(:)
C
      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')
C
      iAOtr = 0
      iAOsq = 1
      Do iSym = 1, nSym
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nBasI = nBas(iSym)
        nCorI = nFroI + nIshI
C
        !! full density matrix
      ! Call SQUARE(WRK1(1+iAOtr),DIA(iAOsq),1,nBasI,nBasI)
      ! !! off-diagonal elements have to be halved
      ! Do Mu = 1, nBasI
      !   Do Nu = 1, nBasI
      !     If (Mu.eq.Nu) Cycle
      !     DIA(iAOsq+Mu-1+nBasI*(Nu-1))
     *!       = 0.5D+00*DIA(iAOsq+Mu-1+nBasI*(Nu-1))
      !   End Do
      ! End Do
C
        !! inactive density matrix
        Call DGEMM_('N','T',nBasI,nBasI,nCorI,
     *              2.0D+00,CMOPT2,nBasI,CMOPT2,nBasI,
     *              0.0D+00,DI(iAOsq),nBasI)
C
        !! inactive+active density matrix
        !! Somehow, the above density matrix obtained by calling
        !! Get_D1AO is incorrect... at least, cannot be used.
        ! 1) inactive part
        Call DCopy_(nBasI**2,DI,1,DIA,1)
        ! 2)  RDMSA is defined in CASSCF orbitals, so transform RDMSA to
        !     CASPT2 orbital basis
        Call DGemm_('T','N',nAshI,nAshI,nAshI,
     *              1.0D+00,Trf(1+nCorI+nBasI*nCorI),nBasI,RDMSA,nAshI,
     *              0.0D+00,WRK2,nAshI)
        Call DGemm_('N','N',nAshI,nAshI,nAshI,
     *              1.0D+00,WRK2,nAshI,
     *                      Trf(1+nCorI+nBasI*nCorI),nBasI,
     *              0.0D+00,WRK1,nAshI)
        ! 3) Finally, add the active part
        Call DGemm_('N','N',nBasI,nAshI,nAshI,
     *              1.0D+00,CMOPT2(1+nBasI*nCorI),nBasI,
     *                      WRK1,nAshI,
     *              0.0D+00,WRK2,nBasI)
        Call DGemm_('N','T',nBasI,nBasI,nAshI,
     *              1.0D+00,WRK2,nBasI,
     *                      CMOPT2(1+nBasI*nCorI),nBasI,
     *              1.0D+00,DIA,nBasI)
C
        iAOtr = iAOtr + nBasI*(nBasI+1)/2
        iAOsq = iAOsq + nBasI*nBasI
      End Do
C
      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
C
      Return
C
      End Subroutine OLagFroD
C
C-----------------------------------------------------------------------
C
      Subroutine OLagFro1(DPT2,OLag)
C
      use caspt2_global, only: FIFA_all
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension DPT2(*),OLag(*)
C
      iMO  = 1
      DO iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym)
        nFroI = nFro(iSym)
        If (nOrbI.gt.0.and.nFroI.gt.0) Then
          nIshI = nIsh(iSym)
          nBasI = nBas(iSym)
          !! Make sure that the frozen orbital of the orbital Lagrangian
          !! is zero
          Call DCopy_(nOrbI*nFroI,[0.0D+00],0,OLag,1)
          Do iOrb = 1, nFroI
            Do jOrb = nFroI+1, nFroI+nIshI
C         write(6,*) iorb,jorb,OLag(iMO+iOrb-1+nOrbI*(jOrb-1))
              Tmp = -0.5D+00*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))
     *                       -OLag(iMO+jOrb-1+nOrbI*(iOrb-1)))
     *            /(FIFA_all(iOrb+nBasI*(iOrb-1))
     *             -FIFA_all(jOrb+nBasI*(jOrb-1)))
C         write(6,*) tmp
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     *          = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) + Tmp
              DPT2(iMO+jOrb-1+nOrbI*(iOrb-1))
     *          = DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) + Tmp
            End Do
          End Do
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do
C     write(6,*) "DPT2 after frozen orbital"
C     call sqprt(dpt2,nbast)
C
      End Subroutine OLagFro1
C
C-----------------------------------------------------------------------
C
      Subroutine OLagFro2(DPT2,FPT2,ERI,Scr)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension DPT2(*),FPT2(*),ERI(*),Scr(*)
C
C     write(6,*) "FPT2 before frozen orbital"
C     call sqprt(fpt2,nbast)
      iMO = 1
      isymi = 1
      isymj = 1
      isyma = 1
      isymb = 1
      DO iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        !! Fpq = ((pq|rs)-1/2(pr|qs))*Drs
        Do iOrb = 1, nFroI
          Do jOrb = nFroI+1, nFroI+nIshI
            Scal = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))!*2.0D+00
            Call Coul(iSymA,iSymI,iSymB,iSymJ,
     *                iOrb,jOrb,
     *                ERI,Scr)
            Call DaXpY_(nOrbI*nOrbI,         Scal,ERI,1,FPT2(iMO),1)
C
            Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *                iOrb,jOrb,
     *                ERI,Scr)
            Call DaXpY_(nOrbI*nOrbI,-0.5D+00*Scal,ERI,1,FPT2(iMO),1)
          End Do
        End Do
C
        !! Symmetrize FPT2
        Do iOrb = 1, nOrbI
          Do jOrb = 1, iOrb-1
            Val = (FPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     *            +FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*0.5D+00
            FPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
            FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
          End Do
        End Do
        iMO = iMO + nOrbI*nOrbI
      End Do
C     write(6,*) "FPT2 after frozen orbital"
C     call sqprt(fpt2,nbast)
C
      End Subroutine OLagFro2
C
C-----------------------------------------------------------------------
C
      Subroutine OLagFro3(FIFA,FIMO,WRK1,WRK2)
C
      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension FIFA(*),FIMO(*),WRK1(*),WRK2(*)
      Character(Len=8) Label
      real(kind=wp),allocatable :: WFLT(:)
C
      !! Read H_{\mu \nu}
      call mma_allocate(WFLT,NBTRI,Label='WFLT')
      IRC=-1
      IOPT=6
      ICOMP=1
      ISYLBL=1
      Label='OneHam  '
      CALL RDONE(IRC,IOPT,Label,ICOMP,WFLT,ISYLBL)
C
      !! AO -> MO transformation
      iAO   = 1
      iAOtr = 1
      iCMO  = 1
      iMO   = 1
      DO iSym = 1, nSym
        nBasI = nBas(iSym)
        nOrbI = nBas(iSym)-nDel(iSym)
C
        !! FIFA
        !! WRK1 = G(D)
        Call DCopy_(nBasI*nBasI,FIFA(iAO),1,WRK1,1)
        !! WRK1 = H+G(D)
        Call Square(WFLT(iAOtr),WRK2,1,nBasI,nBasI)
        Call DaXpY_(nBasI*nBasI,1.0D+00,WRK2,1,WRK1,1)
        !! AO -> MO transformation of H+G(D)
        Call OLagTrf(2,iSym,CMOPT2(iCMO),FIFA(iMO),WRK1,WRK2)
C
        !! FIMO
        !! WRK1 = G(D)
        Call DCopy_(nBasI*nBasI,FIMO(iAO),1,WRK1,1)
C     call docpy_nbasI*nbasi,0.0d+00,0,wrk1,1)
        !! WRK1 = H+G(D)
        Call Square(WFLT(iAOtr),WRK2,1,nBasI,nBasI)
        Call DaXpY_(nBasI*nBasI,1.0D+00,WRK2,1,WRK1,1)
        !! AO -> MO transformation of H+G(D)
        Call OLagTrf(2,iSym,CMOPT2(iCMO),FIMO(iMO),WRK1,WRK2)
C       Do iOrb = 1, nfro(1)+nish(1)
C       scal=4.0d+00
C         Call Coul(iSymA,iSymI,iSymB,iSymJ,
C    *              iOrb,iOrb,
C    *              WRK1,WRK2)
C         Call DaXpY_(nOrbI*nOrbI,         Scal,WRK1,1,FIMO(iMO),1)
C
C         Call Exch(iSymA,iSymI,iSymB,iSymJ,
C    *              iOrb,iOrb,
C    *              WRK1,WRK2)
C         Call DaXpY_(nOrbI*nOrbI,-0.5D+00*Scal,WRK1,1,FIMO(iMO),1)
C       End Do
C
        iAO   = iAO   + nBasI*nBasI
        iAOtr = iAOtr + nBasI*(nBasI+1)/2
        iCMO  = iCMO  + nBasI*nOrbI !?
        iMO   = iMO   + nOrbI*nOrbI
      End Do
C     write(6,*) "FIFA"
C     call sqprt(fifa,nbast)
C     write(6,*) "FIMO"
C     call sqprt(fimo,nbast)
C
      call mma_deallocate(WFLT)
C
      End Subroutine OLagFro3
C
C-----------------------------------------------------------------------
C
      Subroutine OLagFroSq(iSym,Ftr,Fsq)
C
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
C
      Implicit Real*8 (A-H,O-Z)
C
#include "caspt2.fh"
C
      Dimension Ftr(*),Fsq(*)
      real(kind=wp),allocatable :: EPS_loc(:)
C
      call mma_allocate(EPS_loc,nBasT,Label='EPS_loc')
      Call Get_dArray('RASSCF OrbE',EPS_loc,nBasT)
C
      nOrbI = nBas(iSym)-nDel(iSym)
      nFroI = nFro(iSym)
      Call DCopy_(nOrbI**2,[0.0D+00],0,Fsq,1)
C
      !! Frozen orbital
      Do iOrb = 1, nFroI
        Fsq(iOrb+nOrbI*(iOrb-1)) = EPS_loc(iOrb)
      End Do
C
      !! Other orbitals
      NSEQ = 0
      Do iOrb = nFroI+1, nOrbI
        Do jOrb = nFroI+1, iOrb
          NSEQ = NSEQ + 1
          Fsq(iOrb+nOrbI*(jOrb-1)) = Ftr(NSEQ)
          Fsq(jOrb+nOrbI*(iOrb-1)) = Ftr(NSEQ)
        End Do
      End Do
C
      call mma_deallocate(EPS_loc)
C
      End Subroutine OLagFroSq
C
C-----------------------------------------------------------------------
C
      !! focktwo.f
      SUBROUTINE OLagFro4(iSym0,iSymI,iSymJ,iSymK,iSymL0,
     *                    DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,WRK1)

      USE CHOVEC_IO
      use Cholesky, only: InfVec, nDimRS
      use EQSOLV
      use ChoCASPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp

      IMPLICIT REAL*8 (A-H,O-Z)

#include "warnings.h"
#include "caspt2.fh"

      Dimension DPT2AO(*),DPT2CAO(*),FPT2AO(*),FPT2CAO(*),WRK1(*)
      real(kind=wp),allocatable :: CHSPC(:),WRK2(:)
      Integer ISTLT(8),ISTSQ(8),iSkip(8),ipWRK(8)

      integer nnbstr(8,3)

      ! INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)

      iSym = iSym0
      call getritrfinfo(nnbstr,maxvec,n2)

      ISTSQ(1)=0
      ISTLT(1)=0
      Do jSym = 2, nSym
        nB  = nBas(jSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(jSym) = ISTSQ(jSym-1) + nB2
        ISTLT(jSym) = ISTLT(jSym-1) + nB3
      End Do
      Do jSym = 1, nSym
        iSkip(jSym) = 1
        ipWRK(jSym) = 1
      End Do

      nBasI  = nBas(iSymI)
      nBasJ  = nBas(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI.EQ.iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
      If (nBasIJ.eq.0) Return

      nBasK  = nBas(iSymK)
      iSMax  = iSymK
      If (iSymK.EQ.iSymI) iSMax = iSymJ
      iSymL  = 1+iEor(iSymIJ-1,iSymK-1)
      IF (iSymL.GT.iSMax) Return !! should not
      nBasL  = nBas(iSymL0)
      nBasKL = nBasK*nBasL
      IF (iSymK.EQ.iSymL0) nBasKL = (nBasK*(nBasK+1))/2
      If (nBasKL.eq.0) Return

      call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')

      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym).EQ.0) Return

      ! ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      ! JRED1=iWork(ipnt)
      ! JRED2=iWork(ipnt-1+NumCho_PT2(iSym))
      JRED1=InfVec(1,2,iSym)
      JRED2=InfVec(NumCho_PT2(iSym),2,iSym)
!     write(6,*) "jred1,jred2 = ", jred1,jred2

* Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED.EQ.0) Cycle

        ILOC=3
        CALL CHO_X_SETRED(IRC,ILOC,JRED)
* For a reduced set, the structure is known, including
* the mapping between reduced index and basis set pairs.
* The reduced set is divided into suitable batches.
* First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.
        ! JEND=JSTART+NVECS_RED-1

* Determine batch length for this reduced set.
* Make sure to use the same formula as in the creation of disk
* address tables, etc, above:
        NBATCH=1+(NVECS_RED-1)/MXNVC

* Loop over IBATCH
        JV1=JSTART
        DO IBATCH=1,NBATCH
C         write(6,*) "ibatch,nbatch = ", ibatch,nbatch
          IBATCH_TOT=IBATCH_TOT+1

          JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
          JV2=JV1+JNUM-1

          JREDC=JRED
* Read a batch of reduced vectors
          CALL CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,
     &                            NUMV,JREDC,MUSED)
C
          IF(NUMV.ne.JNUM) THEN
            write(6,*)' Rats! CHO_VECRD was called, assuming it to'
            write(6,*)' read JNUM vectors. Instead it returned NUMV'
            write(6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
            write(6,*)' Back to the drawing board?'
            CALL QUIT(_RC_INTERNAL_ERROR_)
          END IF
          IF(JREDC.NE.JRED) THEN
            write(6,*)' Rats! It was assumed that the Cholesky vectors'
            write(6,*)' in HALFTRNSF all belonged to a given reduced'
            write(6,*)' set, but they don''t!'
            write(6,*)' JRED, JREDC:',JRED,JREDC
            write(6,*)' Back to the drawing board?'
            write(6,*)' Let the program continue and see what happens.'
          END IF
C
          ipVecL = 1
          Do iVec = 1, NUMV
            !! (strange) reduced form -> squared AO vector (mu nu|iVec)
            jVref = 1 !! only for iSwap=1
C           lscr  = nBasI*(nBasI+1)/2
            ! If (l_NDIMRS.LT.1) Then
            If (size(nDimRS).lt.1) Then
              lscr  = NNBSTR(iSym,3)
            Else
              JREDL = INFVEC(iVec,2,iSym)
              ! lscr  = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
              lscr  = nDimRS(iSym,JREDL)
            End If
            JVEC1 = 1
            iSwap = 2
            WRK2(:) = 0.0d+00
            Call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,jVref,
     *                      JVEC1,1,1,iSym,JREDC,iSwap,ipWRK,WRK2,
     *                      iSkip)
            ipVecL = ipVecL + lscr
C
C           ----- Fock-like transformations -----
C
            Call FDGTRF_RI(WRK2,DPT2AO ,FPT2AO )
            Call FDGTRF_RI(WRK2,DPT2CAO,FPT2CAO)
          End Do
          JV1=JV1+JNUM
        End Do
      End Do
C
      call mma_deallocate(CHSPC)
      call mma_deallocate(WRK2)
C
      !! Have to symmetrize Fock-transformed matrices
      Do i = 1, nBasI
        Do j = 1, i-1
          tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*0.5d+00
          FPT2AO(i+nBasI*(j-1)) = Tmp
          FPT2AO(j+nBasI*(i-1)) = Tmp
          tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*0.5d+00
          FPT2CAO(i+nBasI*(j-1)) = Tmp
          FPT2CAO(j+nBasI*(i-1)) = Tmp
        End Do
      End Do
C
      Return
C
      Contains
C
      Subroutine FDGTRF_RI(ChoVec,DD,FF)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension ChoVec(*),DD(*),FF(*)
C
      !! Coulomb
      Scal = DDot_(nBasI**2,DD,1,ChoVec,1)
      Call DaXpY_(nBasI**2,Scal,ChoVec,1,FF,1)
C
      !! Exchange
      Call DGEMM_('T','N',nBasI,nBasI,nBasI,
     *            1.0D+00,ChoVec,nBasI,DD,nBasI,
     *            0.0D+00,WRK1,nBasI)
      Call DGEMM_('T','T',nBasI,nBasI,nBasI,
     *           -0.5D+00,ChoVec,nBasI,WRK1,nBasI,
     *            1.0D+00,FF,nBasI)
C
      End Subroutine FDGTRF_RI
C
      End Subroutine OLagFro4
