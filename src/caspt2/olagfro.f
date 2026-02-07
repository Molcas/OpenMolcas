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

      use caspt2_module, only: NSYM, NFRO, NORB, NDEL, NBAS
      use definitions, only: wp, iwp

      implicit none

      real(kind=wp), intent(in) :: DPT2_ori(*)
      real(kind=wp), intent(inout) :: DPT2(*)

      integer(kind=iwp) :: iMO1, iMO2, iSym, nOrbI1, nOrbI2, nFroI,
     &  iOrb, iOrb1, iOrb2, jOrb, jOrb1, jOrb2

      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI1 > 0) Then
          nFroI = nFro(iSym)
          !! Do for all orbitals
          Do iOrb = 1, nOrbI1
            iOrb1 = iOrb
            iOrb2 = iOrb+nFroI
            Do jOrb = 1, nOrbI1
              jOrb1 = jOrb
              jOrb2 = jOrb+nFroI
              DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     &          = DPT2_ori(iMO1+iOrb1-1+nOrbI1*(jOrb1-1))
              DPT2(iMO2+jOrb2-1+nOrbI2*(iOrb2-1))
     &          = DPT2_ori(iMO1+jOrb1-1+nOrbI1*(iOrb1-1))
            End Do
          End Do
        End If
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do

      End Subroutine OLagFro0
!
!-----------------------------------------------------------------------
!
      Subroutine OLagFroD(DIA,DI,RDMSA,Trf)

      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NFRO, NISH, NASH, NBAS, NBSQT
      use Constants, only: Zero, One, Two

      implicit none

#include "intent.fh"

      real(kind=wp), intent(_OUT_) :: DIA(*), DI(*)
      real(kind=wp), intent(in) :: RDMSA(*), Trf(*)

      real(kind=wp),allocatable :: WRK1(:),WRK2(:)

      integer(kind=iwp) :: iAOtr, iAOsq, iSym, nFroI, nIshI, nAshI,
     &  nBasI, nCorI

      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')

      iAOtr = 0
      iAOsq = 1
      Do iSym = 1, nSym
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nBasI = nBas(iSym)
        nCorI = nFroI + nIshI

        !! full density matrix
      ! Call SQUARE(WRK1(1+iAOtr),DIA(iAOsq),1,nBasI,nBasI)
      ! !! off-diagonal elements have to be halved
      ! Do Mu = 1, nBasI
      !   Do Nu = 1, nBasI
      !     If (Mu == Nu) Cycle
      !     DIA(iAOsq+Mu-1+nBasI*(Nu-1))
     &!       = Half*DIA(iAOsq+Mu-1+nBasI*(Nu-1))
      !   End Do
      ! End Do

        !! inactive density matrix
        Call DGEMM_('N','T',nBasI,nBasI,nCorI,
     &              Two,CMOPT2,nBasI,CMOPT2,nBasI,
     &              Zero,DI(iAOsq),nBasI)

        !! inactive+active density matrix
        !! Somehow, the above density matrix obtained by calling
        !! Get_D1AO is incorrect... at least, cannot be used.
        ! 1) inactive part
        DIA(1:nBasI**2) = DI(1:nBasI**2)
        ! 2)  RDMSA is defined in CASSCF orbitals, so transform RDMSA to
        !     CASPT2 orbital basis
        Call DGemm_('T','N',nAshI,nAshI,nAshI,
     &              One,Trf(1+nCorI+nBasI*nCorI),nBasI,RDMSA,nAshI,
     &              Zero,WRK2,nAshI)
        Call DGemm_('N','N',nAshI,nAshI,nAshI,
     &              One,WRK2,nAshI,Trf(1+nCorI+nBasI*nCorI),nBasI,
     &              Zero,WRK1,nAshI)
        ! 3) Finally, add the active part
        Call DGemm_('N','N',nBasI,nAshI,nAshI,
     &              One,CMOPT2(1+nBasI*nCorI),nBasI,WRK1,nAshI,
     &              Zero,WRK2,nBasI)
        Call DGemm_('N','T',nBasI,nBasI,nAshI,
     &              One,WRK2,nBasI,CMOPT2(1+nBasI*nCorI),nBasI,
     &              One,DIA,nBasI)

        iAOtr = iAOtr + nBasI*(nBasI+1)/2
        iAOsq = iAOsq + nBasI*nBasI
      End Do

      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)

      Return

      End Subroutine OLagFroD
!
!-----------------------------------------------------------------------
!
      Subroutine OLagFro1(DPT2,OLag)

      use caspt2_global, only: FIFA_all
      use caspt2_module, only: NSYM, NFRO, NISH, NBAS, NDEL
      use Constants, only: Zero, Half
      use definitions, only: wp, iwp

      implicit none

      real(kind=wp), intent(inout) :: DPT2(*), OLag(*)

      integer(kind=iwp) :: iMO, iSym, nOrbI, nFroI, nIshI, nBasI, iOrb,
     &  jOrb
      real(kind=wp) :: Tmp

      iMO  = 1
      DO iSym = 1, nSym
        nOrbI = nBas(iSym)-nDel(iSym)
        nFroI = nFro(iSym)
        If (nOrbI > 0 .and. nFroI > 0) Then
          nIshI = nIsh(iSym)
          nBasI = nBas(iSym)
          !! Make sure that the frozen orbital derivative is zero
          !! (it does not appear in the PT2 energy)
          OLag(1:nOrbI*nFroI) = Zero
          Do iOrb = 1, nFroI
            Do jOrb = nFroI+1, nFroI+nIshI
              Tmp = -Half*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))
     &                    -OLag(iMO+jOrb-1+nOrbI*(iOrb-1)))
     &            /(FIFA_all(iOrb+nBasI*(iOrb-1))
     &             -FIFA_all(jOrb+nBasI*(jOrb-1)))
              DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     &          = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) + Tmp
              DPT2(iMO+jOrb-1+nOrbI*(iOrb-1))
     &          = DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) + Tmp
            End Do
          End Do
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do
!     write(u6,*) 'DPT2 after frozen orbital'
!     call sqprt(dpt2,nbast)

      End Subroutine OLagFro1
!
!-----------------------------------------------------------------------
!
      Subroutine OLagFro2(DPT2,FPT2,ERI,Scr)

      use caspt2_module, only: NSYM, NFRO, NISH, NDEL, NBAS
      use Constants, only: Half
      use definitions, only: wp, iwp

      implicit none

#include "intent.fh"

      real(kind=wp), intent(in) :: DPT2(*)
      real(kind=wp), intent(inout) :: FPT2(*)
      real(kind=wp), intent(_OUT_) :: ERI(*), Scr(*)

      integer(kind=iwp) :: iMO, iSymI, iSymJ, iSymA, iSymB, iSym, nOrbI,
     &  nFroI, nIshI, iOrb, jOrb
      real(kind=wp) :: Scal, Val

!     write(u6,*) 'FPT2 before frozen orbital'
!     call sqprt(fpt2,nbast)
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
            Scal = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,ERI,Scr)
            FPT2(iMO:iMO+nOrbI*nOrbI-1) = FPT2(iMO:iMO+nOrbI*nOrbI-1)
     &        + Scal*ERI(1:nOrbI*nOrbI)
            Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,ERI,Scr)
            FPT2(iMO:iMO+nOrbI*nOrbI-1) = FPT2(iMO:iMO+nOrbI*nOrbI-1)
     &        - Half*Scal*ERI(1:nOrbI*nOrbI)
          End Do
        End Do

        !! Symmetrize FPT2
        Do iOrb = 1, nOrbI
          Do jOrb = 1, iOrb-1
            Val = (FPT2(iMO+iOrb-1+nOrbI*(jOrb-1))
     &            +FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*Half
            FPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
            FPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
          End Do
        End Do
        iMO = iMO + nOrbI*nOrbI
      End Do
!     write(u6,*) 'FPT2 after frozen orbital'
!     call sqprt(fpt2,nbast)

      End Subroutine OLagFro2
!
!-----------------------------------------------------------------------
!
      Subroutine OLagFro3(FIFA,FIMO,WRK1,WRK2)

      use caspt2_global, only: CMOPT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NDEL, NBAS, NBTRI

      implicit none

#include "intent.fh"

      real(kind=wp), intent(inout) :: FIFA(*), FIMO(*)
      real(kind=wp), intent(_OUT_) :: WRK1(*), WRK2(*)

      Character(Len=8) :: Label
      real(kind=wp), allocatable :: WFLT(:)
      integer(kind=iwp) :: IRC, IOPT, ICOMP, ISYLBL, iAO, iAOtr, iCMO,
     &  iMO, iSym, nBasI, nOrbI

      !! Read H_{\mu \nu}
      call mma_allocate(WFLT,NBTRI,Label='WFLT')
      IRC=-1
      IOPT=6
      ICOMP=1
      ISYLBL=1
      Label='OneHam  '
      CALL RDONE(IRC,IOPT,Label,ICOMP,WFLT,ISYLBL)

      !! AO -> MO transformation
      iAO   = 1
      iAOtr = 1
      iCMO  = 1
      iMO   = 1
      DO iSym = 1, nSym
        nBasI = nBas(iSym)
        nOrbI = nBas(iSym)-nDel(iSym)

        !! FIFA
        !! WRK1 = G(D)
        WRK1(1:nBasI*nBasI) = FIFA(iAO:iAO+nBasI*nBasI-1)
        !! WRK1 = H+G(D)
        Call Square(WFLT(iAOtr),WRK2,1,nBasI,nBasI)
        WRK1(1:nBasI*nBasI) = WRK1(1:nBasI*nBasI) + WRK2(1:nBasI*nBasI)
        !! AO -> MO transformation of H+G(D)
        Call OLagTrf(2,iSym,CMOPT2(iCMO),FIFA(iMO),WRK1,WRK2)

        !! FIMO
        !! WRK1 = G(D)
        WRK1(1:nBasI*nBasI) = FIMO(iAO:iAO+nBasI*nBasI-1)
        !! WRK1 = H+G(D)
        Call Square(WFLT(iAOtr),WRK2,1,nBasI,nBasI)
        WRK1(1:nBasI*nBasI) = WRK1(1:nBasI*nBasI) + WRK2(1:nBasI*nBasI)
        !! AO -> MO transformation of H+G(D)
        Call OLagTrf(2,iSym,CMOPT2(iCMO),FIMO(iMO),WRK1,WRK2)

        iAO   = iAO   + nBasI*nBasI
        iAOtr = iAOtr + nBasI*(nBasI+1)/2
        iCMO  = iCMO  + nBasI*nOrbI !?
        iMO   = iMO   + nOrbI*nOrbI
      End Do
!     write(u6,*) 'FIFA'
!     call sqprt(fifa,nbast)
!     write(u6,*) 'FIMO'
!     call sqprt(fimo,nbast)

      call mma_deallocate(WFLT)

      End Subroutine OLagFro3
!
!-----------------------------------------------------------------------
!
      Subroutine OLagFroSq(iSym,Ftr,Fsq)

      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NFRO, NDEL, NBAS, NBAST
      use Constants, only: Zero

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: iSym
      real(kind=wp), intent(in) :: Ftr(*)
      real(kind=wp), intent(_OUT_) :: Fsq(*)

      real(kind=wp),allocatable :: EPS_loc(:)
      integer(kind=iwp) :: nOrbI, nFroI, iOrb, NSEQ, jOrb

      call mma_allocate(EPS_loc,nBasT,Label='EPS_loc')
      Call Get_dArray('RASSCF OrbE',EPS_loc,nBasT)

      nOrbI = nBas(iSym)-nDel(iSym)
      nFroI = nFro(iSym)
      Fsq(1:nOrbI**2) = Zero

      !! Frozen orbital
      Do iOrb = 1, nFroI
        Fsq(iOrb+nOrbI*(iOrb-1)) = EPS_loc(iOrb)
      End Do

      !! Other orbitals
      NSEQ = 0
      Do iOrb = nFroI+1, nOrbI
        Do jOrb = nFroI+1, iOrb
          NSEQ = NSEQ + 1
          Fsq(iOrb+nOrbI*(jOrb-1)) = Ftr(NSEQ)
          Fsq(jOrb+nOrbI*(iOrb-1)) = Ftr(NSEQ)
        End Do
      End Do

      call mma_deallocate(EPS_loc)

      End Subroutine OLagFroSq
!
!-----------------------------------------------------------------------
!
      !! focktwo.f
      SUBROUTINE OLagFro4(iSym0,iSymI,iSymJ,iSymK,iSymL0,
     &                    DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,WRK1)

      USE CHOVEC_IO, only: NVLOC_CHOBATCH
      use Cholesky, only: InfVec, nDimRS
      use ChoCASPT2, only: NUMCHO_PT2, NCHSPC, MXNVC
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use Constants, only: Zero, One, Half
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_module, only: NSYM, NBAS, NBSQT, NBTCHES

      implicit none

#include "warnings.h"
#include "intent.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#endif

      integer(kind=iwp), intent(in) :: iSym0, iSymI, iSymJ, iSymK,
     &  iSymL0
      real(kind=wp), intent(inout) :: DPT2AO(*), DPT2CAO(*)
      real(kind=wp), intent(_OUT_) :: FPT2AO(*), FPT2CAO(*), WRK1(*)

      real(kind=wp), allocatable :: CHSPC(:), WRK2(:)
      integer(kind=iwp) :: ISTLT(8), ISTSQ(8), iSkip(8), ipWRK(8),
     &  nnbstr(8,3), iSym, maxvec, n2, jSym, nB, nB2, nB3, nBasI, nBasJ,
     &  iSymIJ, nBasIJ, nBasK, iSMax, iSymL, nBasL, nBasKL, IBATCH_TOT,
     &  JRED1, JRED2, JSTART, NVECS_RED, ILOC, IRC, JRED, NBATCH, JV1,
     &  IBATCH, JNUM, JV2, JREDC, NUMV, MUSED, ipVecL, iVec, jVref,
     &  lscr, JREDL, JVEC1, iSwap, i, j
      real(kind=wp) :: tmp

      ! INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)

      !! It shoudl be zero, but just in case
      FPT2AO(1:NBSQT) = Zero
      FPT2CAO(1:NBSQT) = Zero

#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        !! To broadcast DPT2AO and DPT2CAO
        If (.not.King()) Then
          DPT2AO(1:NBSQT) = Zero
          DPT2CAO(1:NBSQT) = Zero
        End If
        CALL GADSUM (DPT2AO,NBSQT)
        CALL GADSUM (DPT2CAO,NBSQT)
      End If
#endif

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
      If (iSymI == iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
      If (nBasIJ == 0) Return

      nBasK  = nBas(iSymK)
      iSMax  = iSymK
      If (iSymK == iSymI) iSMax = iSymJ
      iSymL  = 1+iEor(iSymIJ-1,iSymK-1)
      IF (iSymL > iSMax) Return !! should not
      nBasL  = nBas(iSymL0)
      nBasKL = nBasK*nBasL
      IF (iSymK == iSymL0) nBasKL = (nBasK*(nBasK+1))/2
      If (nBasKL == 0) Return

      call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')

      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym) == 0) Return

      ! ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      ! JRED1=iWork(ipnt)
      ! JRED2=iWork(ipnt-1+NumCho_PT2(iSym))
      JRED1=InfVec(1,2,iSym)
      JRED2=InfVec(NumCho_PT2(iSym),2,iSym)
!     write(u6,*) 'jred1,jred2 = ', jred1,jred2

* Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED == 0) Cycle

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
          IBATCH_TOT=IBATCH_TOT+1

          JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
          JV2=JV1+JNUM-1

          JREDC=JRED
* Read a batch of reduced vectors
          CALL CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,
     &                            NUMV,JREDC,MUSED)

          IF(NUMV /= JNUM) THEN
            write(u6,*)' Rats! CHO_VECRD was called, assuming it to'
            write(u6,*)' read JNUM vectors. Instead it returned NUMV'
            write(u6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
            write(u6,*)' Back to the drawing board?'
            CALL QUIT(_RC_INTERNAL_ERROR_)
          END IF
          IF(JREDC /= JRED) THEN
            write(u6,*)' Rats! It was assumed that the Cholesky vectors'
            write(u6,*)' in HALFTRNSF all belonged to a given reduced'
            write(u6,*)' set, but they don''t!'
            write(u6,*)' JRED, JREDC:',JRED,JREDC
            write(u6,*)' Back to the drawing board?'
            write(u6,*)' Let the program continue and see what happens.'
          END IF

          ipVecL = 1
          Do iVec = 1, NUMV
            !! (strange) reduced form -> squared AO vector (mu nu|iVec)
            jVref = 1 !! only for iSwap=1
!           lscr  = nBasI*(nBasI+1)/2
            ! If (l_NDIMRS < 1) Then
            If (size(nDimRS) < 1) Then
              lscr  = NNBSTR(iSym,3)
            Else
              JREDL = INFVEC(iVec,2,iSym)
              ! lscr  = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
              lscr  = nDimRS(iSym,JREDL)
            End If
            JVEC1 = 1
            iSwap = 2
            WRK2(:) = Zero
            Call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,jVref,
     &                      JVEC1,1,1,iSym,JREDC,iSwap,ipWRK,WRK2,
     &                      iSkip)
            ipVecL = ipVecL + lscr
!
!           ----- Fock-like transformations -----
!
            Call FDGTRF_RI(WRK2,DPT2AO ,FPT2AO )
            Call FDGTRF_RI(WRK2,DPT2CAO,FPT2CAO)
          End Do
          JV1=JV1+JNUM
        End Do
      End Do

      call mma_deallocate(CHSPC)
      call mma_deallocate(WRK2)

      !! Have to symmetrize Fock-transformed matrices
      Do i = 1, nBasI
        Do j = 1, i-1
          tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*Half
          FPT2AO(i+nBasI*(j-1)) = Tmp
          FPT2AO(j+nBasI*(i-1)) = Tmp
          tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*Half
          FPT2CAO(i+nBasI*(j-1)) = Tmp
          FPT2CAO(j+nBasI*(i-1)) = Tmp
        End Do
      End Do

#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        CALL GADSUM (FPT2AO,NBSQT)
        CALL GADSUM (FPT2CAO,NBSQT)
      End If
#endif

      Return

      Contains

      Subroutine FDGTRF_RI(ChoVec,DD,FF)

      implicit none

      real(kind=wp), intent(in) :: ChoVec(*), DD(*)
      real(kind=wp), intent(inout) :: FF(*)

      real(kind=wp) :: Scal
      real(kind=wp), external :: ddot_

      !! Coulomb
      Scal = DDot_(nBasI**2,DD,1,ChoVec,1)
      FF(1:nBasI**2) = FF(1:nBasI**2) + Scal*ChoVec(1:nBasI**2)

      !! Exchange
      Call DGEMM_('T','N',nBasI,nBasI,nBasI,
     &            One,ChoVec,nBasI,DD,nBasI,
     &            Zero,WRK1,nBasI)
      Call DGEMM_('T','T',nBasI,nBasI,nBasI,
     &           -Half,ChoVec,nBasI,WRK1,nBasI,
     &            One,FF,nBasI)

      End Subroutine FDGTRF_RI

      End Subroutine OLagFro4
