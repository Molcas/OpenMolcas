!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      SUBROUTINE VVVOX2(nAux,KEEP,iSym,iSymI,iSymJ,iSymK,iSymL,nBasT,   &
     &                  vLag,CMO,WRK,                                   &
     &                  DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,                  &
     &                  FIFA,FIMO)

      use ChoVec_io, only: NVLOC_CHOBATCH
      use Cholesky, only: InfVec
      use caspt2_global, only: LuGAMMA
      use ChoCASPT2, only: numcho_pt2, NCHSPC, MXNVC
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: IFMSCOUP, NSYM, NFROT, NISH, NASH, NSSH, &
     &                         NBAS, JSTATE, iRlxRoot, NBTCHES
      use Constants, only: Zero, One, Half, Two

      implicit none

#include "warnings.h"

      integer(kind=iwp), intent(in) :: nAux(8), KEEP(8), iSym, iSymI,   &
     &  iSymJ, iSymK, iSymL, nBasT
      real(kind=wp), intent(inout) :: vLag(nBasT,nBasT),                &
     &  WRK(nBasT,nBasT), FPT2AO(nBasT**2), FPT2CAO(nBasT**2),          &
     &  FIFA(nBasT**2), FIMO(nBasT**2)
      real(kind=wp), intent(in) :: CMO(nBasT,nBasT), DPT2AO(nBasT**2),  &
     &  DPT2CAO(nBasT**2)

      real(kind=wp), allocatable :: CHSPC(:), HTSPC(:), HTVec(:)
      integer(kind=iwp) :: ISTLT(8), ISTSQ(8), iSkip(8), ipWRK(8), jSym,&
     &  nAuxT, nB, nB2, nB3, nBasI, KEEPI, nBasJ, KEEPJ, iSymIJ, nBasIJ,&
     &  nBasK, KEEPK, iSMax, iSymL_, nBasL, KEEPL, nBasKL, nIshI, nAshI,&
     &  nSshI, nOrbI, IBATCH_TOT, JRED1, JRED2, JRED, JSTART, NVECS_RED,&
     &  ILOC, IRC, NBATCH, JV1, IBATCH, JNUM, JV2, JREDC, NUMV, MUSED,  &
     &  iVec, i, j
      real(kind=wp) :: tmp

      Do jSym = 1, nSym
        iSkip(jSym) = 1
        ipWRK(jSym) = 1
      End Do

      ISTSQ(1)=0
      ISTLT(1)=0
      nAuxT = 0
      Do jSym = 2, nSym
        nB  = nBas(jSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(jSym) = ISTSQ(jSym-1) + nB2
        ISTLT(jSym) = ISTLT(jSym-1) + nB3
        nAuxT = nAuxT + nAux(jSym)
      End Do

      nBasI  = nBas(iSymI)
      KEEPI  = KEEP(iSymI)
      ! nAuxI  = nAux(iSymI)
      nBasJ  = nBas(iSymJ)
      KEEPJ  = KEEP(iSymJ)
      ! nAuxJ  = nAux(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI == iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
      If (nBasIJ == 0) Return

      nBasK  = nBas(iSymK)
      KEEPK  = KEEP(iSymK)
      ! nAuxK  = nAux(iSymK)
      iSMax  = iSymK
      If (iSymK == iSymI) iSMax = iSymJ
      iSymL_ = 1+iEor(iSymIJ-1,iSymK-1)
      IF (iSymL_ > iSMax) Return !! should not
      nBasL  = nBas(iSymL_)
      KEEPL  = KEEP(iSymL_)
      ! nAuxL  = nAux(iSymL_)
      nBasKL = nBasK*nBasL
      IF (iSymK == iSymL) nBasKL = (nBasK*(nBasK+1))/2
      If (nBasKL == 0) Return

      ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
      IF (KEEPI+KEEPJ+KEEPK+KEEPL /= 0) Return
      !! This will not work when the number of the inactive orbital is 0
!     IF (nAuxI+nAuxJ+nAuxK+nAuxL == 0) Return ! frozen orbitals

      jSym = iSymJ
      ! kSym = iSymK
      ! lSym = iSymL
      nIshI = nIsh(iSym)
      ! nIshJ = nIsh(jSym)
      ! nIshK = nIsh(kSym)
      ! nIshL = nIsh(lSym)
      nAshI = nAsh(iSym)
      ! nAshJ = nAsh(jSym)
      ! nAshK = nAsh(kSym)
      ! nAshL = nAsh(lSym)
      nSshI = nSsh(iSym)
      ! nSshJ = nSsh(jSym)
      ! nSshK = nSsh(kSym)
      ! nSshL = nSsh(lSym)
      nOrbI = nIshI+nAshI+nSshI

      call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
      call mma_allocate(HTSPC,NCHSPC,Label='HISPC')
      call mma_allocate(HTVec,nBasT*nBasT,Label='HTVEC')

      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym) == 0) Return

      JRED1=InfVec(1,2,jSym)
      JRED2=InfVec(NumCho_PT2(jSym),2,jSym)

! Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED == 0) Cycle

        ILOC=3
        CALL CHO_X_SETRED(IRC,ILOC,JRED)
! For a reduced set, the structure is known, including
! the mapping between reduced index and basis set pairs.
! The reduced set is divided into suitable batches.
! First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.
        ! JEND=JSTART+NVECS_RED-1

! Determine batch length for this reduced set.
! Make sure to use the same formula as in the creation of disk
! address tables, etc, above:
        NBATCH=1+(NVECS_RED-1)/MXNVC

! Loop over IBATCH
        JV1=JSTART
        DO IBATCH=1,NBATCH
!         Write(u6,*) 'ibatch,nbatch = ', ibatch,nbatch
          IBATCH_TOT=IBATCH_TOT+1

          JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
          JV2=JV1+JNUM-1
!
!         ----- Construct orbital Lagrangian -----
!
          !! CHSPC       :: (mu nu|P)
          !! HTSPC       :: ( q nu|P)
          !! Bra and Ket :: (ia|P) = T_{ij}^{ab}*(jb|P)
          !! In 1), HT_{i mu,P} = C_{mu a}*(ia|P)
          !!        (ia|P) read from disk
          !! In 3), L_{i mu} = HT_{i nu,P} * (mu nu|P)
          !! Then, L_{pi} = C_{mu p} * L_{i mu}^T
          !! Transpose is done in VVVO_Drv2,
          !! and C_{mu p} is in OLagVVVO

          !! 1) Half back-transformation of Bra and Ket density
          !! Read the 3c-2e pseudo-density (in MO), and half transform
          CALL VVVOTRA_RI(CMO,CHSPC,SIZE(CHSPC),HTSPC,                  &
     &                    JNUM,IBATCH_TOT,IBATCH_TOT,nOrbI)

          !! 2) read AO Cholesky vectors,
          !!    then, (strange) reduced form -> squared AO (mu nu|iVec)
          JREDC=JRED
! Read a batch of reduced vectors
          CALL CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,                     &
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

          !! (strange) reduced form -> squared AO (mu nu|iVec)
          !! is it possible to avoid this transformation?
      ! choptr.fh
          Call R2FIP(CHSPC,SIZE(CHSPC),WRK,ipWRK,NUMV,                  &
     &               nBasT,iSym,iSkip,irc,JREDC)
!
!           ----- Fock-like transformations (if needed) -----
!
          If (nFroT == 0) Then
            Do iVec = 1, NUMV
              Call FDGTRF(CHSPC(1+nBasT**2*(iVec-1)),                   &
     &                    DPT2AO,FPT2AO)
              Call FDGTRF(CHSPC(1+nBasT**2*(iVec-1)),                   &
     &                    DPT2CAO,FPT2CAO)
            End Do
          End If

          !! 3) Contract with Cholesky vectors
          Call DGemm_('N','T',nOrbI,nBasI,nBasI*JNUM,                   &
     &                Two,HTSPC,nOrbI,CHSPC,nBasI,                      &
     &                One,vLag,nBasI)

          !! 4) Construct the 3c-2e pseudo-density in AO
          !! D_{p nu} -> D_{mu nu}
          !! i.e., construct B_PT2, used in ALASKA
          Call DGemm_('N','N',nBasI,nBasI*JNUM,nOrbI,                   &
     &                One,CMO(1,1),nBasI,HTSPC,nOrbI,                   &
     &                Zero,CHSPC,nBasI)

          !! 5) Save the 3c-2e pseudo-density in the disk
          !! it may be replaced with ddafile
          Do iVec = 1, NUMV
            If (IFMSCOUP .and. jState /= 1) Then
              Read (LuGamma,Rec=iVec+JV1-1) HTVec(1:nBasI**2)
              Call DaXpY_(nBasI**2,One,                                 &
     &                    CHSPC(1+nBasI**2*(iVec-1)),1,                 &
     &                    HTVec,1)
              Write (LuGamma,Rec=iVec+JV1-1) HTVec(1:nBasI**2)
            Else
             if (jState == iRlxRoot .or. IFMSCOUP) then
              Write (LuGamma,Rec=iVec+JV1-1)                            &
     &        CHSPC(1+nBasI**2*(iVec-1):nBasI**2*iVec)
             end if
            End If
          End Do

          JV1=JV1+JNUM
        End Do
      End Do

      call mma_deallocate(CHSPC)
      call mma_deallocate(HTSPC)
      call mma_deallocate(HTVec)

      !! Have to (?) symmetrize Fock-transformed matrices
      If (nFroT == 0) Then
        Do i = 1, nBasI
          Do j = 1, i-1
            tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*Half
            FPT2AO(i+nBasI*(j-1)) = Tmp
            FPT2AO(j+nBasI*(i-1)) = Tmp
            tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*Half
            FPT2CAO(i+nBasI*(j-1)) = Tmp
            FPT2CAO(j+nBasI*(i-1)) = Tmp
            If (nFroT /= 0) Then
              tmp = (FIFA(i+nBasI*(j-1))+FIFA(j+nBasI*(i-1)))*Half
              FIFA(i+nBasI*(j-1)) = Tmp
              FIFA(j+nBasI*(i-1)) = Tmp
              tmp = (FIMO(i+nBasI*(j-1))+FIMO(j+nBasI*(i-1)))*Half
              FIMO(i+nBasI*(j-1)) = Tmp
              FIMO(j+nBasI*(i-1)) = Tmp
            End If
          End Do
        End Do
      End If

      Return

      Contains

      Subroutine FDGTRF(ChoVec,DD,FF)

      implicit none

      real(kind=wp), intent(in) :: ChoVec(nBasI**2), DD(nBasI**2)
      real(kind=wp), intent(inout) :: FF(nBasI**2)

      real(kind=wp) :: Scal
      real(kind=wp), external :: ddot_

      !! Coulomb
      Scal = DDot_(nBasI**2,DD,1,ChoVec,1)
      Call DaXpY_(nBasI**2,Scal,ChoVec,1,FF,1)

      !! Exchange
      Call DGEMM_('T','N',nBasI,nBasI,nBasI,                            &
     &            One,ChoVec,nBasI,DD,nBasI,                            &
     &            Zero,WRK,nBasI)
      Call DGEMM_('T','T',nBasI,nBasI,nBasI,                            &
     &           -Half,ChoVec,nBasI,WRK,nBasI,                          &
     &            One,FF,nBasI)

      End Subroutine FDGTRF

      Subroutine VVVOTRA_RI(CMO,CHSPC_,NCHSPC,HTSPC_,NVEC,IBSTA,IBEND,  &
     &                      nOrbI)

      implicit none

      integer(kind=iwp), intent(in) :: NCHSPC, NVEC, IBSTA, IBEND, nOrbI
      real(kind=wp), intent(in) :: CMO(nBasI,nOrbI)
      !! CHSPC is used as a temporary array
      real(kind=wp), intent(inout) :: CHSPC_(NCHSPC),                   &
     &                                HTSPC_(nOrbI,nBasT,NVEC)

      integer(kind=iwp), parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) :: IPQ, jVec, nBra

      !! BraAI
      Call Cholesky_Vectors(2,Inactive,Active,JSYM,CHSPC_,NCHSPC,       &
     &                      nBra,IBSTA,IBEND)
      IPQ = nAshI*nIshI
      Do jVec = 1, NVEC
        ! a. AI -> mu I
        Call DGemm_('T','T',nIshI,nBasI,nAshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,                   &
     &                  CMO(1,1+nIshI),nBasI,                           &
     &              Zero,HTSPC_(1,1,jVec),nOrbI)
        ! a. AI -> A mu
        Call DGemm_('N','T',nAshI,nBasI,nIshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,                   &
     &                  CMO(1,1),nBasI,                                 &
     &              Zero,HTSPC_(1+nIshI,1,jVec),nOrbI)
      End Do

      !! BraSI
      Call Cholesky_Vectors(2,Inactive,Virtual,JSYM,CHSPC_,NCHSPC,      &
     &                      nBra,IBSTA,IBEND)
      IPQ = nIshI*nSshI
      Do jVec = 1, NVEC
        ! b. SI -> mu I
        Call DGemm_('T','T',nIshI,nBasI,nSshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,                   &
     &                  CMO(1,1+nIshI+nAshI),nBasI,                     &
     &              One,HTSPC_(1,1,jVec),nOrbI)
        ! a. SI -> S mu
        Call DGemm_('N','T',nSshI,nBasI,nIshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,                   &
     &                  CMO(1,1),nBasI,                                 &
     &              Zero,HTSPC_(1+nIshI+nAshI,1,jVec),nOrbI)
      End Do

      !! BraSA
      Call Cholesky_Vectors(2,Active,Virtual,JSYM,CHSPC_,NCHSPC,        &
     &                      nBra,IBSTA,IBEND)
      IPQ = nAshI*nSshI
      Do jVec = 1, NVEC
        ! d. SA -> mu A
        Call DGemm_('T','T',nAshI,nBasI,nSshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,                   &
     &                  CMO(1,1+nIshI+nAshI),nBasI,                     &
     &              One,HTSPC_(1+nIshI,1,jVec),nOrbI)
        ! b. SA -> S mu
        Call DGemm_('N','T',nSshI,nBasI,nAshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nSshI,                   &
     &                  CMO(1,1+nIshI),nBasI,                           &
     &              One,HTSPC_(1+nIshI+nAshI,1,jVec),nOrbI)
      End Do

      !! BraAA
      Call Cholesky_Vectors(2,Active,Active,JSYM,CHSPC_,NCHSPC,         &
     &                      nBra,IBSTA,IBEND)
      IPQ = nAshI*nAshI
      Do jVec = 1, NVEC
        ! b. AA -> mu A
        Call DGemm_('T','T',nAshI,nBasI,nAshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,                   &
     &                  CMO(1,1+nIshI),nBasI,                           &
     &              One,HTSPC_(1+nIshI,1,jVec),nOrbI)
        ! c. AA -> A mu
        Call DGemm_('N','T',nAshI,nBasI,nAshI,                          &
     &              One,CHSPC_(1+IPQ*(jVec-1)),nAshI,                   &
     &                  CMO(1,1+nIshI),nBasI,                           &
     &              One,HTSPC_(1+nIshI,1,jVec),nOrbI)
      End Do

      End Subroutine VVVOTRA_RI

      End Subroutine VVVOX2
