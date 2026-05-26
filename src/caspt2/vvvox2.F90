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

subroutine VVVOX2(KEEP,iSym,iSymI,iSymJ,iSymK,iSymL,nBasT,vLag,CMO,WRK,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,FIFA,FIMO)

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use ChoVec_io, only: NVLOC_CHOBATCH
use Cholesky, only: InfVec
use caspt2_global, only: LuGAMMA
use ChoCASPT2, only: MXNVC, NCHSPC, numcho_pt2
use caspt2_module, only: IFMSCOUP, iRlxRoot, JSTATE, NASH, NBAS, NBTCHES, NFROT, NISH, NSSH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: KEEP(8), iSym, iSymI, iSymJ, iSymK, iSymL, nBasT
real(kind=wp), intent(inout) :: vLag(nBasT,nBasT), WRK(nBasT,nBasT), FPT2AO(nBasT**2), FPT2CAO(nBasT**2), FIFA(nBasT**2), &
                                FIMO(nBasT**2)
real(kind=wp), intent(in) :: CMO(nBasT,nBasT), DPT2AO(nBasT**2), DPT2CAO(nBasT**2)
integer(kind=iwp) :: i, IBATCH, IBATCH_TOT, ILOC, ipWRK(8), IRC, iSkip(8), iSMax, ISTLT(8), ISTSQ(8), iSymIJ, iSymL_, iVec, j, &
                     JNUM, JRED, JRED1, JRED2, JREDC, JSTART, jSym, JV1, JV2, KEEPI, KEEPJ, KEEPK, KEEPL, MUSED, nAshI, nB, nB2, &
                     nB3, nBasI, nBasIJ, nBasJ, nBasK, nBasKL, nBasL, NBATCH, nIshI, nOrbI, nSshI, NUMV, NVECS_RED
real(kind=wp) :: tmp
real(kind=wp), allocatable :: CHSPC(:), HTSPC(:), HTVec(:)
#include "warnings.h"

iSkip(1:nSym) = 1
ipWRK(1:nSym) = 1

ISTSQ(1) = 0
ISTLT(1) = 0
do jSym=2,nSym
  nB = nBas(jSym-1)
  nB2 = nB**2
  nB3 = nTri_Elem(nB2)
  ISTSQ(jSym) = ISTSQ(jSym-1)+nB2
  ISTLT(jSym) = ISTLT(jSym-1)+nB3
end do

nBasI = nBas(iSymI)
KEEPI = KEEP(iSymI)
!nAuxI = nAux(iSymI)
nBasJ = nBas(iSymJ)
KEEPJ = KEEP(iSymJ)
!nAuxJ = nAux(iSymJ)
iSymIJ = Mul(iSymI,iSymJ)
nBasIJ = nBasI*nBasJ
if (iSymI == iSymJ) nBasIJ = nTri_Elem(nBasI)
if (nBasIJ == 0) return

nBasK = nBas(iSymK)
KEEPK = KEEP(iSymK)
!nAuxK = nAux(iSymK)
iSMax = iSymK
if (iSymK == iSymI) iSMax = iSymJ
iSymL_ = Mul(iSymIJ,iSymK)
if (iSymL_ > iSMax) return !! should not
nBasL = nBas(iSymL_)
KEEPL = KEEP(iSymL_)
!nAuxL = nAux(iSymL_)
nBasKL = nBasK*nBasL
if (iSymK == iSymL) nBasKL = nTri_Elem(nBasK)
if (nBasKL == 0) return

! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
if (KEEPI+KEEPJ+KEEPK+KEEPL /= 0) return
!! This will not work when the number of the inactive orbital is 0
!if (nAuxI+nAuxJ+nAuxK+nAuxL == 0) return ! frozen orbitals

jSym = iSymJ
!kSym = iSymK
!lSym = iSymL
nIshI = nIsh(iSym)
!nIshJ = nIsh(jSym)
!nIshK = nIsh(kSym)
!nIshL = nIsh(lSym)
nAshI = nAsh(iSym)
!nAshJ = nAsh(jSym)
!nAshK = nAsh(kSym)
!nAshL = nAsh(lSym)
nSshI = nSsh(iSym)
!nSshJ = nSsh(jSym)
!nSshK = nSsh(kSym)
!nSshL = nSsh(lSym)
nOrbI = nIshI+nAshI+nSshI

call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
call mma_allocate(HTSPC,NCHSPC,Label='HISPC')
call mma_allocate(HTVec,nBasT**2,Label='HTVEC')

IBATCH_TOT = NBTCHES(iSym)

if (NUMCHO_PT2(iSym) == 0) return

JRED1 = InfVec(1,2,jSym)
JRED2 = InfVec(NumCho_PT2(jSym),2,jSym)

! Loop over JRED
do JRED=JRED1,JRED2

  call Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
  if (NVECS_RED == 0) cycle

  ILOC = 3
  call CHO_X_SETRED(IRC,ILOC,JRED)
  ! For a reduced set, the structure is known, including
  ! the mapping between reduced index and basis set pairs.
  ! The reduced set is divided into suitable batches.
  ! First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.
  !JEND = JSTART+NVECS_RED-1

  ! Determine batch length for this reduced set.
  ! Make sure to use the same formula as in the creation of disk
  ! address tables, etc, above:
  NBATCH = 1+(NVECS_RED-1)/MXNVC

  ! Loop over IBATCH
  JV1 = JSTART
  do IBATCH=1,NBATCH
    !write(u6,*) 'ibatch,nbatch = ',ibatch,nbatch
    IBATCH_TOT = IBATCH_TOT+1

    JNUM = NVLOC_CHOBATCH(IBATCH_TOT)
    JV2 = JV1+JNUM-1

    ! ----- Construct orbital Lagrangian -----

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
    call VVVOTRA_RI(CMO,CHSPC,size(CHSPC),HTSPC,JNUM,IBATCH_TOT,IBATCH_TOT,nOrbI)

    !! 2) read AO Cholesky vectors,
    !!    then, (strange) reduced form -> squared AO (mu nu|iVec)
    JREDC = JRED
    ! Read a batch of reduced vectors
    call CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,NUMV,JREDC,MUSED)
    if (NUMV /= JNUM) then
      write(u6,*) ' Rats! CHO_VECRD was called, assuming it to'
      write(u6,*) ' read JNUM vectors. Instead it returned NUMV'
      write(u6,*) ' vectors: JNUM, NUMV=',JNUM,NUMV
      write(u6,*) ' Back to the drawing board?'
      call QUIT(_RC_INTERNAL_ERROR_)
    end if
    if (JREDC /= JRED) then
      write(u6,*) ' Rats! It was assumed that the Cholesky vectors'
      write(u6,*) ' in HALFTRNSF all belonged to a given reduced'
      write(u6,*) ' set, but they don''t!'
      write(u6,*) ' JRED, JREDC:',JRED,JREDC
      write(u6,*) ' Back to the drawing board?'
      write(u6,*) ' Let the program continue and see what happens.'
    end if

    !! (strange) reduced form -> squared AO (mu nu|iVec)
    !! is it possible to avoid this transformation?
    ! choptr.fh
    call R2FIP(CHSPC,size(CHSPC),WRK,ipWRK,NUMV,nBasT,iSym,iSkip,irc,JREDC)

    ! ----- Fock-like transformations (if needed) -----

    if (nFroT == 0) then
      do iVec=1,NUMV
        call FDGTRF(CHSPC(1+nBasT**2*(iVec-1)),DPT2AO,FPT2AO)
        call FDGTRF(CHSPC(1+nBasT**2*(iVec-1)),DPT2CAO,FPT2CAO)
      end do
    end if

    !! 3) Contract with Cholesky vectors
    call DGemm_('N','T',nOrbI,nBasI,nBasI*JNUM,Two,HTSPC,nOrbI,CHSPC,nBasI,One,vLag,nBasI)

    !! 4) Construct the 3c-2e pseudo-density in AO
    !! D_{p nu} -> D_{mu nu}
    !! i.e., construct B_PT2, used in ALASKA
    call DGemm_('N','N',nBasI,nBasI*JNUM,nOrbI,One,CMO(1,1),nBasI,HTSPC,nOrbI,Zero,CHSPC,nBasI)

    !! 5) Save the 3c-2e pseudo-density in the disk
    !! it may be replaced with ddafile
    do iVec=1,NUMV
      if (IFMSCOUP .and. (jState /= 1)) then
        read(LuGamma,rec=iVec+JV1-1) HTVec(1:nBasI**2)
        HTVec(1:nBasI**2) = HTVec(1:nBasI**2)+CHSPC(nBasI**2*(iVec-1)+1:nBasI**2*iVec)
        write(LuGamma,rec=iVec+JV1-1) HTVec(1:nBasI**2)
      else if ((jState == iRlxRoot) .or. IFMSCOUP) then
        write(LuGamma,rec=iVec+JV1-1) CHSPC(1+nBasI**2*(iVec-1):nBasI**2*iVec)
      end if
    end do

    JV1 = JV1+JNUM
  end do
end do

call mma_deallocate(CHSPC)
call mma_deallocate(HTSPC)
call mma_deallocate(HTVec)

!! Have to (?) symmetrize Fock-transformed matrices
if (nFroT == 0) then
  do i=1,nBasI
    do j=1,i-1
      tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*Half
      FPT2AO(i+nBasI*(j-1)) = Tmp
      FPT2AO(j+nBasI*(i-1)) = Tmp
      tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*Half
      FPT2CAO(i+nBasI*(j-1)) = Tmp
      FPT2CAO(j+nBasI*(i-1)) = Tmp
      if (nFroT /= 0) then
        tmp = (FIFA(i+nBasI*(j-1))+FIFA(j+nBasI*(i-1)))*Half
        FIFA(i+nBasI*(j-1)) = Tmp
        FIFA(j+nBasI*(i-1)) = Tmp
        tmp = (FIMO(i+nBasI*(j-1))+FIMO(j+nBasI*(i-1)))*Half
        FIMO(i+nBasI*(j-1)) = Tmp
        FIMO(j+nBasI*(i-1)) = Tmp
      end if
    end do
  end do
end if

return

contains

subroutine FDGTRF(ChoVec,DD,FF)

  real(kind=wp), intent(in) :: ChoVec(nBasI**2), DD(nBasI**2)
  real(kind=wp), intent(inout) :: FF(nBasI**2)
  real(kind=wp), external :: ddot_

  !! Coulomb
  FF(:) = FF(:)+DDot_(nBasI**2,DD,1,ChoVec,1)*ChoVec(:)

  !! Exchange
  call DGEMM_('T','N',nBasI,nBasI,nBasI,One,ChoVec,nBasI,DD,nBasI,Zero,WRK,nBasI)
  call DGEMM_('T','T',nBasI,nBasI,nBasI,-Half,ChoVec,nBasI,WRK,nBasI,One,FF,nBasI)

end subroutine FDGTRF

subroutine VVVOTRA_RI(CMO,CHSPC_,NCHSPC,HTSPC_,NVEC,IBSTA,IBEND,nOrbI)
  !! CHSPC is used as a temporary array

  integer(kind=iwp), intent(in) :: NCHSPC, NVEC, IBSTA, IBEND, nOrbI
  real(kind=wp), intent(in) :: CMO(nBasI,nOrbI)
  real(kind=wp), intent(inout) :: CHSPC_(NCHSPC), HTSPC_(nOrbI,nBasT,NVEC)
  integer(kind=iwp) :: IPQ, jVec, nBra
  integer(kind=iwp), parameter :: Inactive = 1, Active = 2, Virtual = 3

  !! BraAI
  call Cholesky_Vectors(2,Inactive,Active,JSYM,CHSPC_,NCHSPC,nBra,IBSTA,IBEND)
  IPQ = nAshI*nIshI
  do jVec=1,NVEC
    ! a. AI -> mu I
    call DGemm_('T','T',nIshI,nBasI,nAshI,One,CHSPC_(1+IPQ*(jVec-1)),nAshI,CMO(1,1+nIshI),nBasI,Zero,HTSPC_(1,1,jVec),nOrbI)
    ! a. AI -> A mu
    call DGemm_('N','T',nAshI,nBasI,nIshI,One,CHSPC_(1+IPQ*(jVec-1)),nAshI,CMO(1,1),nBasI,Zero,HTSPC_(1+nIshI,1,jVec),nOrbI)
  end do

  !! BraSI
  call Cholesky_Vectors(2,Inactive,Virtual,JSYM,CHSPC_,NCHSPC,nBra,IBSTA,IBEND)
  IPQ = nIshI*nSshI
  do jVec=1,NVEC
    ! b. SI -> mu I
    call DGemm_('T','T',nIshI,nBasI,nSshI,One,CHSPC_(1+IPQ*(jVec-1)),nSshI,CMO(1,1+nIshI+nAshI),nBasI,One,HTSPC_(1,1,jVec),nOrbI)
    ! a. SI -> S mu
    call DGemm_('N','T',nSshI,nBasI,nIshI,One,CHSPC_(1+IPQ*(jVec-1)),nSshI,CMO(1,1),nBasI,Zero,HTSPC_(1+nIshI+nAshI,1,jVec),nOrbI)
  end do

  !! BraSA
  call Cholesky_Vectors(2,Active,Virtual,JSYM,CHSPC_,NCHSPC,nBra,IBSTA,IBEND)
  IPQ = nAshI*nSshI
  do jVec=1,NVEC
    ! d. SA -> mu A
    call DGemm_('T','T',nAshI,nBasI,nSshI,One,CHSPC_(1+IPQ*(jVec-1)),nSshI,CMO(1,1+nIshI+nAshI),nBasI,One,HTSPC_(1+nIshI,1,jVec), &
                nOrbI)
    ! b. SA -> S mu
    call DGemm_('N','T',nSshI,nBasI,nAshI,One,CHSPC_(1+IPQ*(jVec-1)),nSshI,CMO(1,1+nIshI),nBasI,One,HTSPC_(1+nIshI+nAshI,1,jVec), &
                nOrbI)
  end do

  !! BraAA
  call Cholesky_Vectors(2,Active,Active,JSYM,CHSPC_,NCHSPC,nBra,IBSTA,IBEND)
  IPQ = nAshI*nAshI
  do jVec=1,NVEC
    ! b. AA -> mu A
    call DGemm_('T','T',nAshI,nBasI,nAshI,One,CHSPC_(1+IPQ*(jVec-1)),nAshI,CMO(1,1+nIshI),nBasI,One,HTSPC_(1+nIshI,1,jVec),nOrbI)
    ! c. AA -> A mu
    call DGemm_('N','T',nAshI,nBasI,nAshI,One,CHSPC_(1+IPQ*(jVec-1)),nAshI,CMO(1,1+nIshI),nBasI,One,HTSPC_(1+nIshI,1,jVec),nOrbI)
  end do

end subroutine VVVOTRA_RI

end subroutine VVVOX2
