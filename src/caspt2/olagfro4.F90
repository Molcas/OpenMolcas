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

subroutine OLagFro4(NBSQT,iSym0,iSymI,iSymJ,iSymK,iSymL0,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,WRK1)

use Symmetry_Info, only: Mul
use CHOVEC_IO, only: NVLOC_CHOBATCH
use Cholesky, only: InfVec, nDimRS, nnBstR
use ChoCASPT2, only: MXNVC, NCHSPC, NUMCHO_PT2
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use caspt2_module, only: NBAS, NBTCHES, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NBSQT, iSym0, iSymI, iSymJ, iSymK, iSymL0
real(kind=wp), intent(inout) :: DPT2AO(NBSQT), DPT2CAO(NBSQT)
real(kind=wp), intent(out) :: FPT2AO(NBSQT), FPT2CAO(NBSQT), WRK1(NBSQT)
integer(kind=iwp) :: i, IBATCH, IBATCH_TOT, ILOC, ipVecL, ipWRK(8), IRC, iSkip(8), iSMax, ISTLT(8), ISTSQ(8), iSwap, iSym, iSymIJ, &
                     iSymL, iVec, j, JNUM, JRED, JRED1, JRED2, JREDC, JREDL, JSTART, jSym, JV1, JV2, JVEC1, jVref, lscr, MUSED, &
                     nB, nB2, nB3, nBasI, nBasIJ, nBasJ, nBasK, nBasKL, nBasL, NBATCH, NUMV, NVECS_RED
real(kind=wp) :: tmp
real(kind=wp), allocatable :: CHSPC(:), WRK2(:)
#include "warnings.h"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#endif

!! It shoudl be zero, but just in case
FPT2AO(1:NBSQT) = Zero
FPT2CAO(1:NBSQT) = Zero

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  !! To broadcast DPT2AO and DPT2CAO
  if (.not. King()) then
    DPT2AO(1:NBSQT) = Zero
    DPT2CAO(1:NBSQT) = Zero
  end if
  call GADGOP(DPT2AO,NBSQT,'+')
  call GADGOP(DPT2CAO,NBSQT,'+')
end if
#endif

iSym = iSym0

ISTSQ(1) = 0
ISTLT(1) = 0
do jSym=2,nSym
  nB = nBas(jSym-1)
  nB2 = nB*nB
  nB3 = (nB2+nB)/2
  ISTSQ(jSym) = ISTSQ(jSym-1)+nB2
  ISTLT(jSym) = ISTLT(jSym-1)+nB3
end do
do jSym=1,nSym
  iSkip(jSym) = 1
  ipWRK(jSym) = 1
end do

nBasI = nBas(iSymI)
nBasJ = nBas(iSymJ)
iSymIJ = Mul(iSymI,iSymJ)
nBasIJ = nBasI*nBasJ
if (iSymI == iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
if (nBasIJ == 0) return

nBasK = nBas(iSymK)
iSMax = iSymK
if (iSymK == iSymI) iSMax = iSymJ
iSymL = Mul(iSymIJ,iSymK)
if (iSymL > iSMax) return !! should not
nBasL = nBas(iSymL0)
nBasKL = nBasK*nBasL
if (iSymK == iSymL0) nBasKL = (nBasK*(nBasK+1))/2
if (nBasKL == 0) return

call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
call mma_allocate(WRK2,NBSQT,Label='WRK2')

IBATCH_TOT = NBTCHES(iSym)

if (NUMCHO_PT2(iSym) == 0) return

!ipnt = ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
!JRED1 = iWork(ipnt)
!JRED2 = iWork(ipnt-1+NumCho_PT2(iSym))
JRED1 = InfVec(1,2,iSym)
JRED2 = InfVec(NumCho_PT2(iSym),2,iSym)
!write(u6,*) 'jred1,jred2 = ',jred1,jred2

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
    IBATCH_TOT = IBATCH_TOT+1

    JNUM = NVLOC_CHOBATCH(IBATCH_TOT)
    JV2 = JV1+JNUM-1

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

    ipVecL = 1
    do iVec=1,NUMV
      !! (strange) reduced form -> squared AO vector (mu nu|iVec)
      jVref = 1 !! only for iSwap=1
      !lscr = nBasI*(nBasI+1)/2
      !if (l_NDIMRS < 1) then
      if (size(nDimRS) < 1) then
        lscr = NNBSTR(iSym,3)
      else
        JREDL = INFVEC(iVec,2,iSym)
        !lscr = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
        lscr = nDimRS(iSym,JREDL)
      end if
      JVEC1 = 1
      iSwap = 2
      WRK2(:) = Zero
      call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,jVref,JVEC1,1,1,iSym,JREDC,iSwap,ipWRK,WRK2,iSkip)
      ipVecL = ipVecL+lscr

      ! ----- Fock-like transformations -----

      call FDGTRF_RI(WRK2,DPT2AO,FPT2AO)
      call FDGTRF_RI(WRK2,DPT2CAO,FPT2CAO)
    end do
    JV1 = JV1+JNUM
  end do
end do

call mma_deallocate(CHSPC)
call mma_deallocate(WRK2)

!! Have to symmetrize Fock-transformed matrices
do i=1,nBasI
  do j=1,i-1
    tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*Half
    FPT2AO(i+nBasI*(j-1)) = Tmp
    FPT2AO(j+nBasI*(i-1)) = Tmp
    tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*Half
    FPT2CAO(i+nBasI*(j-1)) = Tmp
    FPT2CAO(j+nBasI*(i-1)) = Tmp
  end do
end do

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GADGOP(FPT2AO,NBSQT,'+')
  call GADGOP(FPT2CAO,NBSQT,'+')
end if
#endif

return

contains

subroutine FDGTRF_RI(ChoVec,DD,FF)

  real(kind=wp), intent(in) :: ChoVec(nBasI**2), DD(nBasI**2)
  real(kind=wp), intent(inout) :: FF(nBasI**2)
  real(kind=wp) :: Scal
  real(kind=wp), external :: ddot_

  !! Coulomb
  Scal = DDot_(nBasI**2,DD,1,ChoVec,1)
  FF(1:nBasI**2) = FF(1:nBasI**2)+Scal*ChoVec(1:nBasI**2)

  !! Exchange
  call DGEMM_('T','N',nBasI,nBasI,nBasI,One,ChoVec,nBasI,DD,nBasI,Zero,WRK1,nBasI)
  call DGEMM_('T','T',nBasI,nBasI,nBasI,-Half,ChoVec,nBasI,WRK1,nBasI,One,FF,nBasI)

end subroutine FDGTRF_RI

end subroutine OLagFro4
