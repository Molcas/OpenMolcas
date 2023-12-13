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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!               1998, Roland Lindh                                     *
!***********************************************************************

subroutine TraDrv(IPR,lSquare,iSym,jSym,kSym,lSym,iBas,jBas,kBas,lBas,iOrb,jOrb,kOrb,lOrb,iFro,jFro,kFro,lFro,iIsh,jIsh,kIsh,lIsh, &
                  iAsh,jAsh,kAsh,lAsh,ij_Bas_pairs,kl_Bas_pairs,ij_Orb_pairs,kl_Orb_pairs,off_PUVX,off_sqMat,off_ltMat,mxSym,CMO, &
                  PUVX,D1I,FI,D1A,FA,ExFac)
!***********************************************************************
!                                                                      *
!     control subsection for:                                          *
!     - transformation of ERIs from AO to MO basis                     *
!     - Fock matrix generation                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!     Modified to only need unique symmetry blocks, R. Lindh, March '98*
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IPR, iSym, jSym, kSym, lSym, iBas, jBas, kBas, lBas, iOrb, jOrb, kOrb, lOrb, iFro, jFro, kFro, &
                                 lFro, iIsh, jIsh, kIsh, lIsh, iAsh, jAsh, kAsh, lAsh, ij_Bas_pairs, kl_Bas_pairs, ij_Orb_pairs, &
                                 kl_Orb_pairs, mxSym, off_PUVX(mxSym,mxSym,mxSym), off_sqMat(*), off_ltMat(*)
logical(kind=iwp), intent(in) :: lSquare
real(kind=wp), intent(in) :: CMO(*), D1I(*), D1A(*), ExFac
real(kind=wp), intent(inout) :: PUVX(*), FI(*), FA(*)
#include "timers.fh"
integer(kind=iwp) :: case1, case2, i, i1, i2, iiiOff, iiOff, ij_pair, iOff, iOpt, iRc, j, jjjOff, jjOff, jMax, kkkOff, kkOff, &
                     kl_pair, klBas, lllOff, llOff, nBuf2, nBuf3, nInBuf, nOff, nPairs, nPQVX, nScrt1, nTURS
real(kind=wp) :: dum1, dum2, dum3
logical(kind=iwp) :: Process_Twice
real(kind=wp), allocatable :: Buf2(:), Buf3(:)
real(kind=wp), allocatable, target :: InBuf(:), PQVX(:), Scrt1(:), TURS(:)
real(kind=wp), pointer :: Buf9(:) => null(), PQRS(:) => null(), PQRS_(:) => null()

! generate offsets
iiOff = off_sqMat(iSym)+iFro*iBas+1
jjOff = off_sqMat(jSym)+jFro*jBas+1
kkOff = off_sqMat(kSym)+kFro*kBas+1
llOff = off_sqMat(lSym)+lFro*lBas+1
iiiOff = iiOff+iIsh*iBas
jjjOff = jjOff+jIsh*jBas
kkkOff = kkOff+kIsh*kBas
lllOff = llOff+lIsh*lBas

! Generate logical flag for unique symmetry blocks

Process_Twice = (iTri(iSym,jSym) /= iTri(kSym,lSym)) .and. (.not. lSquare)

! find cases
case1 = 4
if (iSym == jSym) case1 = case1-2
if (iSym == kSym) case1 = case1-1
case2 = iAsh*jAsh*kAsh*lAsh

! quit, if this integral block is not used
if ((case1 == 4) .and. (case2 == 0)) return

! allocate memory
nScrt1 = 0
if (iSym == jSym) nScrt1 = iBas*jBas
if (kSym == lSym) nScrt1 = max(nScrt1,kBas*lBas)
nBuf2 = max(jBas*iAsh,kBas*lAsh,iBas*jAsh)
nBuf3 = max(jOrb*iAsh,kAsh*lAsh,iOrb*jAsh)

! PQVX and TURS are the halftransformed integrals

nPQVX = ij_Bas_pairs*kl_Orb_pairs
if (Process_Twice) then
  nTURS = ij_Orb_pairs*kl_Bas_pairs
  nBuf2 = max(nBuf2,ij_Orb_pairs,lBas*kAsh)
  nBuf3 = max(nBuf3,kAsh*lOrb,kOrb*lAsh)
else
  nTURS = 0
end if

if (nScrt1 /= 0) call mma_allocate(Scrt1,nScrt1,Label='Scrt1')
call mma_allocate(Buf2,nBuf2,Label='Buf2')
call mma_allocate(Buf3,nBuf3,Label='Buf3')

call mma_allocate(PQVX,nPQVX,Label='PQVX')
PQVX(:) = Zero
if (nTURS > 0) then
  call mma_allocate(TURS,nTURS,Label='TURS')
  TURS(:) = Zero
end if

call mma_maxDBLE(nInBuf)
nInBuf = min(nInBuf,(ij_Bas_pairs*kl_Bas_pairs+1))
nInBuf = max(nInBuf,(kl_Bas_pairs+1))
call mma_allocate(InBuf,nInBuf,Label='InBuf')

if ((IPR >= 5) .and. (case2 /= 0)) then
  write(u6,'(1X,4I2,2X,4I4,2X,4I4,2X,4I4)') iSym,jSym,kSym,lSym,iBas,jBas,kBas,lBas,iOrb,jOrb,kOrb,lOrb,iAsh,jAsh,kAsh,lAsh
end if

nPairs = 0
ij_pair = 0
nOff = 0
do i=1,iBas
  jMax = jBas
  if (iSym == jSym) jMax = i
  do j=1,jMax
    ij_pair = ij_pair+1

    ! read a block of electron repulsion integrals
    if (nPairs == 0) then
      iOpt = 1
      if (ij_pair /= 1) iOpt = 2
      call RdOrd(iRc,iOpt,iSym,jSym,kSym,lSym,InBuf,nInBuf,nPairs)
      nOff = 0
    end if

    ! unwrap triangular matrix of electron repulsion integrals
    if (kSym == lSym) then
      klBas = kBas*(kBas+1)/2
      call Square(InBuf(nOff+1),Scrt1,1,kBas,kBas)
      PQRS(1:kBas**2) => Scrt1(1:kBas**2)
    else
      klBas = kBas*lBas
      PQRS(1:klBas) => InBuf(nOff+1:nOff+klBas)
    end if
    PQRS_(1:klBas) => InBuf(nOff+1:nOff+klBas)

    ! generate Fock matrices
    if (case1 /= 4) then
      call Timing(Piaget_1,dum1,dum2,dum3)
      if (case1 == 2) then
        call Ftwo(case1,ExFac,iSym,kSym,i,j,off_sqMat,off_ltMat,D1I,FI,D1A,FA,PQRS)
      else
        call Ftwo(case1,ExFac,iSym,jSym,i,j,off_sqMat,off_ltMat,D1I,FI,D1A,FA,PQRS)
      end if
      call Timing(Piaget_2,dum1,dum2,dum3)
      Piaget_2 = Piaget_2-Piaget_1
      Piaget_3 = Piaget_3+Piaget_2
    end if

    ! first half transformation of electron repulsion integrals:
    ! (ij!kl) --> (vx!ij)
    if (case2 /= 0) then
      call Timing(Candino_1,dum1,dum2,dum3)
      call Tra2A(ij_pair,ij_Bas_pairs,kl_Orb_pairs,kSym,lSym,kBas,lBas,kAsh,lAsh,CMO(kkkOff),CMO(lllOff),PQRS,Buf2,Buf3,PQVX)
      if (Process_Twice) then
        call Tra2C(i,iSym,iBas,iAsh,j,jSym,jBas,jAsh,kl_Bas_pairs,ij_Orb_pairs,CMO(iiiOff),CMO(jjjOff),PQRS_,Buf2,TURS)
      end if
      call Timing(Candino_2,dum1,dum2,dum3)
      Candino_2 = Candino_2-Candino_1
      Candino_3 = Candino_3+Candino_2
    end if
    PQRS => null()
    PQRS_ => null()

    nPairs = nPairs-1
    nOff = nOff+kl_Bas_pairs

  end do
end do
call Timing(Candino_1,dum1,dum2,dum3)
if (IPR >= 99) then
  call RecPrt('PQVX',' ',PQVX,ij_Bas_Pairs,kl_Orb_Pairs)
  if (Process_Twice) call RecPrt('TURS',' ',TURS,kl_Bas_Pairs,ij_Orb_Pairs)
end if

! Do second half transformation, skip if block is zero long

if (case2 /= 0) then
  i1 = off_PUVX(iSym,jSym,kSym)
  i2 = off_PUVX(jSym,iSym,kSym)
  !call FZero(PUVX(1+i1),iOrb*jAsh*kl_Orb_pairs)
  !call FZero(PUVX(1+i2),jOrb*iAsh*kl_Orb_pairs)
  do kl_pair=1,kl_Orb_pairs
    iOff = (kl_pair-1)*ij_Bas_pairs

    ! unwrap triangular matrix of electron repulsion integrals
    if (iSym == jSym) then
      call Square(PQVX(iOff+1),Scrt1,1,iBas,iBas)
      Buf9(1:iBas**2) => Scrt1(1:iBas**2)
    else
      Buf9(1:iBas*jBas) => PQVX(iOff+1:iOff+iBas*jBas)
    end if

    ! second half transformation of electron repulsion integrals:
    ! (vx!ij) --> (pu!vx)
    call Tra2B(iSym,jSym,iBas,jBas,iAsh,jAsh,iOrb,jOrb,kl_pair,kl_Orb_pairs,CMO(iiOff),CMO(jjOff),CMO(iiiOff),CMO(jjjOff),Buf9, &
               Buf2,Buf3,Buf3,PUVX(1+i1),PUVX(1+i2))
    Buf9 => null()

  end do
  if (IPR >= 99) then
    call RecPrt('PUVX(i,j,k)',' ',PUVX(1+i1),iOrb,jAsh*kl_Orb_pairs)
    call RecPrt('PUVX(j,i,k)',' ',PUVX(1+i2),jOrb,iAsh*kl_Orb_pairs)
  end if

  if (Process_Twice) then
    i1 = off_PUVX(kSym,lSym,iSym)
    i2 = off_PUVX(lSym,kSym,iSym)
    !call FZero(PUVX(1+i1),kOrb*iAsh*ij_Orb_pairs)
    !call FZero(PUVX(1+i2),lOrb*kAsh*ij_Orb_pairs)
    do ij_pair=1,ij_Orb_pairs
      iOff = (ij_pair-1)*kl_Bas_pairs

      ! unwrap triangular matrix of electron repulsion integrals
      if (kSym == lSym) then
        call Square(TURS(iOff+1),Scrt1,1,kBas,kBas)
        Buf9(1:kbas**2) => Scrt1(1:kbas**2)
      else
        Buf9(1:kBas*lBas) => TURS(iOff+1:iOff+kBas*lBas)
      end if

      ! econd half transformation of electron repulsion integrals
      ! vx!ij) --> (pu!vx)
      call Tra2B(kSym,lSym,kBas,lBas,kAsh,lAsh,kOrb,lOrb,ij_pair,ij_Orb_pairs,CMO(kkOff),CMO(llOff),CMO(kkkOff),CMO(lllOff),Buf9, &
                 Buf2,Buf3,Buf3,PUVX(1+i1),PUVX(1+i2))
      Buf9 => null()
    end do
    if (IPR >= 99) then
      call RecPrt('PUVX(k,l,i)',' ',PUVX(1+i1),kOrb,lAsh*ij_Orb_pairs)
      call RecPrt('PUVX(l,k,i)',' ',PUVX(1+i2),lOrb,kAsh*ij_Orb_pairs)
    end if
  end if

end if

! deallocate memory
call mma_deallocate(InBuf)
if (nTURS /= 0) call mma_deallocate(TURS)
call mma_deallocate(PQVX)
call mma_deallocate(Buf3)
call mma_deallocate(Buf2)
if (allocated(Scrt1)) call mma_deallocate(Scrt1)
call Timing(Candino_2,dum1,dum2,dum3)
Candino_2 = Candino_2-Candino_1
Candino_3 = Candino_3+Candino_2

return

end subroutine TraDrv
