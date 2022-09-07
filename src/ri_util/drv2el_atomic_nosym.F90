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
! Copyright (C) 1990,1991,1993,1998, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv2El_Atomic_NoSym(Integral_WrOut,ThrAO,iCnttp,jCnttp,TInt,nTInt,In_Core,ADiag,LuA,ijS_req,Keep_Shell)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals.                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!***********************************************************************

use Basis_Info, only: nBas
use iSD_data, only: iSD
use Wrj12, only: SO2Ind
use k2_arrays, only: Sew_Scr
use Basis_Info, only: dbsc
use Gateway_global, only: force_out_of_core, iWROpt
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
external :: Integral_WrOut
real(kind=wp) :: ThrAO
integer(kind=iwp) :: iCnttp, jCnttp, nTInt, LuA, ijS_req, Keep_Shell
real(kind=wp), allocatable :: TInt(:), ADiag(:)
logical(kind=iwp) :: In_Core
#include "setup.fh"
#include "iTOffs.fh"
integer(kind=iwp) :: iAddr, iBfn, ij, ijAng, ijS, iS, iSeed, iTInt, iTOff, iWROpt_Save, ji, jS, jTInt, klAng, klS, kS, lS, MaxMem, &
                     MemLow, MemSew, MemT, mTInt, mTInt2, nBfn, nBfn_i, nBfn_j, nBfn_k, nBfn_l, nij, nIrrep_Save, nSkal, nTInt2
logical(kind=iwp) :: Do_ERIs, Do_RI_Basis, DoFock, DoGrad, FreeK2, Indexation, Only_DB, Out_of_Core, Verbose
integer(kind=iwp), allocatable :: IJInd(:,:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
! Temporary modifications to facilitate atomic calculations

nIrrep_Save = nIrrep
nIrrep = 1
iWROpt_Save = iWROpt
iWROpt = 1

Do_RI_Basis = dbsc(iCnttp)%Aux

call Set_Basis_Mode_Atomic(iCnttp,jCnttp)
call Setup_iSD()

if (Do_RI_Basis .and. (ijS_req == 0)) then
  call WarningMessage(2,'Do_RI_Basis .and. (ijS_req == 0)')
  call Abend()
end if
!write(u6,*) 'Do_RI_Basis=',Do_RI_Basis
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation. Do not generate
! tables for indexation.

DoGrad = .false.
DoFock = .false.
Indexation = .false.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of pairs

call mma_allocate(IJInd,2,nSkal*(nSkal+1)/2,Label='IJInd')
nij = 0
nBfn = 0
if (Do_RI_Basis) then
  iS = nSkal   ! Dummy shell
  do jS=1,nSkal-1
    nij = nij+1
    IJInd(1,nij) = iS
    IJInd(2,nij) = jS
    nBfn = nBfn+iSD(2,jS)*iSD(3,jS)
  end do
else
  do iS=1,nSkal
    do jS=1,iS
      nij = nij+1
      IJInd(1,nij) = iS
      IJInd(2,nij) = jS
    end do
  end do
end if
!write(u6,*) 'nij=',nij
!                                                                      *
!***********************************************************************
!                                                                      *
! Preallocate some core for Seward!

call mma_MaxDBLE(MemSew)
MemLow = min(MemSew/2,1024*128)
MemSew = max(MemSew/10,MemLow)
call mma_allocate(Sew_Scr,MemSew,Label='Sew_Scr')
!                                                                      *
!***********************************************************************
!                                                                      *
! Determine if only diagonal block should be computed.
! This option is forced when the auxiliary basis set is transformed
! to the Cholesky basis!

Only_DB = (ijS_req /= 0) .or. Do_RI_Basis
!                                                                      *
!***********************************************************************
!                                                                      *
! Choose between in-core and out-of-core options

call mma_MaxDBLE(MemT)
MemT = MemT/2

if (Only_DB) then

  ! Only diagonal block

  nTInt = 0
  if (Do_RI_Basis) then
    do iS=1,nSkal-1 ! Skip the dummy shell
      nBfn_i = iSD(2,iS)*iSD(3,iS)
      if (iS == ijS_req) nTInt = nBfn_i
    end do
    if (nTInt == 0) then
      call WarningMessage(2,'Drv2el_atomic_nosym: nTInt == 0')
      call Abend()
    end if

    call mma_Allocate(SO2Ind,nBfn,Label='SO2Ind')
    do iBfn=1,nBfn
      SO2Ind(iBfn) = iBfn
    end do
    nSOs = nBfn
  else
    do iS=1,nSkal
      nBfn_i = iSD(2,iS)*iSD(3,iS)
      do jS=1,iS-1
        ijS = iS*(iS-1)/2+jS
        nBfn_j = iSD(2,jS)*iSD(3,jS)
        if (ijS == ijS_req) nTInt = nBfn_i*nBfn_j
      end do
      ijS = iS*(iS+1)/2
      if (ijS == ijS_req) nTInt = nBfn_i*(nBfn_i+1)/2
    end do
  end if
  mTInt = nTInt
  nTInt2 = mTInt*nTInt
  if (nTInt2 > MemT) then
    call WarningMessage(2,'Not enough memory!')
    call Abend()
  end if

  iTOffs(1) = 0

  In_Core = .true.       ! no out-of-core option needed.
  Out_of_core = .false.
  mTInt2 = nTInt2

else

  ! All blocks

  nTInt = nBas(0)*(nBas(0)+1)/2
  mTInt = nTInt
  nTInt2 = nTInt**2
  In_Core = nTInt2 <= MemT
  if (Force_out_of_Core) In_Core = .false.
  Out_of_Core = .not. In_Core

  ! Compute the size of the array, TInt, to write the integrals to.

  if (Out_of_Core) then

    ! Find the larges block for a fixed shell pair.

    mTInt = 0
    do iS=1,nSkal
      nBfn_i = iSD(2,iS)*iSD(3,iS)
      do jS=1,iS-1
        nBfn_j = iSD(2,jS)*iSD(3,jS)
        mTInt = max(mTInt,nBfn_i*nBfn_j)
      end do
      mTInt = max(mTInt,nBfn_i*(nBfn_i+1)/2)
    end do
    nTInt2 = mTInt*nTInt

    ! Open file for the A-vectors

    iAddr = 0
    iSeed = 63
    LuA = IsFreeUnit(iSeed)
    call DaName_MF_WA(LuA,'AVEC0')

    call mma_allocate(ADiag,nTInt,label='ADiag')
    call FZero(ADiag,nTInt)

  else

    ! In-core option

    mTInt2 = nTInt2
    iTOffs(1) = 0  ! Offset permanently set to zero!

  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
iTOffs(2) = nTInt  ! # of rows in TInt
iTOffs(3) = mTInt  ! # of colums in TInt
call mma_allocate(TInt,nTInt2,label='TInt')
if (In_Core) call FZero(TInt,nTInt2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Now do a quadruple loop over shells

iTOff = 0
iTOffs(4) = 0    ! Offset to the ij set
do ijS=1,nij
  iS = IJInd(1,ijS)
  jS = IJInd(2,ijS)

  nBfn_i = iSD(2,iS)*iSD(3,iS)
  nBfn_j = iSD(2,jS)*iSD(3,jS)
  ijAng = iSD(1,iS)+iSD(1,jS)

  if (Out_of_Core) then
    if (iS == jS) then
      mTInt = nBfn_i*(nBfn_i+1)/2
    else
      mTInt = nBfn_i*nBfn_j
    end if
    mTInt2 = mTInt*nTInt
    call FZero(TInt,mTInt2)
    iTOffs(1) = iTOff
    iTOffs(3) = mTInt
  end if

  iTOffs(5) = 0    ! Offset to the kl set
  do klS=1,ijS
    kS = IJInd(1,klS)
    lS = IJInd(2,klS)

    nBfn_k = iSD(2,kS)*iSD(3,kS)
    nBfn_l = iSD(2,lS)*iSD(3,lS)
    klAng = iSD(1,kS)+iSD(1,lS)

    ! For Only_DB compute the shell quadruplet (ijS_req|ijS_req)

    Do_ERIs = (.not. Only_DB) .or. ((ijS == ijS_req) .and. (klS == ijS_req))

    ! If high angular combination of the product basis is skipped
    ! do not compute the contributions.

    Do_ERIs = Do_ERIs .and. (ijAng <= Keep_Shell) .and. (klAng <= Keep_Shell)

    if (Do_ERIs) then
      call Eval_IJKL(iS,jS,kS,lS,TInt,mTInt2,Integral_WrOut)
    end if

    if (.not. Only_DB) then
      if (kS == lS) then
        iTOffs(5) = iTOffs(5)+nBfn_k*(nBfn_k+1)/2
      else
        iTOffs(5) = iTOffs(5)+nBfn_k*nBfn_l
      end if
    end if

  end do

  ! For out-of-core version write the integrals to disk!
  ! Pick up the diagonal elements.

  if (Out_of_Core) then
    call dDaFile(LuA,1,TInt,mTInt2,iAddr)
    call dcopy_(mTInt,TInt(1+iTOff),nTInt+1,ADiag(1+iTOff),1)
    iTOff = iTOff+mTInt
  end if

  if (.not. Only_DB) then
    if (iS == jS) then
      iTOffs(4) = iTOffs(4)+nBfn_i*(nBfn_i+1)/2
    else
      iTOffs(4) = iTOffs(4)+nBfn_i*nBfn_j
    end if
  end if

end do      ! ijS

if (Do_RI_Basis) call mma_deallocate(SO2Ind)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************

if (Out_of_Core) call mma_deallocate(TInt)
call mma_deallocate(IJInd)
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate integral environment.

Verbose = .false.
FreeK2 = .true.
call Term_Ints(Verbose,FreeK2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Square TInt from upper triangular to full.

if (In_Core .and. (.not. Do_RI_Basis)) then

  do iTInt=1,nTInt
    do jTInt=1,iTInt-1
      ij = (jTInt-1)*nTInt+iTInt
      ji = (iTInt-1)*nTInt+jTInt
      TInt(ij) = TInt(ji)
    end do
  end do
  !call RecPrt('Drv2El_atomic: TInt',' ',TInt,nTInt,nTInt)

else if (.not. Do_RI_Basis) then

  nij = nBas(0)*(nBas(0)+1)/2
  call mma_MaxDBLE(MaxMem)
  call Square_A(LuA,nij,MaxMem,Force_Out_of_Core)

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()
nIrrep = nIrrep_Save
iWROpt = iWROpt_Save
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Drv2El_Atomic_NoSym
