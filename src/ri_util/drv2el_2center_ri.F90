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
! Copyright (C) 1990,1991,1993,1998,2005, Roland Lindh                 *
!               1990, IBM                                              *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUGPRINT_
subroutine Drv2El_2Center_RI(ThrAO,A_Diag,MaxCntr)
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
!             Modified to 2-center ERIs for RI June '05                *
!***********************************************************************

use setup, only: nSOs
use Basis_Info, only: nBas_Aux
use iSD_data, only: iSO2Sh, nShBF
use RI_glob, only: iOffA, Lu_A, SO2Ind
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Int_Options, only: iTOffs
use Integral_interfaces, only: Int_PostProcess, int_wrout
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: ThrAO
real(kind=wp), allocatable, intent(out) :: A_Diag(:)
integer(kind=iwp), intent(out) :: MaxCntr
integer(kind=iwp) :: iAddr, iAddr_AQ(0:7), iIrrep, ip_A_n, ipAs_Diag, iS, iSeed, jS, kCol, kCol_Irrep(0:7), kS, lJ, lS, mB, &
                     MemLow, MemSew, nA_Diag, nB, nBfn2, nBfnTot, nSkal, nTInt, nTInt_, nZero
real(kind=wp) :: A_int, TCpu1, TCpu2, TMax_all, TWall1, TWall2
logical(kind=iwp) :: DoFock, DoGrad, Indexation
character(len=6) :: Name_Q
real(kind=wp), allocatable :: Scr(:), TInt(:), TMax(:), Tmp(:,:)
procedure(int_wrout) :: Integral_RI_2
integer(kind=iwp), external :: IsFreeUnit, nMemAm

!                                                                      *
!***********************************************************************
!                                                                      *
call StatusLine('Seward: ','Computing 2-center RI integrals')
!                                                                      *
!***********************************************************************
!                                                                      *
! Handle only the auxiliary basis set

call Set_Basis_Mode('Auxiliary')
call SetUp_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation.

DoGrad = .false.
DoFock = .false.
Indexation = .true.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)

call mma_Allocate(SO2Ind,nSOs,Label='SO2Ind')
call Mk_iSO2Ind(iSO2Sh,SO2Ind,nSOs,nSkal)

MaxCntr = 0

nBfn2 = 0
nBfnTot = 0
do iIrrep=0,nIrrep-1
  iTOffs(iIrrep+1) = nBfn2

  lJ = nBas_Aux(iIrrep)
  if (iIrrep == 0) lJ = lJ-1
  nBfn2 = nBfn2+lJ**2
  nBfnTot = nBfnTot+lJ
end do
nA_Diag = nBfnTot
call mma_allocate(A_Diag,nA_Diag,Label='A_Diag')
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,Label='TMax')
call mma_allocate(Tmp,nSkal,nSkal,Label='Tmp')
call Shell_MxSchwz(nSkal,Tmp)

TMax(:) = Tmp(:,nSkal)
TMax_all = Zero
do iS=1,nSkal
  TMax_all = max(TMax_all,TMax(iS))
end do

call mma_deallocate(Tmp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Preallocate some core for Seward!

call mma_maxDBLE(MemSew)
MemLow = min(MemSew/2,1024*128)
MemSew = max(MemSew/10,MemLow)
call xSetMem_Ints(MemSew)
!                                                                      *
!***********************************************************************
!                                                                      *
! Temporary buffer for computed integrals, compute the largest
! required buffer size and set up iOffA.

nTInt = 0
do jS=1,nSkal-1
  nTInt = max(nTInt,nMemAm(nShBF,nIrrep,nSkal-1,jS,iOffA,.true.))
end do
call mma_allocate(TInt,nTInt,Label='TInt')
call mma_allocate(Scr,nTInt,Label='Scr')
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)

! Open files for the A-vectors, set iAddr_AQ, kCol_iIrrep and iOffA(3,iIrrep).

nBfnTot = 0
do iIrrep=0,nIrrep-1
  iOffA(3,iIrrep) = nBfnTot
  mB = nBas_Aux(iIrrep)
  if (iIrrep == 0) mB = mB-1
  nBfnTot = nBfnTot+mB

  iSeed = 63+iIrrep
  Lu_A(iIrrep) = IsFreeUnit(iSeed)
  write(Name_Q,'(A4,I2.2)') 'AVEC',iIrrep
  if (mB /= 0) call DaName_MF_WA(Lu_A(iIrrep),Name_Q)

  iAddr_AQ(iIrrep) = 0
  kCol_Irrep(iIrrep) = 0
end do

Int_PostProcess => Integral_RI_2
iS = nSkal
kS = nSkal

do jS=1,nSkal-1
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! Initialize the buffer

  nTInt_ = nMemAm(nShBF,nIrrep,nSkal-1,jS,iOffA,.true.)
  TInt(1:nTInt_) = Zero
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! First compute the A matrix

  do lS=1,jS

    A_int = TMax(jS)*TMax(lS)
    if (A_Int >= CutInt) then
      call Eval_IJKL(iS,jS,kS,lS,Scr,nTInt_)
      TInt(1:nTInt_) = TInt(1:nTInt_)+Scr(1:nTInt_)
    end if

  end do ! lS
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! Write the A-vectors to disk

  do iIrrep=0,nIrrep-1
    mB = iOffA(2,iIrrep)   ! # of bf of shell jS
    if (mB /= 0) then

      ip_A_n = 1+iOffA(1,iIrrep)
      iAddr = iAddr_AQ(iIrrep) ! Disk address

      nB = nBas_Aux(iIrrep)
      if (iIrrep == 0) nB = nB-1 ! subtract dummy af
      do kCol=1+kCol_Irrep(iIrrep),mB+kCol_Irrep(iIrrep)

        ! Write the A-vector to file

        call dDaFile(Lu_A(iIrrep),1,TInt(ip_A_n),kCol,iAddr)

        ipAs_Diag = 1+iOffA(3,iIrrep)+kCol-1
        A_Diag(ipAs_Diag) = TInt(ip_A_n+kCol-1)
        nZero = nB-kCol
        if (nZero /= 0) call dDaFile(Lu_A(iIrrep),0,TInt(ip_A_n),nZero,iAddr)

        ip_A_n = ip_A_n+kCol
      end do

      kCol_Irrep(iIrrep) = kCol_Irrep(iIrrep)+mB
      iAddr_AQ(iIrrep) = iAddr
    end if

  end do ! iIrrep
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
end do   ! jS
!----------------------------------------------------------------------*
!                                                                      *
! Release the Seward core memory, the buffer, etc.
!
call Free_iSD()
call xRlsMem_Ints()
call mma_deallocate(Scr)
call mma_deallocate(TInt)
call mma_deallocate(TMax)
call mma_deallocate(SO2Ind)
nullify(Int_PostProcess)
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate integral environment.

call Term_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *

#undef _DEBUGPRINT_
end subroutine Drv2El_2Center_RI
