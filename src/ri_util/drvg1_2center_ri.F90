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
! Copyright (C) 1990,1991,1992,2000,2007, Roland Lindh                 *
!               1990, IBM                                              *
!***********************************************************************

!#define _CD_TIMING_
subroutine Drvg1_2Center_RI(Grad,Temp,nGrad,ij2,nij_Eff)
!***********************************************************************
!                                                                      *
!  Object: driver for 2-center two-electron integrals in the RI scheme.*
!                                                                      *
!   The integral derivative is formulated as                           *
!   -Sum(ML) X_ij^K   V_LM^(1) X_kl^L  where                           *
!                                                                      *
!  X_ij^K = Sum(L) R_ij_L  Q_L^K                                       *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified for gradient calculation. January '92           *
!             Modified for SetUp_Ints. January '00                     *
!             Modified for 2-center RI gradients, January '07          *
!***********************************************************************

use setup, only: mSkal, MxPrm
use Index_Functions, only: nTri_Elem
use iSD_data, only: iSD
use pso_stuff, only: A_PT2
use k2_arrays, only: Sew_Scr
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use RICD_Info, only: Do_RI, RI_2C
use Symmetry_Info, only: nIrrep
use Para_Info, only: King, nProcs
use RI_glob, only: A, AMP2, CijK, DoCholExch, iMP2prpt, MxChVInShl, nIJR, nKvec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Eight, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad, nij_Eff, ij2(2,nij_Eff)
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
integer(kind=iwp) :: i, iAng, id, ij, iS, iSym1, iSym2, jDen, jlS, jS, jS_, kS, lA, lA_MP2, lS, lS_, mij, nij, nIJRMax, nSkal
real(kind=wp) :: A_int, ThrAO, TMax_all
logical(kind=iwp) :: DoFock, DoGrad, Indexation
integer(kind=iwp), allocatable :: Pair_Index(:,:)
real(kind=wp), allocatable :: TMax1(:), TMax2(:,:)
logical(kind=iwp), external :: Rsv_Tsk

!                                                                      *
!***********************************************************************
!                                                                      *
Temp(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Handle only the auxiliary basis set.

if (Do_RI) then
  call Set_Basis_Mode('Auxiliary')
  RI_2C = .true.
else
  call Set_Basis_Mode('Valence')
end if
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Precompute k2 entities.

Indexation = .false.
DoFock = .false.
DoGrad = .true.
ThrAO = Zero
call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
mSkal = nSkal
!                                                                      *
!***********************************************************************
!                                                                      *
MxPrm = 0
do iAng=0,S%iAngMx
  MxPrm = max(MxPrm,S%MaxPrm(iAng))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax2,nSkal,nSkal,Label='TMax2')
call Shell_MxSchwz(nSkal,TMax2)
TMax_all = Zero

if (Do_RI) then

  call mma_allocate(TMax1,nSkal,Label='TMax1')
  TMax1(:) = TMax2(:,nSkal)

  do iS=1,nSkal-1
    TMax_all = max(TMax_all,TMax1(iS))
  end do

else

  do ij=1,nij_Eff
    iS = ij2(1,ij)
    jS = ij2(2,ij)
    TMax_all = max(TMax_all,TMax2(iS,jS))
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate some scratch arrays to be used by the pget routines.
! In particular we will have temporary arrays for A_IJ and C_ijK.

! Lower case: valence basis set
! Upper case: auxiliary basis sets

if (DoCholExch) then

  ! Find the largest number of contractions in any given shell
  ! of auxiliary functions.

  MxChVInShl = 1
  if (Do_RI) then
    do i=1,nSkal
      MxChVInShl = max(MxChVInShl,iSD(3,i))
    end do
  else
    write(u6,*) 'Not Implemented for Cholesky yet!'
    call Abend()
  end if

  ! Scratch for A_IJ

  lA = MxChVInShl*MxChVInShl
  call mma_allocate(A,lA,Label='A')
  if (iMP2Prpt == 2) then
    lA_MP2 = MxChVInShl
    call mma_allocate(AMP2,lA_MP2,2,Label='AMP2')
  end if

  ! Find the largest set of ij. The basis i and j is due to the
  ! CD of the one-particle density.

  nIJRMax = 0
  do jDen=1,nKvec
    do iSym1=1,nIrrep
      do iSym2=1,nIrrep
        nIJRMax = max(nIJRMax,nIJR(iSym1,iSym2,jDen))
      end do
    end do
  end do

  ! Get scratch for C_kl^I and C_kl^J.
  ! Note that we need nDen arrays for C_kl^I and one for C_kl^J
  ! A_IJ = Sum(kl) C_kl^I x C_kl^J

  call mma_allocate(CijK,nIJRMax*MxChVInShl*(nKvec+1),Label='CijK')

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

if (Do_RI) then

  mij = nSkal-1
  call mma_allocate(Pair_Index,2,mij,Label='Shij')
  nij = 0
  do iS=1,nSkal-1
    if (TMax_All*TMax1(iS) >= CutInt) then
      nij = nij+1
      Pair_Index(1,nij) = nSkal
      Pair_Index(2,nij) = iS
    end if
  end do
else
  mij = nij_Eff
  call mma_allocate(Pair_Index,2,mij,Label='Shij')
  nij = 0
  do ij=1,nij_Eff
    iS = ij2(1,ij)
    jS = ij2(2,ij)
    if (TMax_All*TMax2(iS,jS) >= CutInt) then
      nij = nij+1
      Pair_Index(1,nij) = iS
      Pair_Index(2,nij) = jS
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! For a parallel implementation the iterations over shell-pairs
! are parallelized.

call Init_Tsk(id,nTri_Elem(nij))
!                                                                      *
!***********************************************************************
!                                                                      *
! In MPP case dispatch one processor to do 1-el gradients first

if ((nProcs > 1) .and. King()) then
  if (Do_RI) call Free_iSD()
  call Drvh1(Grad,Temp,nGrad)
# ifdef _DEBUGPRINT_
  call PrGrad(' Gradient excluding two-electron contribution',Grad,nGrad)
# endif
  Temp(:) = Zero
  if (Do_RI) then
    call Set_Basis_Mode('Auxiliary')
    call Setup_iSD()
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
do while (Rsv_Tsk(id,jlS))
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.

  ! Now do a quadruple loop over shells

  jS_ = int((One+sqrt(Eight*real(jlS,kind=wp)-Three))*Half)
  iS = Pair_Index(1,jS_)
  jS = Pair_Index(2,jS_)
  lS_ = int(real(jlS,kind=wp)-real(jS_,kind=wp)*(real(jS_,kind=wp)-One)*Half)
  kS = Pair_Index(1,lS_)
  lS = Pair_Index(2,lS_)

  if (Do_RI) then
    A_int = TMax1(jS)*TMax1(lS)
  else
    A_int = TMax2(iS,jS)*TMax2(kS,lS)
  end if
  if (A_Int < CutInt) cycle

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Eval_g1_ijkl(iS,jS,kS,lS,Temp,nGrad,A_Int)

end do
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
RI_2C = .false.
call mma_deallocate(Sew_Scr,safe='*')
call Free_Tsk(id)
call mma_deallocate(A_PT2,safe='*')
call mma_deallocate(Pair_Index)
call mma_deallocate(TMax1,safe='*')
call mma_deallocate(TMax2,safe='*')
!                                                                      *
!***********************************************************************
!                                                                      *
call Term_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
if (DoCholExch) then
  call mma_deallocate(CijK)
  call mma_deallocate(A)
end if
call mma_deallocate(AMP2,safe='*')

call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Drvg1_2Center_RI
