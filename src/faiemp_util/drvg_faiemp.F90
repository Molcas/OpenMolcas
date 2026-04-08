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
! Copyright (C) Ben Swerts                                             *
!***********************************************************************

subroutine Drvg_FAIEMP(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
!  Object: driver for the derivatives of central-fragment              *
!          two-electron integrals. The four outermost loops            *
!          will controll the type of the two-electron integral, eg.    *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis func-  *
!          tions of the requested type.                                *
!                                                                      *
!     Author: Ben Swerts                                               *
!                                                                      *
!     based on Drvg1                                                   *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use setup, only: mSkal, MxPrm
use k2_arrays, only: Sew_Scr
use Basis_Info, only: nBas, nBas_Frag
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Eight, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
integer(kind=iwp) :: i, iAng, ijS, iOpt, iS, jS, klS, kS, lS, nBas_Valence(0:7), nBT, nBVT, nij, nSkal, nSkal_Valence
real(kind=wp) :: A_int, Cnt, P_Eff, ThrAO, TMax_all, TskHi, TskLw
logical(kind=iwp) :: DoFock, DoGrad, Indexation, lDummy, lNoSkip, Triangular
integer(kind=iwp), allocatable :: Pair_Index(:,:)
real(kind=wp), allocatable :: TMax(:,:)
logical(kind=iwp), external :: Rsv_GTList

!                                                                      *
!***********************************************************************
!                                                                      *
! Handle both the valence and the fragment basis set

call Set_Basis_Mode('Valence')
call Nr_Shells(nSkal_Valence)
call Set_Basis_Mode('WithFragments')
call SetUp_iSD()
nBT = 0
nBVT = 0
do i=0,nIrrep-1
  nBas_Valence(i) = nBas(i)
  nBVT = nBVT+nTri_Elem(nBas(i))
  nBas(i) = nBas(i)+nBas_Frag(i)
  nBT = nBT+nTri_Elem(nBas(i))
end do
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
#ifdef _DEBUGPRINT_
write(u6,*) 'nSkal, nSkal_Valence, nSkal_Frag = ',nSkal,nSkal_Valence,nSkal-nSkal_Valence
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Prepare handling of the combined two-particle density.

call PrepP_FAIEMP(nBas_Valence,nBT,nBVT)
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

call mma_allocate(TMax,nSkal,nSkal,label='TMax')
call Shell_MxSchwz(nSkal,TMax)
TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(Pair_Index,2,nTri_Elem(nSkal),label='Pair_Index')
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      Pair_Index(1,nij) = iS
      Pair_Index(2,nij) = jS
    end if
  end do
end do
P_Eff = real(nij,kind=wp)
!                                                                      *
!***********************************************************************
!                                                                      *
Triangular = .true.
call Init_TList(Triangular,P_Eff)
call Init_PPList()
call Init_GTList()
iOpt = 0
Temp(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
do
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.
  if (.not. Rsv_GTList(TskLw,TskHi,iOpt,lDummy)) exit

  ! Now do a quadruple loop over shells

  ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
  klS = int(TskLw-real(ijS,kind=wp)*(real(ijS,kind=wp)-One)/Two)
  Cnt = TskLw

  Cnt = Cnt-One
  klS = klS-1
  do
    Cnt = Cnt+One
    if (Cnt-TskHi > 1.0e-10_wp) exit
    klS = klS+1
    if (klS > ijS) then
      ijS = ijS+1
      klS = 1
    end if
    iS = Pair_Index(1,ijS)
    jS = Pair_Index(2,ijS)
    kS = Pair_Index(1,klS)
    lS = Pair_Index(2,klS)

    A_int = TMax(iS,jS)*TMax(kS,lS)

    lNoSkip = A_Int >= CutInt
    ! only calculate needed integrals and only update the valence part of the
    ! Fock matrix (iS > nSkal_Valence, lS <= nSkal_Valence, jS and kS
    ! belonging to different regions)
    if (jS <= nSkal_Valence) then
      lNoSkip = lNoSkip .and. (kS > nSkal_Valence)
    else
      lNoSkip = lNoSkip .and. (kS <= nSkal_Valence)
    end if
    lNoSkip = lNoSkip .and. (lS <= nSkal_Valence)
    if (.not. lNoSkip) cycle

    call Eval_g1_ijkl(iS,jS,kS,lS,Temp,nGrad,A_Int)

  end do

end do
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Sew_Scr,safe='*')
call Free_GTList()
call Free_PPList()
call Free_TList()
call mma_deallocate(Pair_Index)
call mma_deallocate(TMax)
!                                                                      *
!***********************************************************************
!                                                                      *
call Term_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *

! Accumulate the final results
call DScal_(nGrad,Half,Temp,1)
#ifdef _DEBUGPRINT_
call PrGrad('The FAIEMP 2-electron Contribution',Temp,nGrad)
#endif
call daxpy_(nGrad,One,Temp,1,Grad,1)

call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Drvg_FAIEMP
