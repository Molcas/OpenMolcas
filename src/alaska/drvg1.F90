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
! Copyright (C) 1990-1992,2000, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drvg1(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals. The four outermost loops *
!          will control the type of the two-electron integral, eg.     *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis        *
!          functions of the requested type.                            *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified for gradient calculation. January '92           *
!             Modified for SetUp_Ints. January '00                     *
!***********************************************************************

use setup, only: mSkal, MxPrm
use k2_arrays, only: Sew_Scr
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use Para_Info, only: nProcs, King
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
integer(kind=iwp) :: i, iAng, ijMax, ijS, iOpt, iS, j, jAng, jS, kls, kS, lS, nHrrab, nij, nSkal
real(kind=wp) :: A_int, Cnt, P_Eff, ThrAO, TMax_all, TskHi, TskLw
logical(kind=iwp) :: DoFock, DoGrad, Indexation, lDummy, Triangular
character(len=8) :: Method_chk
integer(kind=iwp), allocatable :: Pair_Index(:,:)
real(kind=wp), allocatable :: TMax(:,:)
logical(kind=iwp), external :: Rsv_GTList
!*********** columbus interface ****************************************
integer(kind=iwp) :: Columbus
!                                                                      *
!***********************************************************************
!                                                                      *

Temp(:) = Zero

call StatusLine('Alaska: ','Computing 2-electron gradients')
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
!-- Precompute k2 entities.

call Get_iScalar('Columbus',Columbus)
Indexation = .false.
call Get_cArray('Relax Method',Method_chk,8)
! MP2 gradients:
if (Method_chk == 'MBPT2   ') Indexation = .true.
!*********** columbus interface ****************************************
! in order to access the half-sorted density matrix file
! some index arrays must be generated
if (Columbus == 1) Indexation = .true.
DoFock = .false.
DoGrad = .true.
ThrAO = Zero
call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
mSkal = nSkal
!                                                                      *
!***********************************************************************
!                                                                      *
!-- Prepare handling of two-particle density.

call PrepP()

!                                                                      *
!***********************************************************************
!                                                                      *
MxPrm = 0
do iAng=0,S%iAngMx
  MxPrm = max(MxPrm,S%MaxPrm(iAng))
end do
!
!***********************************************************************
!                                                                      *
!-- Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
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

call mma_allocate(Pair_Index,2,nskal*(nSkal+1)/2,Label='Ind_ij')
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
!-- Compute FLOPs for the transfer equation.

do iAng=0,S%iAngMx
  do jAng=0,iAng
    nHrrab = 0
    do i=0,iAng+1
      do j=0,jAng+1
        if (i+j <= iAng+jAng+1) then
          ijMax = min(iAng,jAng)+1
          nHrrab = nHrrab+ijMax*2+1
        end if
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
Triangular = .true.
call Init_TList(Triangular,P_Eff)
call Init_PPList()
call Init_GTList()
iOpt = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! In MPP case dispatch one processor to do 1-el gradients first

if ((nProcs > 1) .and. King()) then
  call Drvh1(Grad,Temp,nGrad)
  !if (nPrint(1) >= 15) call PrGrad(' Gradient excluding two-electron contribution',Grad,lDisp(0))
  Temp(:) = Zero
end if
!                                                                      *
!***********************************************************************
!                                                                      *
ijS = 0
klS = 0
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
  klS = kls-1
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
    if (A_Int < CutInt) cycle

    call Eval_g1_ijkl(iS,jS,kS,lS,Temp,nGrad,A_Int)

  end do

  ! Task endpoint
end do
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _MOLCAS_MPP_
call GADGOP(Temp,nGrad,'+')
#endif
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
call CloseP()
call Term_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()

!write(6,*) 'Exit Drvg1'

end subroutine Drvg1
