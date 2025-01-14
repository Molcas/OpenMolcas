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
! Copyright (C) 1990,1991,1993,1998,2025, Roland Lindh                 *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv2El_ijij(Pair_Index,nPairs,TMax,nSkal)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals.                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!             Modified driver for (ij|ij) integrals. Jan '25           *
!***********************************************************************

use iSD_data, only: iSD
use setup, only: mSkal
use Basis_Info, only: dbsc
use Gateway_Info, only: CutInt
use stdalloc, only: mma_allocate, mma_deallocate
use Integral_interfaces, only: Int_PostProcess, int_wrout
use Basis_Info, only: Shells
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in):: nPairs, nSkal
integer(kind=iwp), intent(in):: Pair_Index(2,nPairs)
real(kind=wp), intent(inout) :: TMax(nSkal,nSkal)

integer(kind=iwp) :: iCnttp, ijS, iS, jCnttp, jS, id_Tsk, iShll, jShll
real(kind=wp) :: A_int
character(len=72) :: SLine
real(kind=wp), allocatable :: TInt(:)
integer(kind=iwp), parameter :: nTInt = 1
logical(kind=iwp), external :: Rsv_Tsk
procedure(int_wrout) :: Integral_ijij

!                                                                      *
!***********************************************************************
!                                                                      *
SLine = 'Computing 2-electron integrals'
call StatusLine('Seward: ',SLine)
!                                                                      *
!***********************************************************************
!                                                                      *
Int_PostProcess =>  Integral_ijij
call mma_allocate(TInt,nTint,Label='TInt')
!                                                                      *
!***********************************************************************
!                                                                      *
Call Init_Tsk(id_Tsk,nPairs)

do
   If(.Not.Rsv_Tsk(id_Tsk,ijS)) exit
   iS = Pair_Index(1,ijS)
   jS = Pair_Index(2,ijS)

    iCnttp = iSD(13,iS)
    jCnttp = iSD(13,jS)
    if (dbsc(iCnttp)%fMass /= dbsc(jCnttp)%fMass) Cycle

    iShll = iSD(0,iS)

    ! In case of auxiliary basis sets we want the iS index to point at the dummay shell.
    if (Shells(iShll)%Aux .and. (iS /= mSkal)) cycle

    jShll = iSD(0,jS)

    ! In case the first shell is the dummy auxiliary basis shell make sure that
    ! the second shell, jS, also is a auxiliary basis shell.
    if (Shells(iShll)%Aux .and. (.not. Shells(jShll)%Aux)) cycle

    ! Make sure that the second shell never is the dummy auxiliary basis shell.
    if (Shells(jShll)%Aux .and. (jS == mSkal)) cycle

    A_int = TMax(iS,jS)*TMax(iS,jS)
    if (A_Int < CutInt) Cycle

    call Eval_IJKL(iS,jS,iS,jS,TInt,nTInt)
!   Write (6,*) iS, jS, Tmax(iS,jS), Sqrt(Abs(TInt(1)))
    TMax(iS,jS)=Sqrt(Abs(TInt(1)))
    TMax(jS,iS)=TMax(iS,jS)

end do

call Free_Tsk(id_Tsk)
call mma_deallocate(TInt)
nullify(Int_PostProcess)

end subroutine Drv2El_ijij
