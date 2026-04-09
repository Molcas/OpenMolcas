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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2016,2017, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine GrdClc(Do_All)
!***********************************************************************
!                                                                      *
!     purpose: Compute gradients and write on disk.                    *
!                                                                      *
!                                                                      *
!     input:                                                           *
!       Do_All  : variable telling what gradients compute: .true. -    *
!                 all gradients, .False. - last gradient               *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!***********************************************************************

use InfSCF, only: CMO_Ref, FockMO, Iter, Iter_Start, kOV, mOV, nBO, nBT, nD, nOO, OneHam, Ovrlp
use LnkLst, only: LLGrad, LLlGrd, PutVec
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
logical(kind=iwp), intent(inout) :: Do_All
integer(kind=iwp) :: iOpt, LpStrt
real(kind=wp), allocatable :: GrdOO(:,:), GrdOV(:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

! Allocate memory for gradients
call mma_allocate(GrdOO,nOO,nD,Label='GrdOO')
call mma_allocate(GrdOV,mOV,Label='GrdOV')

! Find the beginning of the loop
if (Do_All) then
  LpStrt = Iter_Start
  Do_All = .false.
else
  LpStrt = Iter
end if

! Compute all gradients / last gradient

do iOpt=LpStrt,Iter

  if (iOpt == Iter) call Mk_FockMO(OneHam,Ovrlp,nBT,CMO_Ref,nBO,FockMO,nOO,nD,iOpt)

  call EGrad(OneHam,Ovrlp,nBT,CMO_Ref,nBO,GrdOO,nOO,nD,iOpt)

  call vOO2OV(GrdOO,nOO,GrdOV,mOV,nD,kOV)

  ! Write Gradient to linked list

  call PutVec(GrdOV,mOV,iOpt,'OVWR',LLGrad)
  if (iOpt == Iter) call PutVec(GrdOV,mOV,iOpt,'OVWR',LLlGrd)

# ifdef _DEBUGPRINT_
  write(u6,*) 'GrdClc: Put Gradient iteration:',iOpt
  write(u6,*) 'iOpt,mOV=',iOpt,mOV
  call NrmClc(GrdOO,nOO*nD,'GrdClc','GrdOO')
  call NrmClc(GrdOV,mOV,'GrdClc','GrdOV')
# endif
end do

! Deallocate memory

call mma_deallocate(GrdOV)
call mma_deallocate(GrdOO)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine GrdClc
