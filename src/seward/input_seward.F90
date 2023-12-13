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
! Copyright (C) 1990,1991,1993, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

subroutine Input_Seward(lOPTO)
!***********************************************************************
!                                                                      *
!     Object: to read the input to the integral package.               *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!                                                                      *
!             January '91 additional input for property calculations.  *
!             October '93 split up to RdCtl and SoCtl.                 *
!***********************************************************************

use Sizes_of_Seward, only: S
use Basis_Info, only: nBas
use Gateway_global, only: Test, PrPrt, Primitive_Pass
use Gateway_Info, only: Do_GuessOrb
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: lOPTO
#include "Molcas.fh"
#include "print.fh"
integer(kind=iwp) :: iRout
logical(kind=iwp), save :: Show_Save
logical(kind=iwp), external :: Reduce_Prt
character(len=LenIn), allocatable :: Mamn(:)
integer(kind=iwp), parameter :: nMamn = MaxBfn+MaxBfn_Aux

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
!                                                                      *
!***********************************************************************
!                                                                      *
if (Primitive_Pass) then
  Show_Save = Show
else
  Show = Show_Save
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Adjust the print level and some other parameters depending on
! if we are iterating or not.
!
! Set Show to false if Seward is run in property mode.
Show = Show .and. (.not. Prprt)

if (Reduce_Prt() .and. (nPrint(iRout) < 6) .and. (.not. Prprt)) then
  Show = .false.
  Do_GuessOrb = .false.
end if

Show = Show .and. (.not. Primitive_Pass)
Show = Show .or. Test
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify storage of basis functions to be in accordance with a
! calculation in the primitive or contracted basis.

call Flip_Flop(Primitive_Pass)
!                                                                      *
!***********************************************************************
!                                                                      *
! Start of output, collect all output to this routine!

if (Show) call Output1_Seward(lOPTO)
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Generate the SO or AO basis set

call mma_allocate(Mamn,nMamn,label='Mamn')
call SOCtl_Seward(Mamn,nMamn)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Test) then
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write information on the run file.

if (Primitive_Pass) then
  call Put_iArray('nBas_Prim',nBas,nIrrep)
  call Info2Runfile()
end if
call Put_cArray('Unique Basis Names',Mamn(1),(LENIN8)*S%nDim)
call Put_iArray('nBas',nBas,nIrrep)
call mma_deallocate(Mamn)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Input_Seward
