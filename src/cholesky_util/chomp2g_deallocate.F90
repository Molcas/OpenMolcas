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
! Copyright (C) 2004,2005, Thomas Bondo Pedersen                       *
!               2010, Jonas Bostrom                                    *
!               2021, Roland Lindh                                     *
!***********************************************************************

subroutine ChoMP2g_deallocate(irc)
!
! Purpose: Deallocate memory needed for
!          MP2-gradients or properties.

use ChoMP2, only: AdrR1, AdrR2, ChoMP2g_allocated, EFrozT, EOccuT, EVirtT, MP2D, MP2D_e, MP2D_e_full, MP2D_full, MP2W, MP2W_e, &
                  MP2W_e_full, MP2W_full
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: iSym

irc = 0

if (ChoMP2g_allocated) then

  call mma_deallocate(MP2D_full)
  call mma_deallocate(MP2W_full)
  call mma_deallocate(MP2D_e_full)
  call mma_deallocate(MP2W_e_full)
  do iSym=1,8
    nullify(MP2D(iSym)%A)
    nullify(MP2W(iSym)%A)
    nullify(MP2D_e(iSym)%A)
    nullify(MP2W_e(iSym)%A)
  end do
  call mma_deallocate(AdrR2)
  call mma_deallocate(AdrR1)
  call mma_deallocate(EVirtT)
  call mma_deallocate(EOccuT)
  call mma_deallocate(EFrozT)

  ChoMP2g_allocated = .false.

end if

end subroutine ChoMP2g_deallocate
