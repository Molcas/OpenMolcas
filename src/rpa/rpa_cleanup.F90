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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_Cleanup(irc)

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Clean up after RPA run (deallocate etc.)

use RPA_globals, only: CMO, EMO, l_CMO, l_EMO, l_OccEn, l_VirEn, OccEn, RPAModel, VirEn
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc

irc = 0

! Set "Relax Method" on Runfile
call Put_cArray('Relax Method',RPAModel,8)

! Deallocate memory
call mma_deallocate(CMO)
l_CMO = 0
call mma_deallocate(EMO)
l_EMO = 0
call mma_deallocate(OccEn)
l_OccEn(:) = 0
call mma_deallocate(VirEn)
l_VirEn(:) = 0

end subroutine RPA_Cleanup
