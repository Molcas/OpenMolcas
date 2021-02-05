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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

subroutine Finish_Kriging()

use kriging_mod

implicit none
#include "stdalloc.fh"

!write(6,*) 'Deallocating all kriging variables'
call mma_deallocate(Index_PGEK)
call Deallocate_protected()
call mma_deallocate(x0)
call mma_deallocate(full_R)
call mma_deallocate(full_RInv)
call mma_deallocate(dl)
call mma_deallocate(rl)
call mma_deallocate(Rones)
call mma_deallocate(kv)
call mma_deallocate(gpred)
call mma_deallocate(hpred)
call mma_deallocate(l)
call mma_deallocate(ll)
call mma_deallocate(cv)
call mma_deallocate(cvMatFder)
call mma_deallocate(cvMatSder)
call mma_deallocate(cvMatTder)
if (allocated(layer_U)) call mma_deallocate(layer_U)

return

end subroutine Finish_Kriging
