!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Faro_CITrans_Free()

use Faroald_CITrans, only: Faro_graph_free, Faro_sym_free, Faro_ndet_group, Faro_ncsf_group, Faro_ndoc_group, Faro_nsoc_group, &
                           Faro_spintabs_free
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(Faro_ndoc_group,safe='*')
call mma_deallocate(Faro_nsoc_group,safe='*')
call mma_deallocate(Faro_ndet_group,safe='*')
call mma_deallocate(Faro_ncsf_group,safe='*')
call Faro_spintabs_free()
call Faro_graph_free()
call Faro_sym_free()

end subroutine Faro_CITrans_Free
