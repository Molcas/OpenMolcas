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

subroutine Faro_Free()
! The finalization subroutine lives outside of the
! Faroald module so that it can be called separately.

use Faroald, only: Faro_ex1_a, Faro_ex1_b, Faro_gtuvx, Faro_htu, Faro_ipat_of_det, Faro_occpat, mma_deallocate, Faro_ndoub, &
                   Faro_nsing, Faro_nSpin_Comb, Faro_nComb, Faro_ibComb, Faro_ictsdt

implicit none

call mma_deallocate(Faro_ex1_a,safe='*')
call mma_deallocate(Faro_ex1_b,safe='*')
call mma_deallocate(Faro_htu,safe='*')
call mma_deallocate(Faro_gtuvx,safe='*')
call mma_deallocate(Faro_ipat_of_det,safe='*')
call mma_deallocate(Faro_occpat,safe='*')
call mma_deallocate(Faro_ndoub,safe='*')
call mma_deallocate(Faro_nsing,safe='*')
call mma_deallocate(Faro_nSpin_Comb,safe='*')
call mma_deallocate(Faro_nComb,safe='*')
call mma_deallocate(Faro_ibComb,safe='*')
call mma_deallocate(Faro_ictsdt,safe='*')

end subroutine Faro_Free
