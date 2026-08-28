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

subroutine faroald_free()
! The finalization subroutine lives outside of the
! faroald module so that it can be called separately.

use faroald, only: ex1_a, ex1_b, gtuvx, htu, ipat_of_det, occpat, mma_deallocate, ndoub, nsing, nSpin_Comb, nComb, ibComb, ictsdt
use faroald, only: conf_by_nopen, icnf_out, conf_reo, conf_arcw

implicit none

call mma_deallocate(ex1_a,safe='*')
call mma_deallocate(ex1_b,safe='*')
call mma_deallocate(htu,safe='*')
call mma_deallocate(gtuvx,safe='*')
call mma_deallocate(ipat_of_det,safe='*')
call mma_deallocate(occpat,safe='*')
call mma_deallocate(ndoub,safe='*')
call mma_deallocate(nsing,safe='*')
call mma_deallocate(nSpin_Comb,safe='*')
call mma_deallocate(nComb,safe='*')
call mma_deallocate(ibComb,safe='*')
call mma_deallocate(ictsdt,safe='*')
call mma_deallocate(conf_by_nopen,safe='*')
call mma_deallocate(icnf_out,safe='*')
call mma_deallocate(conf_reo,safe='*')
call mma_deallocate(conf_arcw,safe='*')

end subroutine faroald_free
