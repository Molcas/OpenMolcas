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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine hess_cvb(vec)

use casvb_global, only: civb2, civb3, civb4, civb5, dvbdet, grad1, grad2, gradx, hessorb, hesst, icrit, iorts, n_hess, nfr, npr, &
                        orbinv, orbs, owrk2, sorbs, vec1
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: vec(nfr)
real(kind=wp), allocatable :: hess1(:), hess2(:)
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

n_hess = n_hess+1
if (.not. up2date_cvb('OOHESS')) then
  call make_cvb('OOHESS')
  call oohess_cvb(orbs,civb2,civb3,civb4,orbinv,sorbs,owrk2,grad2,gradx,hessorb,hesst)
end if
call mma_allocate(hess1,npr,label='hess1')
call mma_allocate(hess2,npr,label='hess2')
call free2all_cvb(vec,hess1,1)
if (icrit == 1) then
  call hess_svb1_cvb(orbs,civb2,civb3,civb4,civb5,orbinv,sorbs,owrk2,dvbdet,grad1,grad2,hessorb,vec1,iorts,hess1,hess2)
else if (icrit == 2) then
  call hess_evb1_cvb(orbs,civb2,civb3,civb4,sorbs,owrk2,dvbdet,grad1,grad2,hessorb,vec1,iorts,hess1,hess2)
end if
call all2free_cvb(hess2,vec,1)
call mma_deallocate(hess1)
call mma_deallocate(hess2)

return

end subroutine hess_cvb
