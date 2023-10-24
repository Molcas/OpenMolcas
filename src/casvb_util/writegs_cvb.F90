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

subroutine writegs_cvb()

use casvb_global, only: cvb, cvbdet, iapr, ixapr, nalf, nbet, nda, ndetvb, norb, orbs, recn_tmp04
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ia, ib, idetvb, ioffs, ixa
integer(kind=iwp), allocatable :: iabind(:)

call mma_allocate(iabind,ndetvb,label='iabind')

call str2vbc_cvb(cvb,cvbdet)
ioffs = 0
call wris_cvb([ndetvb],1,recn_tmp04,ioffs)
call wris_cvb([norb],1,recn_tmp04,ioffs)
call wris_cvb([nalf],1,recn_tmp04,ioffs)
call wris_cvb([nbet],1,recn_tmp04,ioffs)
call wrrs_cvb(orbs,norb*norb,recn_tmp04,ioffs)
idetvb = 0
do ia=1,nda
  do ixa=ixapr(ia),ixapr(ia+1)-1
    idetvb = idetvb+1
    ib = iapr(ixa)
    iabind(idetvb) = ia+(ib-1)*nda
  end do
end do
call wris_cvb(iabind,ndetvb,recn_tmp04,ioffs)
call wrrs_cvb(cvbdet,ndetvb,recn_tmp04,ioffs)
call make_cvb('WRITEGS')

call mma_deallocate(iabind)

return

end subroutine writegs_cvb
