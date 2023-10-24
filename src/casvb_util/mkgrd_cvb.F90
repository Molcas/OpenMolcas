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

subroutine mkgrd_cvb(civb,civb2,grad,dvbdet,np,doorb)

use casvb_global, only: ndet, ndetvb, npr, nprorb, nvb, strucopt, vbdet
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: civb(0:ndet), civb2(0:ndet), grad(npr), dvbdet(ndetvb)
integer(kind=iwp) :: np
logical(kind=iwp) :: doorb
real(kind=wp), allocatable :: tmp(:)

call fzero(grad,nprorb)
if (doorb) call onedens_cvb(civb,civb2,grad,.false.,1)
if (strucopt) then
  call ci2vbg_cvb(civb2,dvbdet)
  if (np-nprorb == nvb) then
    call vb2strg_cvb(dvbdet,grad(nprorb+1))
  else if (np-nprorb < nvb) then
    call mma_allocate(tmp,nvb,label='tmp')
    call vb2strg_cvb(dvbdet,tmp)
    call fmove_cvb(tmp,vbdet,np-nprorb)
    call mma_deallocate(tmp)
  else
    write(u6,*) ' Error in mkgrd - np-nprorb > nvb :',np,nprorb,nvb
  end if
end if

return

end subroutine mkgrd_cvb
