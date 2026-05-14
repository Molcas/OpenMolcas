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

subroutine dafupd_cvb(lu,ioffset)

use Definitions, only: iwp, RtoI

implicit none
integer(kind=iwp), intent(in) :: lu, ioffset
integer(kind=iwp) :: ibuf(1000), ioff, mxddr, nwrite

mxddr = 1000
nwrite = 1000
ibuf(:) = 0

call iDaFile(lu,8,ibuf,nwrite,mxddr)

if (mxddr < ioffset) then
  ioff = mxddr
  do
    nwrite = min((ioffset-ioff)*RtoI,1000)
    call iDaFile(lu,1,ibuf,nwrite,ioff)
    if (ioff >= ioffset) exit
  end do
end if

return

end subroutine dafupd_cvb
