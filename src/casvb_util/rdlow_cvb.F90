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

subroutine rdlow_cvb(vec,n,fileid,ioffset)

use casvb_global, only: filename
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, ioffset
real(kind=wp), intent(out) :: vec(n)
real(kind=wp), intent(in) :: fileid
integer(kind=iwp) :: ibf, ioffs, lu
logical(kind=iwp) :: newfile
logical(kind=iwp), parameter :: debug = .false.

if (debug) write(u6,*) ' rdlow :',n,fileid,ioffset
call mkfn_cvb(fileid,ibf)
call ibf2unit_cvb(ibf,lu,newfile)

if (newfile) call ioopn_cvb(filename(ibf),lu)

ioffs = ioffset
call dDaFile(lu,2,vec,n,ioffs)

return

end subroutine rdlow_cvb
