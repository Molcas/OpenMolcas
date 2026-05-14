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

subroutine mkfn_cvb(fileid,ibf)

use casvb_global, only: ifilio, fileids, filename, max_rec, nrec, thresh_io
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: fileid
integer(kind=iwp), intent(out) :: ibf
integer(kind=iwp) :: i, ifile, irec
logical(kind=iwp) :: done
character(len=20) :: fn_tmp
logical(kind=iwp), parameter :: debug = .false.

done = .false.
do i=1,nrec
  if (abs(fileid-fileids(i)) < thresh_io) then
    ibf = i
    done = .true.
    exit
  end if
end do
if (.not. done) then
  nrec = nrec+1
  if (nrec > max_rec) then
    write(u6,*) ' nrec > max_rec in mkfn :',nrec,max_rec
    call abend_cvb()
  end if
  ibf = nrec
  ! generate new file name
  ! -> must be at most 8 characters to use daname
  fn_tmp = ' '
  irec = int(fileid)
  ifile = nint(10*(fileid-irec))
  call appendint_cvb(fn_tmp,irec,0)
  call appendint_cvb(fn_tmp,ifile,0)
  filename(ibf) = trim(fn_tmp)
  fileids(ibf) = fileid
  ifilio(ibf) = 0
end if
if (debug) then
  write(u6,*) ' IO information for identifier :',fileid
  write(u6,*) ' IBF is :',ibf
  write(u6,*) ' File name is :',filename(ibf)
end if

return

end subroutine mkfn_cvb
