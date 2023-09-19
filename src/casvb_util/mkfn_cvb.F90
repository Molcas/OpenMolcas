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

implicit real*8(a-h,o-z)
#include "io_cvb.fh"
character*20 fn_tmp
logical debug, done
data debug/.false./

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
    write(6,*) ' nrec > max_rec in mkfn :',nrec,max_rec
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
  filename(ibf) = fn_tmp(1:len_trim_cvb(fn_tmp))
  fileids(ibf) = fileid
  ifilio(ibf) = 0
end if
if (debug) then
  write(6,*) ' IO information for identifier :',fileid
  write(6,*) ' IBF is :',ibf
  write(6,*) ' File name is :',filename(ibf)
end if

return

end subroutine mkfn_cvb
