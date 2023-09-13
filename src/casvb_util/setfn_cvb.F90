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

subroutine setfn_cvb(fileid,fn)

implicit real*8(a-h,o-z)
#include "io_cvb.fh"
character*(*) fn
logical debug
data debug/.false./

lenfn = len_trim_cvb(fn)
do i=1,nrec
  if (fn(1:lenfn) == filename(i)) then
    fileid = fileids(i)
    return
  end if
end do
itry = 0
200 itry = itry+1
fileid_try = dble(itry)
do i=1,nrec
  if (fileid_try == fileids(i)) goto 200
end do
nrec = nrec+1
if (nrec > max_rec) then
  write(6,*) ' nrec > max_rec in setfn :',nrec,max_rec
  call abend_cvb()
end if
! set file name
filename(nrec) = fn
fileids(nrec) = fileid_try
ifilio(nrec) = 0
if (debug) then
  write(6,*) ' IO information for identifier :',fileids(nrec)
  write(6,*) ' IBF is :',nrec
  write(6,*) ' File name is :',filename(nrec)
end if
fileid = fileids(nrec)

return

end subroutine setfn_cvb
