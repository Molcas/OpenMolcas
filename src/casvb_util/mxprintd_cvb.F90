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

subroutine mxprintd_cvb(a,nrow,ncol,itype)

use casvb_global, only: formMXP1, formMXP2, formMXP3, formMXP4

implicit real*8(a-h,o-z)
#include "print_cvb.fh"
parameter(mxbuf=8)
dimension buffer(mxbuf), ibuf(mxbuf), a(*)

nbuf = min((iwidth-4)/(iprec+8),mxbuf)
if (nbuf == 7) nbuf = 6
iform = 1
jin = 1
do
  jend = jin+nbuf-1
  if (ncol <= nbuf) jend = ncol
  if (jend > ncol+nbuf-1) return
  jend = min(ncol,jend)
  k = 0
  do j=jin,jend
    k = k+1
    ibuf(k) = j
  end do
  if (iform == 0) then
    write(6,formMXP1) (ibuf(i),i=1,jend-jin+1)
  else
    write(6,formMXP2) (ibuf(i),i=1,jend-jin+1)
  end if
  do i=1,nrow
    k = 0
    do j=jin,jend
      k = k+1
      if (itype == 0) then
        ind = (j-1)*nrow+i
      else if (itype == 1) then
        if (i >= j) then
          ind = i*(i-1)/2+j
        else
          ind = j*(j-1)/2+i
        end if
      else
        ind = (i-1)*ncol+j
      end if
      buffer(k) = a(ind)
    end do
    if (iform == 0) then
      write(6,formMXP3) i,(buffer(ii),ii=1,jend-jin+1)
    else
      write(6,formMXP4) i,(buffer(ii),ii=1,jend-jin+1)
    end if
  end do
  jin = jend+1
  if (ncol <= nbuf) exit
end do

return

end subroutine mxprintd_cvb
