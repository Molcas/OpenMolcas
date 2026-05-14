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

use Index_Functions, only: iTri
use casvb_global, only: formMXP1, formMXP2, formMXP3, formMXP4, iprec, iwidth
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: a(*)
integer(kind=iwp), intent(in) :: nrow, ncol, itype
integer(kind=iwp), parameter :: mxbuf = 8
integer(kind=iwp) :: i, ibuf(mxbuf), iform, ind, j, jend, jin, k, nbuf
real(kind=wp) :: buffer(mxbuf)

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
    write(u6,formMXP1) ibuf(1:jend-jin+1)
  else
    write(u6,formMXP2) ibuf(1:jend-jin+1)
  end if
  do i=1,nrow
    k = 0
    do j=jin,jend
      k = k+1
      if (itype == 0) then
        ind = (j-1)*nrow+i
      else if (itype == 1) then
        ind = iTri(i,j)
      else
        ind = (i-1)*ncol+j
      end if
      buffer(k) = a(ind)
    end do
    if (iform == 0) then
      write(u6,formMXP3) i,buffer(1:jend-jin+1)
    else
      write(u6,formMXP4) i,buffer(1:jend-jin+1)
    end if
  end do
  jin = jend+1
  if (ncol <= nbuf) exit
end do

return

end subroutine mxprintd_cvb
