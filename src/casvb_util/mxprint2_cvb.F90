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

subroutine mxprint2_cvb(a,nrow,nrow2,ncol,itype)
! Prints matrix A, stored according to ITYPE

use casvb_global, only: formMXP1, formMXP3, iprec, iwidth
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: a(*)
integer(kind=iwp) :: nrow, nrow2, ncol, itype
integer(kind=iwp), parameter :: mxbuf = 8
integer(kind=iwp) :: i, ibuf(mxbuf), ii, ind, j, jend, jin, k, nbuf
real(kind=wp) :: buffer(mxbuf)

nbuf = min((iwidth-4)/(iprec+4),mxbuf)
if (nbuf == 7) nbuf = 6
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
  write(u6,formMXP1) (ibuf(i),i=1,jend-jin+1)
  do i=1,nrow
    k = 0
    do j=jin,jend
      k = k+1
      if (itype == 0) then
        ind = (j-1)*nrow2+i
      else if (itype == 1) then
        if (i >= j) then
          ind = i*(i-1)/2+j
        else
          ind = j*(j-1)/2+i
        end if
      else
        ind = (i-1)*nrow2+j
      end if
      buffer(k) = a(ind)
    end do
    write(u6,formMXP3) i,(buffer(ii),ii=1,jend-jin+1)
  end do
  jin = jend+1
  if (ncol <= nbuf) exit
end do

end subroutine mxprint2_cvb
