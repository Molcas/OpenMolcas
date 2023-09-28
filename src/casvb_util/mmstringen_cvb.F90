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

subroutine mmstringen_cvb(norb,nel,locc,lunocc,nkmin,nkmax)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: norb, nel, locc(*), lunocc(*), nkmin(0:norb), nkmax(0:norb)
#include "WrkSpc.fh"
integer(kind=iwp) :: i_locc, i_lunocc, i_nk, indx, rc
integer(kind=iwp), external :: mstacki_cvb

i_nk = mstacki_cvb(norb+1)
! Spin string loop initialization (use xdet as graph storage):
call imove_cvb(nkmax,iwork(i_nk),norb+1)
! Spin string loop starts here:
indx = 0
do
  indx = indx+1
  i_locc = (indx-1)*nel+1
  i_lunocc = (indx-1)*(norb-nel)+1
  call occupy_cvb(iwork(i_nk),norb,locc(i_locc),lunocc(i_lunocc))
  call loop_cvb(norb,iwork(i_nk),nkmin,nkmax,rc)
  if (rc == 0) exit
end do
call mfreei_cvb(i_nk)

return

end subroutine mmstringen_cvb
