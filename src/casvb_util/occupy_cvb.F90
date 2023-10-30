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

subroutine occupy_cvb(nk,nel,locc,lunocc)

use Definitions, only: iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nel, nk(0:nel)
integer(kind=iwp), intent(_OUT_) :: locc(*), lunocc(*)
integer(kind=iwp) :: iel, iocc, iunocc

iocc = 0
iunocc = 0
do iel=1,nel
  if (nk(iel)-nk(iel-1) == 1) then
    iocc = iocc+1
    locc(iocc) = iel
  else if (nk(iel)-nk(iel-1) == 0) then
    iunocc = iunocc+1
    lunocc(iunocc) = iel
  else
    write(u6,*) ' Error in graphical indexing routine!'
    call abend_cvb()
  end if
end do

return

end subroutine occupy_cvb
