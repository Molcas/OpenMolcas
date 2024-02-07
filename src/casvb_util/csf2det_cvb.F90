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

subroutine csf2det_cvb(vec,detvec,isym_loc,iWay)

use csfbas, only: cts
use GLBBAS, only: DTOC
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: vec(*), detvec(*)
integer(kind=iwp), intent(in) :: isym_loc, iWay
#include "rasdim.fh"
#include "rasscf.fh"
integer(kind=iwp) :: jCopy

if (iWay == 1) then
  if (nac == 0) then
    detvec(1) = vec(1)
    return
  end if

  jCopy = 0
  call csdtvc(vec,detvec,iway,dtoc,cts,isym_loc,jcopy)
else if (iWay == 2) then
  if (nac == 0) then
    vec(1) = detvec(1)
    return
  end if

  jCopy = 0
  call csdtvc(vec,detvec,iway,dtoc,cts,isym_loc,jcopy)
end if

end subroutine csf2det_cvb
