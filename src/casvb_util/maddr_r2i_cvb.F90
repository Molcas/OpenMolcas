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

function maddr_r2i_cvb(m_real_addr)

use casvb_global, only: memdebug
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: maddr_r2i_cvb
integer(kind=iwp) :: m_real_addr
#include "idbl_cvb.fh"
integer(kind=iwp) :: m_int_addr

m_int_addr = (m_real_addr-1)*idbl+1
maddr_r2i_cvb = m_int_addr
if (memdebug) then
  write(u6,*) ' maddr_r2i_cvb: real pointer :',m_real_addr
  write(u6,*) '                int pointer  :',m_int_addr
end if

return

end function maddr_r2i_cvb
