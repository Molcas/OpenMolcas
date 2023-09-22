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

subroutine asonc12_cvb(c,sxc,nvec,citmp,orbs,gjorb,gjorb2,gjorb3,cvbdet)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: nvec
real(kind=wp) :: c(nvb,nvec), sxc(nvb,nvec), citmp(ndet), orbs(norb,norb), gjorb(*), gjorb2(*), gjorb3(*), cvbdet(ndetvb)
integer(kind=iwp) :: ivec

do ivec=1,nvec
  call str2vbf_cvb(c(1,ivec),cvbdet)
  call vb2cif_cvb(cvbdet,citmp)
  call applyts_cvb(citmp,orbs,gjorb,gjorb2,gjorb3)
  call ci2vbg_cvb(citmp,cvbdet)
  call vb2strg_cvb(cvbdet,sxc(1,ivec))
end do

return

end subroutine asonc12_cvb
