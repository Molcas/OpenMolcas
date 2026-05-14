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

subroutine getmo_cvb(cmo,ic)

use casvb_global, only: nbas_mo
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: cmo(*)
integer(kind=iwp), intent(in) :: ic
real(kind=wp), allocatable :: cmo2(:,:)

call mma_allocate(cmo2,nbas_mo,nbas_mo,label='cmo2')
if (ic <= 1) then
  call getmo2_cvb(cmo,cmo2,ic)
else
  call getmo2_cvb(cmo2,cmo,ic)
end if
call mma_deallocate(cmo2)

return

end subroutine getmo_cvb
