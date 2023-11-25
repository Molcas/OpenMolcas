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

subroutine psym1_cvb(civec1,civec2,osym,ientry)

use casvb_global, only: mxirrep, nalf, nbet, nda, ndb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: civec1(nda,ndb), osym(mxirrep)
real(kind=wp), intent(in) :: civec2(nda,ndb)
integer(kind=iwp), intent(in) :: ientry
integer(kind=iwp) :: iasyind(0:mxirrep), ibsyind(0:mxirrep), irpdet(mxirrep)
integer(kind=iwp), allocatable :: isymalf(:), isymbet(:)

call mma_allocate(isymalf,nda,label='isymalf')
call mma_allocate(isymbet,ndb,label='isymbet')

call symgen_cvb(nalf,nbet,nda,ndb,isymalf,isymbet,iasyind,ibsyind,irpdet)

call psym2_cvb(civec1,civec2,isymalf,isymbet,iasyind,ibsyind,osym,ientry)

call mma_deallocate(isymalf)
call mma_deallocate(isymbet)

return

end subroutine psym1_cvb
