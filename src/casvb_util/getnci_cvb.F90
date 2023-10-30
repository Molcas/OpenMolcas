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

!***********************************************************************
!*                                                                     *
!*  GETNCI    := Get number of CASSCF determinants in each irrep.      *
!*                                                                     *
!***********************************************************************
subroutine getnci_cvb(nciloc,nelloc,i2sloc,isymloc)

use casvb_global, only: mxirrep, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: nciloc(*)
integer(kind=iwp), intent(in) :: nelloc, i2sloc, isymloc
integer(kind=iwp) :: iasyind(0:mxirrep), ibsyind(0:mxirrep), irpdet(mxirrep), nalf1, nbet1, nda1, ndb1
integer(kind=iwp), allocatable :: isymalf(:), isymbet(:)

nalf1 = (nelloc+i2sloc)/2
nbet1 = nelloc-nalf1
call icomb_cvb(norb,nalf1,nda1)
call icomb_cvb(norb,nbet1,ndb1)
call mma_allocate(isymalf,nda1,label='isymalf')
call mma_allocate(isymbet,ndb1,label='isymbet')

call symgen_cvb(nalf1,nbet1,nda1,ndb1,isymalf,isymbet,iasyind,ibsyind,irpdet)

if (isymloc == 0) then
  nciloc(1:mxirrep) = irpdet(:)
else
  nciloc(1) = irpdet(isymloc)
end if

call mma_deallocate(isymalf)
call mma_deallocate(isymbet)

return

end subroutine getnci_cvb
