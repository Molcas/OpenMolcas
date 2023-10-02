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

subroutine setiaprtot_cvb()

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp

implicit none
#include "main_cvb.fh"
#include "WrkSpc.fh"
real(kind=wp) :: dum1(1), dum2(1), dum3
real(kind=wp), allocatable :: tmp(:,:)

call mma_allocate(tmp,nda,ndb,label='tmp')
call dpci2vb_cvb(tmp,dum1,dum2,0,dum3,4)
call setiaprtot2_cvb(tmp,iwork(ll(11)),iwork(ll(12)),iwork(ll(13)),iwork(ll(14)),npvb,nda,ndb)
call mma_deallocate(tmp)

return

end subroutine setiaprtot_cvb
