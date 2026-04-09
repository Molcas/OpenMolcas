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
! Copyright (C) 2014, Matthew R. Hennefarth                            *
!***********************************************************************
!  get_charge
!
!> @brief gets total molecular charge
!>
!> @author Matthew R. Hennefarth
!***********************************************************************

function get_charge()

use onedat, only: snoori
use general_data, only: nactel, nfro, nish, ntot1
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u0

implicit none
integer(kind=iwp) :: get_charge
#include "warnings.h"
integer(kind=iwp) :: comp, el_charge, nuc_charge, opt, rc, sylbl
character(len=8) :: label
real(kind=wp), allocatable :: int1e_ovlp(:)

rc = -1
opt = ibset(0,sNoOri)
comp = 1
sylbl = 1
label = 'Mltpl  0'

el_charge = -2*sum(nfro+nish)-nactel

call mma_allocate(int1e_ovlp,ntot1+4,label='int1e_ovlp')
call rdone(rc,opt,label,comp,int1e_ovlp,sylbl)
if (rc /= _RC_ALL_IS_WELL_) then
  call WarningMessage(2,'Error computing system charge')
  write(u0,*) 'Error calling rdone'
  write(u0,*) 'Label = ',label
  write(u0,*) 'RC = ',rc
  call abend()
end if

nuc_charge = nint(int1e_ovlp(size(int1e_ovlp)))
call mma_deallocate(int1e_ovlp)

get_charge = nuc_charge+el_charge

end function get_charge
