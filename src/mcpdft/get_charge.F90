!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2014, Matthew R. Hennefarth                            *
!***********************************************************************

!> @brief gets total molecular charge
!>
!> @author Matthew R. Hennefarth
function get_charge()
  use definitions,only:iwp,wp,u0
  use onedat,only:snoori
  use stdalloc,only:mma_allocate,mma_deallocate
  use general_data,only:nfro,nish,nactel,ntot1
  implicit none

#include "warnings.h"

  integer(kind=iwp) :: get_charge

  character(len=8) :: label
  integer(kind=iwp) :: rc,opt,comp,sylbl,nuc_charge,el_charge
  real(kind=wp),allocatable :: int1e_ovlp(:)

  rc = -1
  opt = ibset(0,sNoOri)
  comp = 1
  sylbl = 1
  label = 'Mltpl  0'

  el_charge = -2*sum(nfro+nish)-nactel

  call mma_allocate(int1e_ovlp,ntot1+4,label='int1e_ovlp')
  call rdone(rc,opt,label,comp,int1e_ovlp,sylbl)
  if(rc /= _RC_ALL_IS_WELL_) then
    call WarningMessage(2, 'Error computing system charge')
    write(u0,*) 'Error calling rdone'
    write(u0,*) 'Label = ', label
    write(u0,*) 'RC = ', rc
    call abend()
  endif

  nuc_charge = nint(int1e_ovlp(size(int1e_ovlp)))
  call mma_deallocate(int1e_ovlp)

  get_charge = nuc_charge+el_charge
endfunction get_charge
