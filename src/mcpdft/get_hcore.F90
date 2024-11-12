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

!> @brief Load bar nuclei Hamiltonian
!>
!> @details
!>   This is h_{mu, nu}, the electron kinetic and nuclear-electron attraction integrals summed in the AO basis.
!>
!> @author Matthew R. Hennefarth
!>
!> @param[out] hcore one electron integrals (kinetic energy + nuclear-electron attraction)
subroutine get_hcore(hcore)
  use definitions,only:iwp,wp,u0
  use onedat,only:snoori,snonuc
  use stdalloc,only:mma_allocate,mma_deallocate
  implicit none

#include "warnings.h"

  real(kind=wp), intent(out) :: hcore(*)

  character(len=8) :: label
  integer(kind=iwp) :: rc,opt,comp,sylbl

  rc = -1
  opt = ibset(ibset(0,sNoOri), sNoNuc)
  comp = 1
  sylbl = 1
  label = 'OneHam  '

  call rdone(rc,opt,label,comp,hcore,sylbl)
  if(rc /= _RC_ALL_IS_WELL_) then
    call WarningMessage(2, 'Error loading hcore integrals')
    write(u0,*) 'Error calling rdone'
    write(u0,*) 'Label = ', label
    write(u0,*) 'RC = ', rc
    call abend()
  endif

endsubroutine get_hcore
