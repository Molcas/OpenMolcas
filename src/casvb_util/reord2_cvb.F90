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

subroutine reord2_cvb(cfrom,cto,imode)
! Front-end routine for molcas reord2, transforms
! from SGA CSFs to split-graph-GUGA CSFs.

use csfbas, only: conf
use glbbas, only: cftp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: cfrom(*)
real(kind=wp), intent(_OUT_) :: cto(*)
integer(kind=iwp), intent(in) :: imode
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
integer(kind=iwp), allocatable :: kcnf(:)

! NAC      rasscf.fh
! NACTEL   general.fh
! STSYM    general.fh
! IPR      rasscf.fh
call mma_allocate(kcnf,nactel,label='kcnf')
call reord2(nac,nactel,stsym,imode,conf,cftp,cfrom,cto,kcnf)
call mma_deallocate(kcnf)


end subroutine reord2_cvb
