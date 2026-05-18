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
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

function RHS_DDOT(NAS,NIS,lg_V1,lg_V2)
!SVC: this routine computes the DDOT of the RHS arrays V1 and V2

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: RHS_DDOT
integer(kind=iwp), intent(in) :: NAS, NIS, lg_V1, lg_V2
real(kind=wp), external :: DDot_
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"

if (Is_Real_Par()) then
  RHS_DDOT = GA_DDOT(lg_V1,lg_V2)
else
#endif
  RHS_DDOT = DDOT_(NAS*NIS,GA_Arrays(lg_V1)%A,1,GA_Arrays(lg_V2)%A,1)
#ifdef _MOLCAS_MPP_
end if
#endif

end function RHS_DDOT
