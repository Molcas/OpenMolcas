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

subroutine RHS_ADD(NAS,NIS,lg_W,W)
!SVC: this routine adds to the local part of a global RHS array the
!matching part of a replicate array.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NAS, NIS, lg_W
real(kind=wp), intent(in) :: W(NAS,*)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iHi, iLo, jHi, jLo, LDW, mW, myRank, NW
#include "global.fh"
#include "mafdecls.fh"

if (Is_Real_Par()) then
  myRank = GA_NodeID()
  call GA_Distribution(lg_W,myRank,iLo,iHi,jLo,jHi)
  if ((iLo /= 0) .and. (jLo /= 0)) then
    NW = (iHi-iLo+1)*(jHi-jLo+1)
    call GA_Access(lg_W,iLo,iHi,jLo,jHi,mW,LDW)
    call DAXPY_(NW,One,W(iLo,jLo),1,DBL_MB(mW),1)
    call GA_Release_Update(lg_W,iLo,iHi,jLo,jHi)
  end if
else
#endif
  call DAXPY_(NAS*NIS,One,W,1,GA_Arrays(lg_W)%A,1)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_ADD
