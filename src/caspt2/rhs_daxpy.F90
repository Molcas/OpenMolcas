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

subroutine RHS_DAXPY(NAS,NIS,ALPHA,lg_V1,lg_V2)
!SVC: this routine computes product ALPHA * V1 and adds to V2

use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: GA_Arrays

implicit none
integer(kind=iwp), intent(in) :: NAS, NIS
real(kind=wp), intent(in) :: ALPHA
integer(kind=iwp), intent(in) :: lg_V1, lg_V2
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
integer(kind=iwp) myRank, iLoV1, iHiV1, jLoV1, jHiV1, iLoV2, iHiV2, jLoV2, jHiV2, NV1, NV2, mV1, LDV1, mV2, LDV2

if (Is_Real_Par()) then
  myRank = GA_NodeID()
  call GA_Distribution(lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
  call GA_Distribution(lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
  if ((iLoV1 /= 0) .and. (iLoV2 /= 0)) then
    NV1 = (iHiV1-iLoV1+1)*(jHiV1-jLoV1+1)
    NV2 = (iHiV2-iLoV2+1)*(jHiV2-jLoV2+1)
    if (NV1 /= NV2) call AbEnd()
    call GA_Access(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
    call GA_Access(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
    ! V2 <- alpha*V1 + V2
    call DAXPY_(NV1,ALPHA,DBL_MB(mV1),1,DBL_MB(mV2),1)
    call GA_Release_Update(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
    call GA_Release(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
  end if
else
#endif
  call DAXPY_(NAS*NIS,ALPHA,GA_Arrays(lg_V1)%A,1,GA_Arrays(lg_V2)%A,1)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine RHS_DAXPY
