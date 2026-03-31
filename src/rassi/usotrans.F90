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
! Copyright (C) 2019, Roland Lindh                                     *
!***********************************************************************

subroutine USOTRANS(USOR,USOI,NSS,EigVec,MSTATE,VSOR,VSOI)

use rassi_global_arrays, only: JBNUM
use Cntrl, only: MLTPLT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NSS, MSTATE
real(kind=wp) :: USOR(NSS,NSS), USOI(NSS,NSS), EigVec(MSTATE,MSTATE), VSOR(NSS,NSS), VSOI(NSS,NSS)
integer(kind=iwp) :: ISS, ISTATE, JOB, JSS, JSS_, KSS, KSS_, MPLET, MSPROJ
real(kind=wp) :: tmp_I, tmp_R
integer(kind=iwp), allocatable :: MAPST(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Before we start we need to backtransform the coefficients of the
! SO states from the basis of the SF states which diagonalize the
! SF Hamiltonian to the basis of the original SF states. This since
! all transition moments, whether or retrived from disk or
! recomputed, are in the basis of the original SF states.

! Mapping from spin states to spin-free state:
call mma_allocate(MAPST,nSS,3,Label='MAPST')
ISS = 0
do ISTATE=1,MSTATE
  JOB = JBNUM(ISTATE)
  MPLET = MLTPLT(JOB)
  do MSPROJ=-MPLET+1,MPLET-1,2
    ISS = ISS+1
    MAPST(ISS,1) = ISTATE
    MAPST(ISS,2) = MPLET
    MAPST(ISS,3) = MSPROJ
  end do
end do

! Let us transform the coefficients in USOR and USOI

do iSS=1,nSS
  do JSS=1,nSS
    tmp_R = Zero
    tmp_I = Zero
    jSS_ = MAPST(JSS,1)
    do kSS=1,nSS
      if (MAPST(kss,2) /= MAPST(jss,2)) cycle
      if (MAPST(kss,3) /= MAPST(jss,3)) cycle
      kSS_ = MAPST(kSS,1)
      tmp_R = tmp_R+USOR(kSS,iSS)*EigVec(jss_,kSS_)
      tmp_I = tmp_I+USOI(kSS,iSS)*EigVec(jss_,kSS_)
    end do
    VSOR(JSS,ISS) = tmp_R
    VSOI(JSS,ISS) = tmp_I
  end do
end do
call mma_deallocate(MAPST)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine USOTRANS
