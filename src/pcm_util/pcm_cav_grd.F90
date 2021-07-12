!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine PCM_Cav_grd(Grad,nGrad)

use PCM_arrays, only: dCntr, dPnt, dRad, dTes, PCM_N, PCM_SQ, PCMiSph, PCMSph, PCMTess
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(out) :: Grad(nGrad)
integer(kind=iwp) :: LcNAtm, MaxAto
real(kind=wp), allocatable :: DerDM(:,:), PCMGrd(:,:)
#include "print.fh"
#include "rctfld.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the geometric contributions to
! derivatives in solution

call mma_allocate(DerDM,nTs,nTs,label='DerDM')
call Get_nAtoms_All(MaxAto)
call mma_allocate(PCMGrd,3,MaxAto,label='PCMGrd')
LcNAtm = ISlPar(42)
call GeoDer(LcNAtm,Conductor,nTs,nS,Eps,PCMSph,PCMiSph,PCM_N,PCMTess,PCM_SQ,DerDM,PCMGrd,dTes,dPnt,dRad,dCntr)
!call RecPrt('PCM_Cav_Grd','(5G20.10)',PCMGrd,3,MaxAto)
call GrdTr_Alaska(PCMGrd,MaxAto,Grad,nGrad)
call mma_deallocate(DerDM)
call mma_deallocate(PCMGrd)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PCM_Cav_grd
