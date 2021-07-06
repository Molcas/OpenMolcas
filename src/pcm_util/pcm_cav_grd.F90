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
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(out) :: Grad(nGrad)
integer(kind=iwp) :: ip_DerDM, ip_PCMGrd, LcNAtm, MaxAto
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the geometric contributions to
! derivatives in solution

call GetMem('DerDM','Allo','Real',ip_DerDM,nTs*nTs)
call Get_nAtoms_All(MaxAto)
call GetMem('PCMGrd','Allo','Real',ip_PCMGrd,3*MaxAto)
LcNAtm = ISlPar(42)
call GeoDer(LcNAtm,Conductor,nTs,nS,Eps,PCMSph,PCMiSph,PCM_N,PCMTess,PCM_SQ,Work(ip_DerDM),Work(ip_PCMGrd),dTes,dPnt,dRad,dCntr)
!call RecPrt('PCM_Cav_Grd','(5G20.10)',Work(ip_PCMGrd),3,MaxAto)
call GrdTr_Alaska(Work(ip_PCMGrd),MaxAto,Grad,nGrad)
call GetMem('PCMGrd','Free','Real',ip_PCMGrd,3*MaxAto)
call GetMem('DerDM','Free','Real',ip_DerDM,nTs*nTs)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PCM_Cav_grd
