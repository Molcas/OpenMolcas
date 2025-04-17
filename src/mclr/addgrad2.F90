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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine AddGrad2(rMat,fact)
! Purpose:
!   Adds the contribution from the gradient to
!    [2]
!   E   Kappa. This is done to insure us about
!   a beautiful convergence of the PCG,
!   which is just the case if E is symmetric.

use MCLR_Data, only: F0SQMO, ipCM, ipMat
use input_mclr, only: nSym, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: rMat(*)
real(kind=wp), intent(in) :: fact
integer(kind=iwp) :: iS
real(kind=wp), allocatable :: Temp(:)

do iS=1,nSym
  if (nOrb(is)*nOrb(is) == 0) cycle
  call mma_allocate(Temp,nOrb(is)**2,Label='Temp')
  ! T=Brillouin matrix
  call DGeSub(F0SQMO(ipCM(is)),nOrb(is),'N',F0SQMO(ipCM(is)),nOrb(is),'T',Temp,nOrb(is),nOrb(is),nOrb(is))
  !             t           t
  ! +1/2 { Kappa T - T kappa  }
  rMat(ipMat(is,is):ipMat(is,is)+nOrb(is)**2-1) = rMat(ipMat(is,is):ipMat(is,is)+nOrb(is)**2-1)-Two*Fact*Temp(1:nOrb(is)**2)
  call mma_deallocate(Temp)
end do

end subroutine AddGrad2
