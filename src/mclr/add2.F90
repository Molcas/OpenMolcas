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
! Copyright (C) 1991, Anders Bernhardsson                              *
!***********************************************************************

subroutine Add2(rMat,fact)
! Purpose:
!   Adds the contribution from the gradient to
!    [2]
!   E   Kappa. This is done to insure us about
!   a beautifull convergence of the PCG,
!   which is just the case if E is symmetric.

use MCLR_Data, only: SFock
use MCLR_data, only: ipCM, ipMat
use input_mclr, only: nSym, nBas, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Four

implicit none
real*8 rMat(*)
real*8 fact
integer iS
real*8, allocatable :: Temp(:)

do iS=1,nSym
  if (nOrb(is)*nOrb(is) == 0) cycle
  call mma_allocate(Temp,nBas(is)**2,Label='Temp')
  ! T=Brillouin matrix
  call DGeSub(SFock(ipCM(is)),nOrb(is),'N',SFock(ipCM(is)),nOrb(is),'T',Temp,nOrb(is),nOrb(is),nOrb(is))
  !             t           t
  ! +1/2 { Kappa T - T kappa  }
  rMat(ipMat(is,is):ipMat(is,is)+nOrb(is)**2-1) = rMat(ipMat(is,is):ipMat(is,is)+nOrb(is)**2-1)-Four*Fact*Temp(1:nOrb(is)**2)
  call mma_deallocate(Temp)
end do

end subroutine Add2
