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

subroutine AddGrad(rKappa,rMat,idsym,fact)
! Purpose:
!   Adds the contribution from the gradient to
!    [2]
!   E   Kappa. This is done to insure us about
!   a beautifull convergence of the PCG,
!   which is just the case if E is symmetric.

use Symmetry_Info, only: Mul
use MCLR_Data, only: F0SQMO
use MCLR_Data, only: ipCM, ipMat
use input_mclr, only: nSym, nOrb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Half

implicit none
real*8 rkappa(*), rMat(*), fact
integer idsym
integer iS, jS
real*8, allocatable :: Tempi(:), Tempj(:)

do iS=1,nSym
  js = Mul(is,idsym)
  if (nOrb(is)*nOrb(js) == 0) cycle
  call mma_allocate(Tempi,nOrb(is)**2,Label='Tempi')
  call mma_allocate(Tempj,nOrb(js)**2,Label='Tempj')
  ! T=Brillouin matrix
  call DGeSub(F0SQMO(ipCM(is)),nOrb(is),'N',F0SQMO(ipCM(is)),nOrb(is),'T',Tempi,nOrb(is),nOrb(is),nOrb(is))
  call DGeSub(F0SQMO(ipCM(js)),nOrb(js),'N',F0SQMO(ipCM(js)),nOrb(js),'T',Tempj,nOrb(js),nOrb(js),nOrb(js))
  !             t           t
  ! +1/2 { Kappa T - T kappa  }
  call DGEMM_('T','N',nOrb(is),nOrb(js),nOrb(js),Half*fact,rkappa(ipMat(js,is)),nOrb(js),Tempj,nOrb(js),One,rMat(ipMat(is,js)), &
              nOrb(is))
  call DGEMM_('N','T',nOrb(is),nOrb(js),nOrb(is),-Half*fact,Tempi,nOrb(is),rKappa(ipmat(js,is)),nOrb(js),One,rMat(ipMat(is,js)), &
              nOrb(is))
  call mma_deallocate(Tempi)
  call mma_deallocate(Tempj)
end do

end subroutine AddGrad
