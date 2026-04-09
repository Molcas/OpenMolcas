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

subroutine AddGrad_sp(rKappa,rMat,F,idsym,r1,r2)
! Purpose:
!   Adds the contribution from the gradient to
!    [2]
!   E   Kappa. This is done to insure us about
!   a beautiful convergence of the PCG,
!   which is just the case if E is symmetric.

use Symmetry_Info, only: Mul
use MCLR_Data, only: ipCM, ipMat, nDens, nDensC
use input_mclr, only: nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: rKappa(nDensC), F(nDens), r1, r2
real(kind=wp), intent(out) :: rMat(nDensC)
integer(kind=iwp), intent(in) :: idSym
integer(kind=iwp) :: iS, jS
real(kind=wp), allocatable :: K(:), M(:), Tempi(:), Tempj(:)

call mma_allocate(K,nDens,Label='K')
call mma_allocate(M,nDens,Label='M')
M(:) = Zero
call Unc(rkappa,K,idsym,r1)

do iS=1,nSym
  js = Mul(is,idsym)
  if (nOrb(is)*nOrb(js) == 0) cycle
  call mma_allocate(Tempi,nOrb(is)**2,Label='Tempi')
  call mma_allocate(Tempj,nOrb(js)**2,Label='Tempj')
  !  T=Brillouin matrix
  call DGeSub(F(ipCM(is)),nOrb(is),'N',F(ipCM(is)),nOrb(is),'T',Tempi,nOrb(is),nOrb(is),nOrb(is))
  call DGeSub(F(ipCM(js)),nOrb(js),'N',F(ipCM(js)),nOrb(js),'T',Tempj,nOrb(js),nOrb(js),nOrb(js))
  !        t             t
  ! { Kappa T  r2 T kappa  }
  call DGEMM_('T','N',nOrb(is),nOrb(js),nOrb(js),One,K(ipMat(js,is)),nOrb(js),Tempj,nOrb(js),Zero,M(ipMat(is,js)),nOrb(is))
  call DGEMM_('N','T',nOrb(is),nOrb(js),nOrb(is),r2,Tempi,nOrb(is),K(ipmat(js,is)),nOrb(js),One,M(ipMat(is,js)),nOrb(is))
  call mma_deallocate(Tempi)
  call mma_deallocate(Tempj)
end do
call COMPRESS(M,rmat,idsym)
call mma_deallocate(K)
call mma_deallocate(M)

end subroutine AddGrad_sp
