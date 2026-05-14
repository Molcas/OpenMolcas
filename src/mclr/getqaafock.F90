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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine GetQaaFock(FOccMO,P2MOt,GDMat,zX,nP2)

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: nDens, nNA
use input_mclr, only: nRoots, ntAsh, ntBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nP2
real(kind=wp), intent(inout) :: FOccMO(nDens), P2MOt(nP2)
real(kind=wp), intent(in) :: GDMat(nTri_Elem(nRoots),nnA,nnA), zX(nTri_Elem(nRoots-1))
integer(kind=iwp) :: IKK, IKL, IKL2, ILL, K, L, nG1, nG1r, nG2, nG2r, NPUVX
logical(kind=iwp) :: debug2
integer(kind=iwp), allocatable :: IndPUVX(:,:,:,:), IndTUVX(:,:,:,:)
real(kind=wp), allocatable :: Fock(:), G1r(:), G2q(:), G2r(:), PQaa(:), PUVX(:), T(:)

!                                                                      *
!***********************************************************************
!                                                                      *
ng1 = nTri_Elem(ntash)
ng2 = nTri_Elem(ng1)

nG1r = ntash**2
nG2r = nTri_Elem(nG1r)
call mma_allocate(Fock,nDens)
call mma_allocate(T,nDens)
call mma_allocate(G1r,nG1r)
call mma_allocate(G2r,nG2r)
call mma_allocate(G2q,ng2)
call mma_allocate(PQaa,ng2)
Debug2 = .false.
PQaa(1:nP2) = Zero
G1r(:) = Zero
G2r(:) = Zero

!*****************

if (Debug2) then
  call Get_PUVXLen(NPUVX)
  call mma_allocate(PUVX,NPUVX)
  call mma_allocate(IndTUVX,ntAsh,ntAsh,ntAsh,ntAsh)
  call mma_allocate(IndPUVX,ntBas,ntAsh,ntAsh,ntAsh)
  call Get_Two_Ind(IndPUVX,IndTUVX)
end if

do K=1,nRoots
  IKK = nTri_Elem(K)
  do L=1,K-1
    ILL = nTri_Elem(L)
    IKL = iTri(K,L)
    IKL2 = nTri_Elem(K-2)+L
    call QaaP2MO(G2q,ng2,GDMat,IKL,IKK,ILL)
    if (Debug2) call QaaVerif(G2q,ng2,PUVX,NPUVX,IndTUVX)
    call G2qtoG2r(G2r,G2q,nG2,nG2r)
    PQaa(:) = PQaa(:)+zX(IKL2)*G2q(:)
  end do
end do

P2MOt(1:nG2) = P2MOt(1:nG2)+PQaa(:)
call Put_dArray('P2MOt',P2MOt,nG2)

call G2qtoG2r(G2r,PQaa,nG2,nG2r)
call FockGen(Zero,G1r,G2r,T,Fock,1)
FOccMO(:) = FOccMO(:)+T(:)
if (Debug2) then
  call mma_deallocate(PUVX)
  call mma_deallocate(IndTUVX)
  call mma_deallocate(IndPUVX)
end if
call mma_deallocate(Fock)
call mma_deallocate(T)
call mma_deallocate(G1r)
call mma_deallocate(G2r)
call mma_deallocate(G2q)
call mma_deallocate(PQaa)

end subroutine GetQaaFock
