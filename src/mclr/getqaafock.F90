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

use MCLR_Data, only: nNA, nDens2
use input_mclr, only: nRoots, ntAsh, ntBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One

implicit none
! Input
integer nP2
real*8, dimension((nRoots-1)*nRoots/2) :: zX
real*8, dimension(nRoots*(nRoots+1)/2,nnA,nnA) :: GDMat
real*8, dimension(nP2) :: P2MOt
! Output
real*8, dimension(nDens2) :: FOccMO
! For Debugging
integer NPUVX
real*8, dimension(:), allocatable :: PUVX
integer, dimension(ntAsh,ntAsh,ntAsh,ntAsh) :: IndTUVX
integer, dimension(ntBas,ntAsh,ntAsh,ntAsh) :: IndPUVX
logical debug2
! Auxiliaries
real*8, dimension(:), allocatable :: G1r, G2r, G2q, Fock, T, PQaa
integer K, L, nG2r, IKL, IKL2, IKK, ILL, nG1, nG2, nG1r
! Statement function
integer i, j, itri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *
ng1 = itri(ntash,ntash)
ng2 = itri(ng1,ng1)

nG1r = ntash**2
nG2r = (nG1r+1)*nG1r/2
call mma_allocate(Fock,nDens2)
call mma_allocate(T,nDens2)
call mma_allocate(G1r,nG1r)
call mma_allocate(G2r,nG2r)
call mma_allocate(G2q,ng2)
call mma_allocate(PQaa,ng2)
Debug2 = .false.
call FZero(PQaa,nP2)
call FZero(G1r,nG1r)
call FZero(G2r,nG2r)

!*****************

if (Debug2) then
  call Get_PUVXLen(NPUVX)
  call mma_allocate(PUVX,NPUVX)
  call Get_Two_Ind(IndPUVX,IndTUVX)
end if

do K=1,nRoots
  IKK = (K+1)*K/2
  do L=1,K-1
    ILL = (L+1)*L/2
    IKL = (K-1)*K/2+L
    IKL2 = (K-1)*(K-2)/2+L
    call QaaP2MO(G2q,ng2,GDMat,IKL,IKK,ILL)
    if (Debug2) call QaaVerif(G2q,ng2,PUVX,NPUVX,IndTUVX)
    call G2qtoG2r(G2r,G2q,nG2,nG2r)
    call daxpy_(ng2,zX(IKL2),G2q,1,PQaa,1)
  end do
end do

call Daxpy_(nG2,One,PQaa,1,P2MOt,1)
call Put_dArray('P2MOt',P2MOt,nG2)

call G2qtoG2r(G2r,PQaa,nG2,nG2r)
call FockGen(Zero,G1r,G2r,T,Fock,1)
call DAxPy_(nDens2,One,T,1,FOccMO,1)
if (Debug2) call mma_deallocate(PUVX)
call mma_deallocate(Fock)
call mma_deallocate(T)
call mma_deallocate(G1r)
call mma_deallocate(G2r)
call mma_deallocate(G2q)
call mma_deallocate(PQaa)

end subroutine GetQaaFock
