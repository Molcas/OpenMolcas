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
! Copyright (C) 1992,1995, Per-Olof Widmark                            *
!               1992,1995, Markus P. Fuelscher                         *
!               1992,1995, Piotr Borowski                              *
!               1992,1995, Martin Schuetz                              *
!               2016,2017, Roland Lindh                                *
!***********************************************************************

subroutine MemAlo()
!***********************************************************************
!                                                                      *
!     purpose: allocate memory for density & fock matrices etc.        *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use LnkLst, only: NodSiz
use InfSCF, only: Aufb, CMO, CMO_Ref, Dens, DSCF, EDFT, EOrb, FockAO, FockMO, HDiag, MaxBas, MemRsv, mOV, MxIter, MxOptm, nBB, &
                  nBO, nBT, nD, nDens, nIter, nMem, nnB, nnOc, nOO, nOV, OccNo, OrbType, TrM, TwoHam, Vxc
use stdalloc, only: mma_allocate, mma_maxDBLE
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: lthCor, lthDii, lthGrd, lthH, lthLiS, lthPMt, lthRst, lthTot, Mx_nIter, MxMem, nIt0

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call Setup_SCF()

! Allocate memory for TrMat, CMO and occupation numbers

call mma_allocate(TrM,nBB,nD,Label='TrM')
call mma_allocate(CMO,nBB,nD,Label='CMO')
call mma_allocate(CMO_Ref,nBB,nD,Label='CMO_Ref')

call mma_allocate(FockAO,nBT,nD,Label='FockAO')
FockAO(:,:) = Zero
call mma_allocate(FockMO,nOO,nD,Label='FockMO')
FockMO(:,:) = Zero

call mma_allocate(OccNo,nnB,nD,Label='OccNo')
OccNo(:,:) = Zero
call mma_allocate(EOrb,nnB,nD,Label='EOrb')
EOrb(:,:) = Zero
call mma_allocate(OrbType,nnB,nD,Label='OrbType')
OrbType(:,:) = 0

nIt0 = 0
Mx_nIter = max(nIter(0),nIter(1)+nIt0)

! Allocate Dens and TwoHam
! a) permanently in core (apart from Dens and TwoHam)
lthCor = 3*nBT+2*nBB+2*nnB+nnOc+MxOptm+1+(MxOptm+1)**2+MxIter+MxIter**2+nTri_Elem(Mx_nIter)+1
! b) space needed by PMat
if (DSCF) then
  lthPMt = 1024*1024+2*nBT
else
  lthPMt = nBB+2*(MaxBas**2)
end if
! c) the biggest scratch space is used by SubRoutine Diis
lthgrd = nOO+nOV+2*nBT+nBT+3*MaxBas**2+nBT+nnB
lthDii = max(2*nOO,lthgrd)
lthLiS = nOV+nBO+lthgrd

lthTot = lthCor+max(lthPMt,max(lthDii,lthLiS))+MxIter*NodSiz*5
!mgs   this has to be fixed once in a more reasonable way...
!MemRsv = lthTot
MemRsv = 0
!mgs
call mma_maxDBLE(MxMem)
lthTot = lthTot+5*nOV
lthRst = MxMem-lthTot
nDens = min(lthRst/(nBT*nD)/2,6)
! We need at least 2 Dens in core at the same time for computing
! the DIIS error vectors
if (nDens < 2) then
  write(u6,*) 'MemAlo: nDens < 2'
  write(u6,*) 'lthTot=',lthTot
  write(u6,*) 'nOV=',nOV
  write(u6,*) 'MxMem=',MxMem
  write(u6,*) 'nDens=',nDens
  write(u6,*) 'lthRst=',lthRst
  write(u6,*) 'nD=',nD
  write(u6,*) 'nBT=',nBT
  call Abend()
end if
! Francesco Aquilante
! for large basis sets reserve some extra memory for the
! auxiliary matrices in egrad.f (limit set to 400 bsf)
! We need at least 2 Dens in core... but the optimal performance
! happens when we can read and store in memory 5 densities
if (nBT >= 80200) nDens = min(nDens,6)

if (nDens > Mx_nIter+1) nDens = Mx_nIter+1
if (nDens < 2) nDens = 2

nMem = nDens-1

call mma_allocate(Dens,nBT,nD,nDens,Label='Dens  ')
Dens(:,:,:) = Zero
call mma_allocate(TwoHam,nBT,nD,nDens,Label='TwoHam')
TwoHam(:,:,:) = Zero
call mma_allocate(Vxc,nBT,nD,nDens,Label='Vxc')
Vxc(:,:,:) = Zero
call mma_allocate(EDFT,MxIter,Label='EDFT')
EDFT(:) = Zero

! Allocate memory for diagonal Hessian with respect to the elements
! of the anti-symmetric matrix kappa which represents the
! rotations. Note that in the aufbau section that the size of this
! matrix is not firmly know - the number of occupied orbitals in
! each irrep is not known. For the unrestricted option the number
! of alpha and beta orbitals may vary. Here we are generous and
! allocate memory with a size that fits all.

!lthH = 0
!do iSym=1,nSym
!  nB = nBas(iSym)
!  mB = nB/2
!  lthH = lthH+(nB-mB)*mB
!end do
if (Aufb) then
  lthH = nBB*nD
else
  lthH = mOV
end if
call mma_allocate(HDiag,lthH,Label='HDiag')

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine MemAlo
