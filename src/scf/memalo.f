************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992,1995, Per-Olof Widmark                            *
*               1992,1995, Markus P. Fuelscher                         *
*               1992,1995, Piotr Borowski                              *
*               1992,1995, Martin Schuetz                              *
*               2016,2017, Roland Lindh                                *
************************************************************************
      SubRoutine MemAlo()
************************************************************************
*                                                                      *
*     purpose: allocate memory for density & fock metrices etc.        *
*                                                                      *
*     called from: SCF                                                 *
*                                                                      *
*     calls to: mma_allocate, FZero                                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher, P. Borowski and M.G. Schuetz       *
*     University of Lund, Sweden, 1992,95                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: split off from SOrb (Martin G. Schuetz, 1996)           *
*                                                                      *
************************************************************************
      use SCF_Arrays
      use Orb_Type
      use LnkLst
      use InfSO
      use InfSCF
      use MxDM
      use Constants
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call Setup()
*
*---- Allocate memory for TrMat, CMO and occupation numbers
*
      nD=(iUHF+1)
      Call mma_allocate(TrM,nBB,nD,Label='TrM')
      Call mma_allocate(CMO,nBB,nD,Label='CMO')
      Call mma_allocate(CMO_Ref,nBB,nD,Label='CMO_Ref')
*
      Call mma_allocate(FockAO,nBT,nD,Label='FockAO')
      FockAO(:,:)=Zero
      Call mma_allocate(FockMO,nOO,nD,Label='FockMO')
      FockMO(:,:)=Zero
*
      Call mma_allocate(OccNo,nnB,nD,Label='OccNo')
      OccNo(:,:)=Zero
      Call mma_allocate(EOrb,nnB,nD,Label='EOrb')
      EOrb(:,:)=Zero
      Call mma_allocate(OrbType,nnB,nD,Label='OrbType')
      OrbType(:,:)=0

      nIt0=0
      Mx_nIter=Max(nIter(0),nIter(1)+nIt0)
*
*---- Allocate Dens and TwoHam
*.... a) permanently in core (apart from Dens and TwoHam)
      lthCor = 3*nBT + 2*nBB + 2*nnB + nnOc +
     &         MxOptm+1 + (MxOptm+1)**2 +
     &         MxIter + MxIter**2 +
     &         Mx_nIter*(Mx_nIter+1)/2+1
*.... b) space needed by PMat
      If (DSCF) Then
         lthPMt = 1024*1024 + 2*nBT
      Else
         lthPMt = nBB + 2*(MaxBas**2)
      End If
*.... c) the biggest scratch space is used by SubRoutine Diis
      lthgrd = nOO+nOV+2*nBT+nBT+3*MaxBas**2+nBT+nnB
      lthDii = Max(2*nOO,lthgrd)
      lthLiS = nOV+nBO+lthgrd
*
      lthTot = lthCor + Max(lthPMt,Max(lthDii,lthLiS))
     &       + MxIter*NodSiz*5
cmgs   this has to be fixed once in a more reasonable way...
c     MemRsv = lthTot
      MemRsv = 0
cmgs
      Call mma_maxDBLE(MxMem)
      lthTot = lthTot + 5*nOV
      lthRst = MxMem - lthTot
      nDens  = Min(lthRst/(nBT*nD)/2,6)
C: We need at least 2 Dens in core at the same time for computing
C: the DIIS error vectors
      If (nDens.lt.2) Then
         Write (6,*) 'MemAlo: nDens.lt.2'
         Write (6,*) 'lthTot=',lthTot
         Write (6,*) 'nOV=',nOV
         Write (6,*) 'MxMem=',MxMem
         Write (6,*) 'nDens=',nDens
         Write (6,*) 'lthRst=',lthRst
         Write (6,*) 'nD=',nD
         Write (6,*) 'nBT=',nBT
         Call Abend()
      End If
C: Francesco Aquilante
C: for large basis sets reserve some extra memory for the
C: auxiliary matrices in egrad.f (limit set to 400 bsf)
C: We need at least 2 Dens in core... but the optimal performance
C: happens when we can read and store in memory 5 densities
      if (nBT.ge.80200)then
         nDens = Min(nDens,6)
      endif

      If (nDens.gt.Mx_nIter + 1) nDens = Mx_nIter + 1
      If (nDens.lt.2) nDens = 2
*
      nMem  = nDens - 1
*
      Call mma_allocate(Dens  ,nBT,nD,nDens,Label='Dens  ')
      Call FZero(Dens,nBT*nD*nDens)
      Call mma_allocate(TwoHam,nBT,nD,nDens,Label='TwoHam')
      Call FZero(TwoHam,nBT*nD*nDens)
      Call mma_allocate(Vxc,nBT,nD,nDens,Label='Vxc')
      Call FZero(Vxc,nBT*nD*nDens)
      Call mma_allocate(EDFT,MxIter,Label='EDFT')
      Call FZero(EDFT,MxIter)

*---  Allocate memory for diagonal Hessian with respect to the elements
*     of the anti-symmetric matrix kappa which represents the
*     rotations. Note that in the aufbau section that the size of this
*     matrix is not firmly know - the number of occupied orbitals in
*     each irrep is not known. For the unrestricted option the number
*     of alpha and beta orbitals may vary. Here we are generous and
*     allocate memory with a size that fits all.
*
*     lthH = 0
*     Do iSym = 1, nSym
*        nB = nBas(iSym)
*        mB = nB/2
*        lthH = lthH + (nB-mB)*mB
*     End Do
      If (Aufb) Then
         lthH = nBB*nD
      Else
         lthH = mOV
      End If
      Call mma_allocate(HDiag,lthH,Label='HDiag')
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
