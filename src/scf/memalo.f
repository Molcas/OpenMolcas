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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "addr.fh"
#include "stdalloc.fh"
#include "lnklst.fh"
#include "infso.fh"
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
*
*     Produce pointers to work for codes that still use GetMem and
*     (i)Work.
*
      mAdCMO    = ip_of_Work(CMO(1,1))
      mAdCMO_ab = mAdCMO + (nD-1)*nBB
*
      Call mma_allocate(Fock,nBT,nD,Label='Fock')
      Call FZero(Fock,nBT*nD)
*
      Call mma_allocate(OccNo,nnB,nD,Label='OccNo')
      Call FZero(OccNo,nnB*nD)
      Call mma_allocate(EOrb,nnB,nD,Label='EOrb')
      Call FZero(EOrb,nnB*nD)
      Call mma_allocate(OrbType,nnB,nD,Label='OrbType')
      Call ICopy(nnB*nD,[0],0,OrbType,1)

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
      Call GetMem('SCF','Max','Real',iDum,MxMem)
      lthTot = lthTot + 5*nOV
      lthRst = MxMem - lthTot
      nDens  = Min(lthRst/(nBT*nD)/2,6)
C: We need at least 2 Dens in core at the same time for computing
C: the DIIS error vectors
      If (nDens.lt.2) Then
         Write (6,*) 'MemAlo: nDens.lt.2'
         Write (6,*) 'nDens=',nDens
         Call QTrace
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
         lthH = nBB
      Else
         lthH = nOV
      End If
      Call mma_allocate(HDiag,lthH,nD,Label='HDiag')
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
