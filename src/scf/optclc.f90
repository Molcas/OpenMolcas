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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2016,2017, Roland Lindh                                *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine OptClc(CInter,nCI,nD,Ind,nInd)
!***********************************************************************
!                                                                      *
! purpose: calculate optimal density matrix and two-electron hamil-    *
!          tonial as well as optimal difference from interpolation     *
!          or extrapolation coefficients.                              *
!                                                                      *
! input:                                                               *
!   Dens    : density matrices (nDT,NumDT)                             *
!   TwoHam  : two-electron hamiltonian matrices (nDT,NumDT)            *
!   Vxc     : external potential       matrices (nDT,NumDT)            *
!   CInter  : interpolation coefficients (nCI)                         *
!                                                                      *
! output:                                                              *
!   Dens and TwoHam                                                    *
!                                                                      *
!***********************************************************************
      use InfSCF, only: kOptim, nBT, nDens, iDisk, MapDns
      use SCF_Arrays, only: Dens, TwoHam, Vxc
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nCI, nD, nInd
      Real*8 CInter(nCI,nD)
      Integer Ind(nInd)

      Real*8, Dimension(:,:), Allocatable:: DnsTmp, TwoTmp, VxcTmp
      Integer Iter_D, iMap, iD, i, MatNO
      Real*8  C
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!                                                                      *
! Allocate memory for matrices that contribute to the optimal one
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
      If (kOptim==1) Return
      Call mma_allocate(DnsTmp,nBT,nD,Label='DnsTmp')
      Call mma_allocate(TwoTmp,nBT,nD,Label='TwoTmp')
      Call mma_allocate(VxcTmp,nBT,nD,Label='VxcTmp')
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Accumulate liner combinations in position iPsLst
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Start with the last iteration
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
      iter_d=Ind(kOptim)
!
      iMap=MapDns(iter_d)
      If (iMap.lt.0) Then
         Call RWDTG(-iMap,DnsTmp,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
         Call RWDTG(-iMap,TwoTmp,nBT*nD,'R','TWOHAM',iDisk,SIZE(iDisk,1))
         Call RWDTG(-iMap,VxcTmp,nBT*nD,'R','dVxcdR',iDisk,SIZE(iDisk,1))
      Else
         Call DCopy_(nBT*nD,Dens  (1,1,iMap),1,DnsTmp,1)
         Call DCopy_(nBT*nD,TwoHam(1,1,iMap),1,TwoTmp,1)
         Call DCopy_(nBT*nD,Vxc   (1,1,iMap),1,VxcTmp,1)
      End If
!
      Do iD = 1, nD
!
         C = CInter(kOptim,iD)
         Call DScal_(nBT,C,DnsTmp(1,iD),1)
         Call DScal_(nBT,C,TwoTmp(1,iD),1)
         Call DScal_(nBT,C,VxcTmp(1,iD),1)
!
      End Do
!
      Call dcopy_(nBT*nD,DnsTmp,1,Dens  (1,1,nDens),1)
      Call dcopy_(nBT*nD,TwoTmp,1,TwoHam(1,1,nDens),1)
      Call dcopy_(nBT*nD,VxcTmp,1,Vxc   (1,1,nDens),1)
!
      Do i = 1, kOptim - 1
         C = CInter(i,1)
         MatNo = Ind(i)
!
         iMap=MapDns(MatNo)
         If (iMap.lt.0) Then
            Call RWDTG(-iMap,DnsTmp,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
            Call RWDTG(-iMap,TwoTmp,nBT*nD,'R','TWOHAM',iDisk,SIZE(iDisk,1))
            Call RWDTG(-iMap,VxcTmp,nBT*nD,'R','dVxcdR',iDisk,SIZE(iDisk,1))
         Else
            Call DCopy_(nBT*nD,Dens  (1,1,iMap),1,DnsTmp,1)
            Call DCopy_(nBT*nD,TwoHam(1,1,iMap),1,TwoTmp,1)
            Call DCopy_(nBT*nD,Vxc   (1,1,iMap),1,VxcTmp,1)
         End If
!
         Do iD = 1, nD
            C = CInter(i,iD)
            call daxpy_(nBT,C,DnsTmp(1,iD),1,Dens  (1,iD,nDens),1)
            call daxpy_(nBT,C,TwoTmp(1,iD),1,TwoHam(1,iD,nDens),1)
            call daxpy_(nBT,C,VxcTmp(1,iD),1,Vxc   (1,iD,nDens),1)
         End Do
!
      End Do ! i
!
! Deallocate memory
!
      Call mma_deallocate(DnsTmp)
      Call mma_deallocate(TwoTmp)
      Call mma_deallocate(VxcTmp)
!
#ifdef _DEBUGPRINT_
      Call NrmClc(Dens  (1,1,nDens),nBT*nD,'OptClc','D in iPsLst. ')
      Call NrmClc(TwoHam(1,1,nDens),nBT*nD,'OptClc','T in iPsLst. ')
      Call NrmClc(Vxc   (1,1,nDens),nBT*nD,'OptClc','V in iPsLst. ')
#endif
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
      Return
      End SubRoutine OptClc
