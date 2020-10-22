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
* Copyright (C) 1991, Roland Lindh                                     *
*               1991, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Drvh2(Hess,Temp,nHess,show)
************************************************************************
*                                                                      *
* Object: driver for computation of gradient with respect to the one-  *
*         electron hamiltonian and the overlap matrix. The former will *
*         be contracted with the "variational" first order density     *
*         matrix and the latter will be contracted with the generalized*
*         Fock matrix.                                                 *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
*             Anders Bernhardsson Dept. of Theoretical Chemistry,      *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      use Basis_Info, only: nCnttp, dbsc, nBas
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External NaHss,OvrHss, KneHss,PrjHss,SROHss,M1Hss,PCMHss
      External NaMmH,OvrMmH, KneMmH,PrjMMH,sroMMH,M1MMH,PCMMMH
#include "real.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
      Character Label*80
      Real*8    Hess(nHess), Temp(nHess)
      Logical DiffOp,show, lECP
      Real*8, Allocatable:: Fock(:), D0(:), Coor(:,:)
      Integer, Allocatable:: lOper(:)
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
      Call StatusLine(' McKinley:',
     &                ' Computing 1-electron 2rd order derivatives')
*                                                                      *
************************************************************************
*                                                                      *
*     Get the variational density matrix and the occupied Fock matrix.
*
      nFock = 0
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nFock = nFock + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
*
*     Read the variational 1st order density matrix
*     density matrix in AO/SO basis
      Call mma_allocate(D0,nDens,Label='D0')
      Call Get_D1ao_Var(D0,nDens)
*     Read the generalized Fock matrix
*     Fock matrix in AO/SO basis
      Call mma_allocate(Fock,nFock,Label='Fock')
      Call Get_Fock_Occ(Fock,nFock)
*                                                                      *
************************************************************************
*                                                                      *
*...  Prologue
      nComp = 1
      Call mma_allocate(Coor,3,nComp,Label='Coor')
      Coor(:,:)=Zero
      Call mma_allocate(lOper,nComp,Label='lOper')
      lOper(:) = 1
************************************************************************
*1)                                                                    *
*     Trace the generalized Fock matrix with the gradient of the       *
*     overlap matrix.                                                  *
*                                                                      *
************************************************************************
*
      DiffOp = .False.
      Temp(:)=Zero
      Label  = ' The Renormalization Contribution'
      Call Dot1El(OvrHss,OvrMmH,Temp,nHess,DiffOp,Coor,
     &           Fock,nFock,lOper,nComp,Label)
      If (show) write(6,*) label
      If (show) Call HssPrt(Hess,nHess)
      Hess(:) = Hess(:) - Temp(:)
*
************************************************************************
*2)                                                                    *
*     Trace the "variational" zero  order density matrix with the      *
*     gradient of the kinetic energy integrals.                        *
*                                                                      *
************************************************************************
*
      DiffOp = .False.
      Temp(:)=Zero
      Label  = ' The Kinetic Energy Contribution'
      Call Dot1El(KneHss,KneMmH,Temp,nHess,DiffOp,Coor,
     &            D0,nFock,lOper,nComp,Label)
      If (show) write(6,*) label
      If (show) Call HssPrt(Temp,nHess)
      Hess(:) = Hess(:) - Temp(:)
*
************************************************************************
*3)                                                                    *
*     Trace the "variational" zero  order density matrix with the      *
*     gradient of the nuclear attraction integrals.                    *
*                                                                      *
************************************************************************
*
      DiffOp = .True.
      Label = ' The Nuclear Attraction Contribution'
      Temp(:)=Zero
      Call Dot1El(NAHss,NAMmH,Temp,nHess,DiffOp,Coor,
     &            D0,nFock,lOper,nComp,Label)
      If (show) write(6,*) label
      if (show) Call HssPrt(Temp,nHess)
      Hess(:) = Hess(:) + Temp(:)
*
************************************************************************
*3)                                                                    *
*     Trace the "variational" zero  order density matrix with the      *
*     gradient of the ECP integrals.                                   *
*                                                                      *
************************************************************************
*
      lECP = .False.
      DO i = 1, nCnttp
         lECP = lECP .or. dbsc(i)%ECP
      End Do
      If (lECP) Then
        DiffOp = .True.
        Label = ' The Projection (ECP) Contribution'
        Temp(:)=Zero
        Call Dot1El(PrjHss,PRJMMH,Temp,nHess,DiffOp,Coor,
     &               D0,nFock,lOper,nComp,Label)
        If (show) write(6,*) label
        if (show) Call HssPrt(Temp,nHess)
        Hess(:) = Hess(:) + Temp(:)
*
        DiffOp = .True.
        Label = ' The Spec. Res. (ECP) Contribution'
        Temp(:)=Zero
        Call Dot1El(SROHss,SROMMH,Temp,nHess,DiffOp,Coor,
     &               D0,nFock,lOper,nComp,Label)
        if (show) Write(6,*) Label,'first part '
        if (show) Call HssPrt(Temp,nHess)
        Hess(:) = Hess(:) + Temp(:)
*
        DiffOp = .True.
        Label = ' The M1 (ECP) Contribution'
        Temp(:)=Zero
        Call Dot1El(m1Hss,m1MMH,Temp,nHess,DiffOp,Coor,
     &               D0,nFock,lOper,nComp,Label)
        if (show) Write(6,*) Label,'second part '
        if (show) Call HssPrt(Temp,nHess)
        Hess(:) = Hess(:) + Temp(:)
      End If
*
************************************************************************
*4)                                                                    *
*     Trace the "variational" zero  order density matrix with the      *
*     gradient of the PCM integrals.                                   *
*                                                                      *
************************************************************************
*
      If (PCM) Then
        DiffOp = .True.
        Label = ' The PCM Contribution'
        Temp(:)=Zero
        Call Dot1El(PCMHss,PCMMMH,Temp,nHess,DiffOp,Coor,
     &               D0,nFock,lOper,nComp,Label)
        If (show) write(6,*) label
        if (show) Call HssPrt(Temp,nHess)
        Hess(:) = Hess(:) + Temp(:)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Epilogue, end
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(lOper)
      Call mma_deallocate(Coor)
      Call mma_deallocate(Fock)
      Call mma_deallocate(D0)
      If (Show) Call HssPrt(Hess,nHess)
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
      Call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)
      Return
      End
