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
* Copyright (C) 2015,2016, Ignacio Fdez. Galvan                        *
************************************************************************
      Subroutine Process_Gradients()
      use Slapaf_Info, only: Gx, Gx0, NAC, Energy, Energy0
      Implicit None
#include "WrkSpc.fh"
#include "real.fh"
#include "info_slapaf.fh"
#include "nadc.fh"
#include "stdalloc.fh"
      Logical Found
      Integer i,nRoots,RC,Read_Grad,Columbus
      Real*8, Allocatable :: Grads(:,:), Ener(:)
      Real*8 E0, E1
      External Read_Grad
*
      Request_Alaska=.False.
      Call mma_Allocate(Grads,3*nsAtom,3)
*
*     First check that all the needed gradients are available
*     and stop to compute them if they are not.
*     For a two-RunFile job, this behaves as if it's one state.
*
      If (TwoRunFiles) Then
        iState(1)=0
        iState(2)=0
      End If
*
      IF (iState(1).NE.0) iState(1)=RootMap(iState(1))
      IF (iState(2).NE.0) iState(2)=RootMap(iState(2))
      i=MAX(iState(1),iState(2))
      iState(2)=MIN(iState(1),iState(2))
      iState(1)=i
*
*     Two states
*
      If ((iState(1).ne.0).and.(iState(2).ne.0)) Then
        Do i=2,1,-1
          RC=Read_Grad(Grads(1,i),3*nsAtom,iState(i),0,0)
          If (RC.eq.0) Then
            Request_Alaska=.True.
            Call Put_iScalar('Relax CASSCF root',iState(i))
            Call Put_iScalar('NumGradRoot',iState(i))
            iState(1)=iState(i)
            iState(2)=0
            Exit
          End If
        End Do
        If ((.Not.Request_Alaska).And.NADC) Then
          RC=Read_Grad(Grads(1,3),3*nsAtom,0,iState(1),iState(2))
          If (RC.eq.0) Request_Alaska=.True.
        End If
*
*     One state
*
      Else
        iState(1)=0
        iState(2)=0
        Call Qpg_iScalar('Relax CASSCF root',Exist)
        If (Exist) Call Get_iScalar('Relax CASSCF root',iState(1))
        If (iState(1).eq.0) iState(1)=1
        RC=Read_Grad(Grads(1,1),3*nsAtom,iState(1),0,0)
        If (RC.eq.0) Request_Alaska=.True.
      End If
*
      If (Request_Alaska) Then
        Call mma_Deallocate(Grads)
        NADC=.False.
        Return
      End If
*
*     Once the gradients are read, convert them to forces
*     and store them as appropriate.
*
*     Energy and gradient of the first (higher) state
*
      nRoots=1
      Call Qpg_iScalar('Number of roots',Found)
      If (Found) Call Get_iScalar('Number of roots',nRoots)
      Call mma_Allocate(Ener,nRoots)
      Call Get_dArray('Last energies',Ener,nRoots)
      If (Max(iState(1),iState(2)).gt.nRoots) Then
        Call WarningMessage(2,'Too few energies in RUNFILE')
        Call Abend()
      End If
      Energy(iter)=Ener(iState(1))
      E1=Ener(iState(1))
      Call dCopy_(3*nsAtom,Grads(1,1),1,Gx(1,1,iter),1)
      Gx(:,:,iter) = -Gx(:,:,iter)
*
*     For a two-RunFile job, read the second (lower) energy
*     and gradient from RUNFILE2
*
      If (TwoRunFiles) Then
        Call NameRun('RUNFILE2')
        iState(2)=0
        Call Qpg_iScalar('Relax CASSCF root',Exist)
        If (Exist) Call Get_iScalar('Relax CASSCF root',iState(2))
        If (iState(1).eq.0) iState(2)=1
        nRoots=1
        Call Qpg_iScalar('Number of roots',Found)
        If (Found) Call Get_iScalar('Number of roots',nRoots)
        Call mma_Deallocate(Ener)
        Call mma_Allocate(Ener,nRoots)
        Call Get_dArray('Last energies',Ener,nRoots)
        Call Get_dArray('GRAD',Grads(1,2),3*nsAtom)
        Call NameRun('RUNFILE')
        RC=-1
      End If
*
      If (iState(2).gt.0) Then
        E0=Ener(iState(2))
*
*       In case of a true conical intersection the Lagrangian is different!
*       In that case we will have that the energy is the average energy
*       and that the constraint is the squared energy difference. We
*       change the energies and gradients here on-the-fly.
*
        If (NADC) Then
          Energy (iter)=(E1+E0)*Half
          Energy0(iter)=E1-E0
          Call daXpY_(3*nsAtom,-One,Grads(1,2),1,Gx(1,1,iter),1)
          Gx(:,:,iter) = Half * Gx(:,:,iter)
          Call dCopy_(3*nsAtom,Grads(1,2),1,Gx0(1,1,iter),1)
          Call daXpY_(3*nsAtom,-One,Grads(1,1),1,Gx0(1,1,iter),1)
          Call Get_iScalar('Columbus',Columbus)
          If (Columbus.ne.1) Then
            Call mma_allocate(NAC,3,nsAtom,Label='NAC')
            Call dCopy_(3*nsAtom,Grads(1,3),1,NAC,1)
*
*           If the coupling derivative vector could not be calculated,
*           use an approximate one.
*
            If (RC.lt.0) Then
              ApproxNADC=.True.
              Call Branching_Plane_Update(Gx,Gx0,NAC,3*nsAtom,iter)
            End If
          End If
        Else
          Energy0(iter)=E0
          Call dCopy_(3*nsAtom,Grads(1,2),1,Gx0(1,1,iter),1)
          Gx0(1,1,iter) = -Gx0(1,1,iter)
        End If
      End If
*
      Call mma_Deallocate(Ener)
      Call mma_Deallocate(Grads)
      Return
*
      End Subroutine Process_Gradients
