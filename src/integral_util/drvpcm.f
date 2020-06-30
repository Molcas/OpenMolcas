************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine DrvPCM(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)
      use Basis_Info
      use PCM_arrays, only: PCMTess, PCMDM
      Implicit Real*8 (A-H,O-Z)
      Real*8 h1(nh1), TwoHam(nh1), D(nh1)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "stdalloc.fh"
      Logical First, Dff, NonEq
      Real*8, Allocatable:: Cord(:,:), Chrg(:), PCM_charge(:,:),
     &                      V_Slow(:), Q_Slow(:), V_Save(:,:),
     &                      V_Tile(:,:)
*
      iRout = 1
      iPrint = nPrint(iRout)
      Call QEnter('DrvPCM')
*                                                                      *
************************************************************************
*                                                                      *
*     Get the total 1st order AO density matrix
*
*     (unused?)
      Call Get_D1ao(ipD1ao,nDens)
      If (nDens.ne.nh1) Then
         Call WarningMessage(2,'DrvPCM: nDens.ne.nh1')
         Write (6,*) nDens,nh1
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate list of all atoms
*
*     Cord: list of all atoms
*
      Call Get_nAtoms_All(MaxAto)
*
      Call mma_allocate(Cord,3,MaxAto,Label='Cord')
      Call mma_allocate(Chrg,  MaxAto,Label='Chrg')
*
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         Z = Charge(jCnttp)
         mCnt = dbsc(jCnttp)%nCntr
         If (AuxCnttp(jCnttp)) mCnt = 0
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            x1 = dbsc(jCnttp)%Coor(1,jCnt)
            y1 = dbsc(jCnttp)%Coor(2,jCnt)
            z1 = dbsc(jCnttp)%Coor(3,jCnt)
            Do i = 0, nIrrep/nStab(ndc) - 1
               iFacx=iPhase(1,iCoset(i,0,ndc))
               iFacy=iPhase(2,iCoset(i,0,ndc))
               iFacz=iPhase(3,iCoset(i,0,ndc))
               Cord(1,nc)  = x1*DBLE(iFacx)
               Cord(2,nc)  = y1*DBLE(iFacy)
               Cord(3,nc)  = z1*DBLE(iFacz)
               Chrg(nc)    = Z
               nc = nc + 1
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(PCM_Charge,2,nTS,Label='PCM_Charge')
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(V_Tile,2,nTs,Label='V_Tile')
      Call mma_allocate(V_Save,2,nTs,Label='V_Save')
      Call mma_allocate(Q_Slow,nTs,Label='Q_Slow')
      Call mma_allocate(V_Slow,nTs,Label='V_Slow')
*
      Call DrvPCM_(h1,TwoHam,D,RepNuc,nh1,First,NonEq,
     &             Chrg,Cord,MaxAto,PCMTess,PCMDM,V_Tile,
     &             V_Save,PCM_Charge,Q_Slow,V_Slow,nTs,Eps,EpsInf)
*
      Call mma_deallocate(V_Slow)
      Call mma_deallocate(Q_Slow)
      Call mma_deallocate(V_Save)
      Call mma_deallocate(V_Tile)
*                                                                      *
************************************************************************
*                                                                      *
*---- Put the current set of PCM charges on the run file.
*
      Call Put_dArray('PCM Charges',PCM_Charge,2*nTs)
      Call mma_deallocate(PCM_Charge)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Chrg)
      Call mma_deallocate(Cord)
      Call GetMem('D1ao','Free','Real',ipD1ao,nDens)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('DrvPCM')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Dff)
      End
