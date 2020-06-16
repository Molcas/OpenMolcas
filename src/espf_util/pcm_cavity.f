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
      Subroutine PCM_Cavity(iPrint,ICharg,NAtm,ToAng,AtmC,IAtm,IsAtMM,
     &                      LcAtmC,LcIAtm,JSurf)
      use PCM_arrays
      Implicit Real*8 (a-h,o-z)
#include "espf.fh"
*
#include "rctfld.fh"
#include "status.fh"
#include "stdalloc.fh"
      Real*8 AtmC(3,NAtm),LcAtmC(3,NAtm)
      Integer IAtm(NAtm),IsAtMM(NAtm),LcIAtm(NAtm)
      Save Rad_Cor,Surf_Inc
      Data Rad_Cor/2.0d0/,Surf_Inc/0.5d0/
      Real*8, Allocatable:: Xs(:), Ys(:), Zs(:), Rs(:)
      Integer, Allocatable:: pNs(:)
*
      Call qEnter('PCM_Cavity')
*
*     Build the cavity.
*
      Call PCMDef(ISlPar,RSlPar,iPrint)
      RSlPar(3) = 5.0d-1
      RSlPar(7) = 2.0d0
      RSlPar(9) = Rad_Cor + Dble(JSurf)*Surf_Inc
*
*     Possibly print parameter values
*
      If(iPrint.ge.99) Then
         write(6,'(''PCM parameters'')')
         Do 10 I = 1, 100
   10      write(6,'(''ISlpar('',i3,'') ='',i6)')I, ISlPar(I)
         Do 20 I = 1, 100
   20      write(6,'(''RSlpar('',i3,'') ='',F8.3)')I, RSlPar(I)
      EndIf
*
*---- Recover solvent data
*
      Call DataSol(ISlPar(15))
*
*     It is necessary to avoid spurious "atoms" which can cause errors
*     in UATM: then let's copy only the real atoms on local arrays
*     of coordinates and atomic numbers.
*
      LcI = 0
      Do 30 I = 1, NAtm
         If(IAtm(I).gt.0 .and. IsAtMM(I).eq.0) then
           LcI = LcI + 1
           LcAtmC(1,LcI) = AtmC(1,I)
           LcAtmC(2,LcI) = AtmC(2,I)
           LcAtmC(3,LcI) = AtmC(3,I)
           LcIAtm(LcI)   = IAtm(I)
         EndIf
   30 Continue
      LcNAtm = LcI
      ISlPar(42) = LcNAtm
*
*---- Define atomic/group spheres
*     Allocate space for X, Y, Z, Radius and NOrd for MxSph spheres
*
      Call mma_allocate(Xs,MxSph,Label='Xs')
      Call mma_allocate(Ys,MxSph,Label='Ys')
      Call mma_allocate(Zs,MxSph,Label='Zs')
      Call mma_allocate(Rs,MxSph,Label='Rs')
      Call mma_allocate(pNs,MxSph,Label='pNs')
*
      NSinit = 0
      Call FndSph(LcNAtm,ICharg,ToAng,LcAtmC,LcIAtm,ISlPar(9),
     &            ISlPar(14),RSlPar(9),Xs,Ys,Zs,Rs,pNs,iPrint)
*
*---- Define surface tesserae
*
      Call FndTess(iPrint,ToAng,LcNAtm,Xs,Ys,Zs,Rs,pNs,MxSph)
*
      Call mma_deallocate(pNs)
      Call mma_deallocate(Rs)
      Call mma_deallocate(Zs)
      Call mma_deallocate(Ys)
      Call mma_deallocate(Xs)
*
*---- If needed compute the geometrical derivatives
*
      If(DoDeriv) then
         RSolv = RSlPar(19)
         LcNAtm = ISlPar(42)
         Call mma_allocate(dTes,nTs,LcNAtm,3,Label='dTes')
         Call mma_allocate(dPnt,nTs,LcNAtm,3,3,Label='dPnt')
         Call mma_allocate(dRad,nS ,LcNAtm,3,Label='dRad')
         Call mma_allocate(dCntr,nS ,LcNAtm,3,3,Label='dCntr')
         Call mma_allocate(PCM_SQ,2,nTs,Label='PCM_SQ')
         Call Deriva(0,ToAng,LcNAtm,nTs,nS,nSInit,RSolv,
     $               PCMTess,Work(ip_Vert),Work(ip_Centr),
     $               PCMSph,iWork(ip_ISph),iWork(ip_IntS),
     $               iWork(ip_N),iWork(ip_NVert),iWork(ip_NewS),
     $               dTes,dPnt,dRad,dCntr)
         If (nPCM_info.eq.0) Then
            Write(6,'(A)') ' GEPOL failed to compute the grid deriv.'
            Write(6,'(A)') ' Reduce the number of surfaces.'
            Call Quit_OnUserError()
         End If
      EndIf
*
      Call qExit('PCM_Cavity')
      Return
      End
