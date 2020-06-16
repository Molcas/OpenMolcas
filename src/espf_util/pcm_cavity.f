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
      Call GetMem('XSph','Allo','Real',ipp_Xs,MxSph)
      Call GetMem('YSph','Allo','Real',ipp_Ys,MxSph)
      Call GetMem('ZSph','Allo','Real',ipp_Zs,MxSph)
      Call GetMem('RSph','Allo','Real',ipp_R ,MxSph)
      Call GetMem('NOrd','Allo','Inte',ipp_N ,MxSph)
*
      NSinit = 0
      Call FndSph(LcNAtm,ICharg,ToAng,LcAtmC,LcIAtm,ISlPar(9),
     &            ISlPar(14),RSlPar(9),Work(ipp_Xs),Work(ipp_Ys),
     &            Work(ipp_Zs),Work(ipp_R),iWork(ipp_N),iPrint)
*
*---- Define surface tesserae
*
      Call FndTess(iPrint,ToAng,LcNAtm,
     &             ipp_Xs,ipp_Ys,ipp_Zs,ipp_R,ipp_N)
*
      Call GetMem('NOrd','Free','Inte',ipp_N ,MxSph)
      Call GetMem('RSph','Free','Real',ipp_R ,MxSph)
      Call GetMem('ZSph','Free','Real',ipp_Zs,MxSph)
      Call GetMem('YSph','Free','Real',ipp_Ys,MxSph)
      Call GetMem('XSph','Free','Real',ipp_Xs,MxSph)
*
*---- If needed compute the geometrical derivatives
*
      If(DoDeriv) then
         RSolv = RSlPar(19)
         LcNAtm = ISlPar(42)
         nDeg=3*LcNAtm
         Call mma_allocate(dTes,nTs,LcNAtm,3,Label='dTes')
         Call GetMem('DerPunt' ,'Allo','Real',ip_DPnt ,3*nTs*NDeg)
         Call GetMem('DerRad'  ,'Allo','Real',ip_DRad ,nS*NDeg)
         Call GetMem('DerCentr','Allo','Real',ip_DCntr,3*nS*NDeg)
         Call GetMem('PCM-Q','Allo','Real',ip_Q,2*nTs)
         Call Deriva(0,ToAng,LcNAtm,nTs,nS,nSInit,RSolv,
     $               Work(ip_Tess),Work(ip_Vert),Work(ip_Centr),
     $               Work(ip_Sph),iWork(ip_ISph),iWork(ip_IntS),
     $               iWork(ip_N),iWork(ip_NVert),iWork(ip_NewS),
     $               dTes,Work(ip_DPnt),Work(ip_DRad),Work(ip_DCntr))
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
