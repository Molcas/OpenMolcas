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
      Subroutine PCM_Init(iPrint,ICharg,NAtm,ToAng,
     &                    AtmC,IAtm,LcAtmC,LcIAtm,nIrrep,NonEq)
      use PCM_arrays
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
#include "status.fh"
      Real*8 AtmC(3,NAtm),LcAtmC(3,NAtm)
      Integer IAtm(NAtm),LcIAtm(NAtm)
      Logical NonEq
      Dimension RJunk(1)
      Real*8, Allocatable:: Xs(:), Ys(:), Zs(:), Rs(:)
      Integer, Allocatable:: pNs(:), VTS(:)
*
*     Build the cavity.
*     Write the input file for GeomView.
*     Form the PCM matrix.
*
*     Possibly print parameter values
      If(iPrint.ge.99) Then
        write(6,'(''PCM parameters'')')
        Do 10 I = 1, 100
   10     write(6,'(''ISlpar('',i3,'') ='',i6)')I, ISlPar(I)
        Do 20 I = 1, 100
   20     write(6,'(''RSlpar('',i3,'') ='',F8.3)')I, RSlPar(I)
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
        If(IAtm(I).gt.0) then
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
      pNs(:)=0
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
*---- Prepare an input file for GeomView visualization tool
*
      Call mma_allocate(VTS,MxVert*nTs,Label='VTS')
      Call GVWrite(1,nTs,NSinit,LcNAtm,LcAtmC,LcIAtm,PCMSph,
     &             PCMTess,iWork(ip_NVert),
     &             Vert,iWork(ip_ISph),RJunk,VTS,MxVert)
      Call mma_deallocate(VTS)

*
*---- If needed compute the geometrical derivatives
*
      If(DoDeriv) then
        RSolv = RSlPar(19)
        Call Deriva(0,ToAng,LcNAtm,nTs,nS,nSInit,RSolv,
     $              PCMTess,Vert,Centr,
     $              PCMSph,iWork(ip_ISph),iWork(ip_IntS),
     $              iWork(ip_N),iWork(ip_NVert),iWork(ip_NewS),
     $              DTes,dPnt,dRad,dCntr)
      EndIf
*
*---- Compute cavitation energy
*
      TAbs = RSlPar(16)
      Call Cavitation(DoDeriv,ToAng,LcNAtm,NS,nTs,RSlPar(46),VMol,TAbs,
     &                TCE,RSolv,PCMSph,PCMTess,iWork(ip_ISph))
*
*---- Define PCM matrix: the inverse is stored in ip_DM
*
      nTs2 = nTs * nTs
c     LenScr = 2 * nTs
      Call GetMem('SMat','Allo','Real',ip_SM,nTs2)
      Call GetMem('SDMat','Allo','Real',ip_SDM,nTs2)
      Call GetMem('TMat','Allo','Real',ip_TM,nTs2)
      Call GetMem('RMat','Allo','Real',ip_RM,nTs2)
      If (NonEq) Then
         Eps_=EpsInf
      Else
         Eps_=Eps
      End If
      Call MatPCM(nTs,Eps_,Conductor,iWork(ip_ISph),
     &            PCMSph, PCMTess,
     &            Work(ip_DM),Work(ip_SM),Work(ip_SDM),
     &            Work(ip_TM),Work(ip_RM))
      Call GetMem('RMat','Free','Real',ip_RM,nTs2)
      Call GetMem('TMat','Free','Real',ip_TM,nTs2)
      Call GetMem('SDMat','Free','Real',ip_SDM,nTs2)
      Call GetMem('SMat','Free','Real',ip_SM,nTs2)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nIrrep)
      End If
      End
