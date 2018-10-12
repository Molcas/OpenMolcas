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
     $                    AtmC,IAtm,LcAtmC,LcIAtm,iOper,nIrrep,NonEq)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "status.fh"
      Real*8 AtmC(3,NAtm),LcAtmC(3,NAtm)
      Integer IAtm(NAtm),LcIAtm(NAtm)
      Logical NonEq
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
      Call GetMem('XSph','Allo','Real',ipp_Xs,MxSph)
      Call GetMem('YSph','Allo','Real',ipp_Ys,MxSph)
      Call GetMem('ZSph','Allo','Real',ipp_Zs,MxSph)
      Call GetMem('RSph','Allo','Real',ipp_R ,MxSph)
      Call GetMem('NOrd','Allo','Inte',ipp_N ,MxSph)
      Call IZero(iWork(ipp_N),MxSph)
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
      Call GetMem('NOrd'   ,'Free','Inte',ipp_N     ,MxSph)
      Call GetMem('RSph'   ,'Free','Real',ipp_R     ,MxSph)
      Call GetMem('ZSph'   ,'Free','Real',ipp_Zs    ,MxSph)
      Call GetMem('YSph'   ,'Free','Real',ipp_Ys    ,MxSph)
      Call GetMem('XSph'   ,'Free','Real',ipp_Xs    ,MxSph)
*
*---- Prepare an input file for GeomView visualization tool
*
      Call GetMem('IVTS','Allo','Inte',ip_VTS,MxVert*nTs)
      Call GVWrite(1,nTs,NSinit,LcNAtm,LcAtmC,LcIAtm,Work(ip_Sph),
     &             Work(ip_Tess),iWork(ip_NVert),
     &             Work(ip_Vert),iWork(ip_ISph),RJunk,iWork(ip_VTS),
     &             MxVert)
      Call GetMem('IVTS','Free','Inte',ip_VTS,MxVert*nTs)

*
*---- If needed compute the geometrical derivatives
*
      If(DoDeriv) then
        RSolv = RSlPar(19)
        Call Deriva(0,ToAng,LcNAtm,nTs,nS,nSInit,RSolv,
     $              Work(ip_Tess),Work(ip_Vert),Work(ip_Centr),
     $              Work(ip_Sph),iWork(ip_ISph),iWork(ip_IntS),
     $              iWork(ip_N),iWork(ip_NVert),iWork(ip_NewS),
     $              Work(ip_DTes),Work(ip_DPnt),Work(ip_DRad),
     $              Work(ip_DCntr))
      EndIf
*
*---- Compute cavitation energy
*
      TAbs = RSlPar(16)
      Call Cavitation(DoDeriv,ToAng,LcNAtm,NS,nTs,RSlPar(46),VMol,TAbs,
     &                TCE,RSolv,Work(ip_Sph),Work(ip_Tess),
     &                iWork(ip_ISph))
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
     &            Work(ip_Sph), Work(ip_Tess),
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
         Call Unused_integer(iOper)
         Call Unused_integer(nIrrep)
      End If
      End
