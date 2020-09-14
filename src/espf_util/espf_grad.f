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
      Subroutine espf_grad (natom,nAtQM,nGrdPt,ipExt,ipGrid,ipB,ipDB,
     &                      ipIsMM,ipGradCl,DoTinker,DoGromacs)
      Implicit Real*8 (A-H,O-Z)
*
*     Gradient due to the external potential
*
#include "espf.fh"
*
#include "disp.fh"
#include "nac.fh"
      Logical Exist,lMMHess,lMMGrd,DoTinker,DoGromacs,isNAC_tmp
      Character*180 Line
      Character*180 Get_Ln
      External Get_Ln
      Dimension FX(4)
      Dimension opnuc(1)
*
      Call QEnter('espf_grad')
      iPL = iPL_espf()
*
      Call Get_Grad(ipGrad,nGrad)
      lMMHess = .False.
      If (nGrad.ne.3*natom) Then
         Write (6,*)
         Write (6,'(/,A)')'nGrad.ne.natom'
         Write (6,'(A,I6)')'nGrad=',nGrad
         Write (6,'(A,I6)')'natom=',natom
         Call Quit_OnUserError()
      End If
      If (iPL.ge.3) Call PrGrad(' Molecular gradients, entering ESPF',
     &                          Work(ipGrad),lDisp(0),ChDisp,4)
*
*     Recover MM gradient and hessian, if any, in QMMM file
*
      Call F_Inquire('QMMM',Exist)
      natMM = 0
      If (((Exist.and.DoTinker).or.DoGromacs).and. .not.isNAC) Then
         Call Allocate_Work(ipMMGrd,6*natom)
         Call Qpg_dArray('MM Grad',lMMGrd,6*natom)
         If (lMMGrd) Then
            Call Get_dArray('MM Grad',Work(ipMMGrd),6*natom)
            call dcopy_(3*natom,Work(ipMMGrd+3*natom),1,Work(ipMMGrd),1)
            call dcopy_(3*natom,[Zero],0,Work(ipMMGrd+3*natom),1)
         Else
            call dcopy_(6*natom,[Zero],0,Work(ipMMGrd),1)
         End If
      End If
*
      If (Exist .and. DoTinker .and. .not.isNAC) Then
         Call GetMem('Hess','Allo','Real',ipHC,3*natom*(3*natom+1)/2)
         call dcopy_((3*natom*(3*natom+1)/2),[Zero],0,Work(ipHC),1)
      End If
*
      If (Exist .and. DoTinker .and. .not.isNAC) Then
         ITkQMMM=IsFreeUnit(15)
         Call Molcas_Open(ITkQMMM,'QMMM')
         Line = ' '
         Do While (Index(Line,'TheEnd ') .eq. 0)
            Line=Get_Ln(ITkQMMM)
            If (Index(Line,'NMM').ne.0) Then
               Call Get_I1(2,natMM)
            Else If (Index(Line,'MMGradient').ne.0) Then
               Call Get_I1(2,iAtom)
               Call Get_F(3,FX,3)
               Do iXYZ = 0, 2
                  iOff = (iAtom-1)*3+iXYZ
                  Work(ipMMGrd+3*natom+iOff) = FX(iXYZ+1) * Angstrom
     &                                                    * ToHartree
                  Work(ipGrad+iOff) = Work(ipGrad+iOff)
     &                              + Work(ipMMGrd+3*natom+iOff)
               EndDo
            Else If (Index(Line,'MMHDiag').ne.0) Then
               lMMHess = .True.
               Call Get_I1(2,iAtom)
               Call Get_F(3,FX,3)
               Do iXYZ = 1, 3
                  Work(ipHC+LHR(iXYZ,iAtom,iXYZ,iAtom)-1) = FX(iXYZ)
     &                              * Angstrom * Angstrom * ToHartree
               EndDo
            Else If (Index(Line,'MMHOff').ne.0) Then
               lMMHess = .True.
               Call Get_I1(2,iAtom)
               Call Get_I1(3,iXYZ)
               Call Get_I1(4,iNumb)
c            Write (6,*) 'HOff read ',iAtom,iXYZ,iNumb
               iCur = 0
               jAtom = iAtom
               jXYZ = iXYZ
51             iStep = Min(4,iNumb-iCur)
               Line = Get_Ln(ITkQMMM)
               Call Get_F(1,FX,iStep)
c            Write (6,'(A,4f10.5)') 'HOff read ',(FX(j),j=1,iStep)
               Do iBla = 1, iStep
                  jXYZ = jXYZ + 1
                  If (jXYZ.eq.4) Then
                     jXYZ = 1
                     jAtom = jAtom + 1
                     If (jAtom.gt.natom) Call Abend
                  End If
                  Work(ipHC+LHR(iXYZ,iAtom,jXYZ,jAtom)-1) = FX(iBla)
     &                             * Angstrom * Angstrom * ToHartree
               End Do
               iCur = iCur + 4
               If (iCur.lt.iNumb) goto 51
            End If
         End Do
         Close (ITkQMMM)
      End If
*
      If (DoGromacs .and. .not.isNAC) Then
         Call dcopy_(3*natom,Work(ipGradCl),1,Work(ipMMGrd+3*natom),1)
         Do i = 0,3*natom-1
            Work(ipGrad+i) = Work(ipGrad+i) + Work(ipMMGrd+3*natom+i)
         End Do
      End If
*
      If (((Exist.and.DoTinker).or.DoGromacs) .and. .not.isNAC) Then
         If (iPL.ge.4) Then
            Call RecPrt('Old MM Grad:',' ',Work(ipMMGrd),3,natom)
            Call RecPrt('New MM Grad:',' ',Work(ipMMGrd+3*natom),3,
     &                  natom)
         End If
         If (natMM.gt.0) Call Put_dArray('MM Grad',Work(ipMMGrd),
     &                                   6*natom)
         Call Free_Work(ipMMGrd)
      End If
*
      If (Exist.and.DoTinker .and. .not.isNAC) Then
         If (lMMHess .and. iPL.ge.4)
     &      Call TriPrt(' In ESPF_grad: MM Hessian','(12f12.7)',
     &      Work(ipHC),3*natom)
         If (lMMHess) Call Put_dArray('MMHessian',Work(ipHC),
     &                                3*natom*(3*natom+1)/2)
         Call GetMem('Hess','Free','Real',ipHC,3*natom*(3*natom+1)/2)
      End If
*
      If (((Exist.and.DoTinker).or.DoGromacs) .and. .not.isNAC) Then
         Call Put_iScalar('No of Internal coordinates',3*natom)
         If (iPL.ge.3) Call PrGrad(' Molecular gradients, after MM',
     &                          Work(ipGrad),lDisp(0),ChDisp,4)
      End If
*
*     External field acting on nuclear charges
*
      If (isNac) Then
         Write(6,*) 'ESPF: Skipping nuclear-external field contribution'
      Else
         Call GetMem('XCharge','Allo','Real',ipXC,nAtom)
         Call Get_dArray('Effective nuclear Charge',Work(ipXC),nAtom)
         Do iAt = 1, nAtom
            iCurXC = ipXC+iAt-1
            iCurG = ipGrad+(iAt-1)*3
            iCurE = ipExt+(iAt-1)*MxExtPotComp
            Work(iCurG  ) = Work(iCurG  ) + Work(iCurXC)*Work(iCurE+1)
            Work(iCurG+1) = Work(iCurG+1) + Work(iCurXC)*Work(iCurE+2)
            Work(iCurG+2) = Work(iCurG+2) + Work(iCurXC)*Work(iCurE+3)
         EndDo
         Call GetMem('XCharge','Free','Real',ipXC,natom)
         If (iPL.ge.3) Call PrGrad(' Molecular grad, after nuc ESPF',
     &                          Work(ipGrad),lDisp(0),ChDisp,4)
      End If
*
*     Here I need the integral derivatives, weighted by B and contracted
*     with the density matrix: B * Sum_mu,nu P_mu,nu*d/dq(<mu|1/R_grid|nu>)
*
*     Basically, it should work like the following lines ... but it can't now.
*
*      opnuc = Dum
*      ncmp = 3
*      iAddPot = -1
*      Call GetMem('dESPF1','Allo','Real',ipD1,nGrdPt*natom*3)
*      Call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipD1),nGrdPt,iAddPot)
*      Call GetMem('dESPF1','Free','Real',ipD1,nGrdPt*natom*3)
*
*     Now I try another thing, following the way derivatives with respect to point
*     charges are computed in alaska. I just copied 2 files from alaska: drvh1 and
*     xfdgrd then renamed them drvespf and bdvgrd.

      Call GetMem('Temp','Allo','Real',ipTemp,3*natom)
      Call GetMem('GridInfo','Allo','Real',ipGrdI,4*nGrdPt)
*     Need to save isNAC here because Prepare calling inisew calling init_seward which resets isNAC .to False.
      isNAC_tmp = isNAC
      Call Prepare(nGrdPt,ipGrid,ipB,ipGrdI)
      Call Drvespf(Work(ipGrad),Work(ipTemp),3*natom,Work(ipGrdI))
      Call GetMem('GridInfo','Free','Real',ipGrdI,4*nGrdPt)
      If (iPL.ge.3) Call PrGrad(' Molecular gradients, after P*B*dV',
     &                          Work(ipGrad),lDisp(0),ChDisp,4)
      Call GetMem('Temp','Free','Real',ipTemp,3*natom)
*
*     Here I need the integrals contracted with the density matrix and weighted
*     by the derivatives of B: dB/dq * Sum_mu,nu P_mu,nu*(<mu|1/R_grid|nu>)
*
      opnuc = Dum
      ncmp = 1
      iAddPot = -1
      Call GetMem('dESPF2','Allo','Real',ipD2,nGrdPt)
      Call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipD2),nGrdPt,iAddPot)
      If (iPL.ge.4) Then
         Write(6,'(/,A,/)') ' PV = '
         Do iPnt = 1, nGrdPt
            Write(6,*) Work(ipD2+iPnt-1)
         EndDo
      End If
      iQM = 0
      Do iAt = 1, natom
         If (iWork(ipIsMM+iAt-1).eq.1) Goto 10
         iQM = iQM + 1
         Do jPnt = 1, NGrdPt
            iCurG  = ipGrad + (iAt-1)*3
            iCurDB = ipDB + (jPnt-1)*nAtQM*3+(iQM-1)*3
            iCurDB1 = ipDB + (jPnt-1)+((iQM-1)*3+0)*nGrdPt
            iCurDB2 = ipDB + (jPnt-1)+((iQM-1)*3+1)*nGrdPt
            iCurDB3 = ipDB + (jPnt-1)+((iQM-1)*3+2)*nGrdPt
            iCurI  = ipD2 + jPnt-1
            Work(iCurG  ) = Work(iCurG  ) + Work(iCurDB1)*Work(iCurI)
            Work(iCurG+1) = Work(iCurG+1) + Work(iCurDB2)*Work(iCurI)
            Work(iCurG+2) = Work(iCurG+2) + Work(iCurDB3)*Work(iCurI)
         EndDo
10       Continue
      EndDo
      isNAC = isNAC_tmp
*
*     Apply Morokuma's scheme if needed
*
      If ((Exist.and.DoTinker).or.DoGromacs)
     &   Call LA_Morok(natom,ipGrad,1)
*
*     Finally
*
      Call Put_Grad(Work(ipGrad),3*natom)
      Call GetMem('dESPF2','Free','Real',ipD2,nGrdPt)
      If (iPL.ge.2) Call PrGrad(' Molecular gradients, after ESPF',
     &                          Work(ipGrad),lDisp(0),ChDisp,4)
      Call Add_Info('Grad',Work(ipGrad),3*natom,6)
      Call GetMem('Grad','Free','Real',ipGrad,3*natom)
*
      Call QExit('espf_grad')
      Return
      End
