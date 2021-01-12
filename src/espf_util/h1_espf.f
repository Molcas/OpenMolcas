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
      Subroutine h1_espf (h1,RepNuc,nh1,First,Do_DFT)
      use Basis_Info, only: nBas
      Implicit Real*8 (A-H,O-Z)
*
* Driver for computing the ESPF one-electron modification
* of the core hamiltonian
*
#include "espf.fh"
*
#include "print.fh"
      Character*180 Line,Get_Ln
      External Get_Ln
      Character*10 ESPFKey
      Real*8 h1(nh1)
      Logical StandAlone,First,Do_DFT
      Logical DynExtPot,Exist,DoTinker,DoGromacs,lMorok,UpdateVMM
      Logical DoDirect
*
      iPL = iPrintLevel(-1)
*
      RealDummy = Zero
*
* Recover some ESPF data
*
      StandAlone = .False.
      DoTinker = .False.
      DoGromacs = .False.
      lMorok = .False.
      DoDirect = .False.
      ipOldMltp = ip_Dummy
      Call F_Inquire('ESPF.DATA',Exist)
      If (Exist) Then
         IPotFl = IsFreeUnit(1)
         Call Molcas_Open(IPotFl,'ESPF.DATA')
10       Line = Get_Ln(IPotFl)
         ESPFKey = Line(1:10)
         If (ESPFKey.eq.'MLTORD    ') Then
            Call Get_I1(2,MltOrd)
            ibla = 0
            Do ii = 0, MltOrd
               ibla = ibla + (ii+2)*(ii+1)/2
            End Do
            MltOrd = ibla
         Else If (ESPFKey.eq.'IRMAX     ') Then
            Call Get_I1(2,iRMax)
         Else If (ESPFKey.eq.'DELTAR    ') Then
            Call Get_F1(2,DeltaR)
         Else If (ESPFKey.eq.'GRIDTYPE  ') Then
            Call Get_I1(2,iGrdTyp)
         Else If (ESPFKey.eq.'TINKER    ') Then
            DoTinker = .True.
         Else If (ESPFKey.eq.'GROMACS   ') Then
            DoGromacs = .True.
         Else If (ESPFKey.eq.'LA_MOROK  ') Then
            lMorok = .True.
         Else If (ESPFKey.eq.'DIRECT    ') Then
            DoDirect = .True.
         Else If (ESPFKey.eq.'MULTIPOLE ') Then
            Call Get_I1(2,nMult)
            Call Allocate_Work(ipOldMltp,nMult)
            Do iMlt = 1, nMult, MltOrd
               Line = Get_Ln(IPotFl)
               Call Get_I1(1,iAt)
               Call Get_F(2,Work(ipOldMltp+iMlt-1),MltOrd)
            End Do
         Else If (ESPFKey.eq.'ENDOFESPF ') Then
            Goto 11
         EndIf
         Goto 10
11       Close (IPotFl)
      Else
         Write(6,*) 'No ESPF.DATA file. Abort'
         Call Quit_OnUserError()
      End If
*
* Is this a fully coupled QM/MM calculation?
*
      DynExtPot = .False.
      iMode = 0
      If (DoTinker) Then
         ITkQMMM = IsFreeUnit(30)
         Call Molcas_Open (ITkQMMM,'QMMM')
         Line = ' '
         Do While (Index(Line,'TheEnd ') .eq. 0)
            Line=Get_Ln(ITkQMMM)
            If (Index(Line,'FullCoupling').ne.0) Then
               DynExtPot = .True.
               iMode = Max(iMode,1)
            Else If (Index(Line,'MMPolar').ne.0) Then
               DynExtPot = .True.
               iMode = Max(iMode,2)
            End If
         End Do
         Close (ITkQMMM)
      End If
      If(.not. DynExtPot) Then
        If (ipOldMltp .ne. ip_Dummy) Call Free_Work(ipOldMltp)
        Return
      End If
*
      ipIsMM = ip_iDummy
      Call espf_init(natom,nAtQM,ipCord,ipIsMM,ipExt)
      nMult = MltOrd * nAtQM
*
* Compute the grid around the molecule
*
      nGrdPt = 0
      ipGrid = ip_Dummy
      ipDGrd = ip_Dummy
      Call StatusLine(' espf:',' Making the grid')
      If (iGrdTyp.eq.1) Then
         Call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.False.,
     &      ipIsMM,-iGrdTyp,ipDGrd,nAtQM)
         Call GetMem('ESPF_Grid','Allo','Real',ipGrid,3*nGrdPt)
         Call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.False.,
     &     ipIsMM,iGrdTyp,ipDGrd,nAtQM)
      Else
         Call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.False.,
     &      ipIsMM,iGrdTyp,ipDGrd,nAtQM)
      End If
*
* Compute the cartesian tensor T, TtT^-1, [TtT^-1]Tt
* and B=ExtPot[TtT^-1]Tt
* Tt means the transpose of T
*
* Warning, at this point ExtPot is not filled (only TTT is needed)
*
      iSize1 = nMult * nGrdPt
      iSize2 = nMult * nMult
      iSize3 = nMult * Max(nMult,nGrdPt)
      Call GetMem('CartTensor','Allo','Real',ipT,iSize1)
      Call GetMem('TT','Allo','Real',ipTT,iSize2)
      Call GetMem('TTT','Allo','Real',ipTTT,iSize3)
      Call GetMem('ExtPot*TTT','Allo','Real',ipB,nGrdPt)
      Call InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,
     &           ipExt,ipB,ipIsMM)
*
* Read the old MM electrostatic potential (and derivatives) from PotFile
*
      IPotFl=IsFreeUnit(IPotFl)
      Call Molcas_Open(IPotFl,'ESPF.EXTPOT')
      Line = Get_Ln(IPotFl)
      Call Get_I1(1,nChg)
      If (nChg .ne. 0) Then
         Write(6,*) 'ESPF: nChg ne 0 in h1_espf'
         Call Abend()
      End If
      Do iAt = 1, natom
         Line = Get_Ln(IPotFl)
         Call Get_I1(1,jAt)
         Call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
      End Do
      Close(IPotFl)
*
* Compute the quantum atomic multipoles
*
      Call GetMem('ESPFMltp','Allo','Real',ipMltp,nMult)
      Call espf_mltp(natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,
     &               ipIsMM,ipExt,iPL-1)
*
* Run Tinker
*
      UpdateVMM = .False.
      If (ipOldMltp .ne. ip_Dummy) Then
         sum1 = Zero
         sum2 = Zero
         sum3 = Zero
         sum4 = Zero
         Do iMlt = 1, nMult, MltOrd
            sum1 = Abs(Work(ipMltp+iMlt-1)-Work(ipOldMltp+iMlt-1))
            UpdateVMM = UpdateVMM .or. (sum1 .gt. 0.001d0)
            If (MltOrd .eq. 4) Then
            sum2 = sum2+(Work(ipMltp+iMlt  )-Work(ipOldMltp+iMlt  ))**2
            sum3 = sum3+(Work(ipMltp+iMlt+1)-Work(ipOldMltp+iMlt+1))**2
            sum4 = sum4+(Work(ipMltp+iMlt+2)-Work(ipOldMltp+iMlt+2))**2
            End If
         End Do
         rms2 = sqrt(sum2/nMult)
         rms3 = sqrt(sum3/nMult)
         rms4 = sqrt(sum4/nMult)
         If (MltOrd .eq. 4) Then
            UpdateVMM = UpdateVMM .or. (rms2 .gt. -0.01d0)
            UpdateVMM = UpdateVMM .or. (rms3 .gt. -0.01d0)
            UpdateVMM = UpdateVMM .or. (rms4 .gt. -0.01d0)
         End If
         Call Free_Work(ipOldMltp)
      Else
         UpdateVMM = .True.
      End If
      iQMchg = 1
      If (First .and. Do_DFT) UpdateVMM = .True.
      If (UpdateVMM) Call RunTinker(natom,Work(ipCord),ipMltp,
     &                              iWork(ipIsMM),MltOrd,
     &                    DynExtPot,iQMchg,natMM,StandAlone,DoDirect)
*
* Read the MM electrostatic potential (and derivatives) from PotFile
*
      IPotFl=IsFreeUnit(IPotFl)
      Call Molcas_Open(IPotFl,'ESPF.EXTPOT')
      Line = Get_Ln(IPotFl)
      Call Get_I1(1,nChg)
      If (nChg .ne. 0) Then
         Write(6,*) 'ESPF: nChg ne 0 in h1_espf'
         Call Abend()
      End If
      Do iAt = 1, natom
         Line = Get_Ln(IPotFl)
         Call Get_I1(1,jAt)
         Call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
      End Do
      Close(IPotFl)
*
* Recompute the B matrix
*
      Call InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,ipTTT,
     &           ipExt,ipB,ipIsMM)
*
* Compute the modification of the core hamiltonian
*
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call StatusLine(' espf:',' Computing energy components')
      Call espf_energy(nBas(0),natom,nGrdPt,ipExt,ipGrid,ipB,h1,nh1,
     &                 RepNuc,RealDummy,DoTinker,DoGromacs,DynExtPot)
*
* Save the modified h1 matrix
*
      Call Put_Temp('h1    XX',h1,nh1)
      Call Put_Temp('PotNucXX',[RepNuc],1)
      If (.not.DynExtPot) Then
         Call Put_Temp('h1_raw  ',h1,nh1)
         Call Put_Temp('PotNuc00',[RepNuc],1)
      End If
*
* Save data in the ESPF.DATA file
*
      Call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,
     &                DoGromacs,lMorok,ipMltp,nMult,ipIsMM,natom,
     &                .False.,.False.,DoDirect)
*
* The end
*
      Call GetMem('ExtPot*TTT','Free','Real',ipB,nGrdPt)
      Call GetMem('ExtPot','Free','Real',ipExt,natom*MxExtPotComp)
      Call GetMem('TTT','Free','Real',ipTTT,iSize3)
      Call GetMem('TT','Free','Real',ipTT,iSize2)
      Call GetMem('CartTensor','Free','Real',ipT,iSize1)
      Call GetMem('ESPFMltp','Free','Real',ipMltp,nMult)
      Call GetMem('ESPF_Grid','Free','Real',ipGrid,3*nGrdPt)
      Call GetMem('IsMM for atoms','Free','Inte',ipIsMM,natom)
      Call GetMem('AtomCoord','Free','Real',ipCord,3*natom)
*
      Return
      End
