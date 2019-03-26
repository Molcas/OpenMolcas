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
      Subroutine espf_analysis(lSave)
      Implicit Real*8 (A-H,O-Z)
#include "espf.fh"
*
      Logical Exist,DoTinker,DoGromacs,lMorok,DoDirect,lSave
      Character*10 ESPFKey
      Character*180 ESPFLine, Get_Ln
      External Get_Ln
*
      Call QEnter('espf_analysis')
      iPL = iPL_espf()
*
      If (iPL.ge.2) Then
         Write(6,*)
         Call CollapseOutput(1,'   ESPF analysis')
         Write(6,'(3X,A)')     '   -------------'
      End If
*
* Recover some ESPF data
*
      MltOrd = 0
      iRMax = 0
      DeltaR = Zero
      iGrdTyp = 0
      nGrdPt = 0
      DoTinker = .False.
      DoGromacs = .False.
      lMorok = .False.
      DoDirect = .False.
      Call F_Inquire('ESPF.DATA',Exist)
      If (Exist) Then
         IPotFl = IsFreeUnit(1)
         Call Molcas_Open(IPotFl,'ESPF.DATA')
10       ESPFLine = Get_Ln(IPotFl)
         ESPFKey = ESPFLine(1:10)
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
         Else If (ESPFKey.eq.'GRID      ') Then
            Call Get_I1(2,nGrdPt)
         Else If (ESPFKey.eq.'TINKER    ') Then
            DoTinker = .True.
         Else If (ESPFKey.eq.'GROMACS   ') Then
            DoGromacs = .True.
         Else If (ESPFKey.eq.'LA_MOROK  ') Then
            lMorok = .True.
         Else If (ESPFKey.eq.'DIRECT    ') Then
            DoDirect = .True.
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
      Call espf_init(natom,nAtQM,ipCord,ipIsMM,ipExt)
      nMult = MltOrd * nAtQM
*
* Read the ESPF potential (and derivatives) from PotFile
*
      IPotFl = IsFreeUnit(33)
      IPotFl=IsFreeUnit(IPotFl)
      Call Molcas_Open(IPotFl,'ESPF.EXTPOT')
      ESPFLine = Get_Ln(IPotFl)
      Call Get_I1(1,nChg)
      If (nChg.ne.0) Then
         Write(6,*) 'ESPF: nChg ne 0 in espf_analysis'
         Call Abend()
      End If
      Do iAt = 1, natom
         ESPFLine = Get_Ln(IPotFl)
         Call Get_I1(1,jAt)
         Call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
      End Do
      Close(IPotFl)
*
* Compute the grid around the molecule
*
*      nGrdPt = 0
      ipGrid = ip_Dummy
      ipDGrd = ip_Dummy
      If (iGrdTyp.eq.1)
     &   Call GetMem('ESPF_Grid','Allo','Real',ipGrid,3*nGrdPt)
      Call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,.False.,
     &   ipIsMM,iGrdTyp,ipDGrd,nAtQM)
*
* Compute the cartesian tensor T, TtT^-1, [TtT^-1]Tt
* and B=ExtPot[TtT^-1]Tt
* Tt means the transpose of T
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
* Now the analysis
*
      Call GetMem('ESPFMltp','Allo','Real',ipMltp,nMult)
      Call espf_mltp(natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,
     &               ipIsMM,ipExt,iPL+1)
      Call Add_Info('ESPF multipoles',Work(ipMltp),nMult,6)
*
* Save some data
*
      If (lSave) Call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,
     &                DoTinker,DoGromacs,lMorok,ipMltp,nMult,ipIsMM,
     &                natom,.False.,.False.,DoDirect)
*
* The end
*
      If (iPL.ge.2) Then
         Call CollapseOutput(0,'   ESPF analysis')
         Write(6,*)
      End If
      Call GetMem('ESPFMltp','Free','Real',ipMltp,nMult)
      Call GetMem('ExtPot*TTT','Free','Real',ipB,nGrdPt)
      Call GetMem('TTT','Free','Real',ipTTT,iSize3)
      Call GetMem('TT','Free','Real',ipTT,iSize2)
      Call GetMem('CartTensor','Free','Real',ipT,iSize1)
      Call GetMem('ESPF_Grid','Free','Real',ipGrid,3*nGrdPt)
      Call GetMem('ExtPot','Free','Real',ipExt,natom*MxExtPotComp)
      Call GetMem('IsMM for atoms','Free','Inte',ipIsMM,natom)
      Call GetMem('AtomCoord','Free','Real',ipCord,3*natom)
*
      Call ClsSew()
*
      Call QExit('espf_analysis')
      Return
      End
