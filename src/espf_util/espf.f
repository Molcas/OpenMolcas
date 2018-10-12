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
      Subroutine espf (ireturn,StandAlone)
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
*
* ESPF Module
*
#include "espf.fh"
*
#include "nsd.fh"
#include "disp.fh"
#include "setup.fh"
#include "status.fh"
#include "lundio.fh"
#include "print.fh"
#include "nac.fh"
      Character Label*8
      Logical Forces,Show_espf,StandAlone,DoTinker,DoGromacs,DynExtPot
      Logical lMorok,DoDirect,isNAC_tmp
*
      Call QEnter('espf')
      iReturn=99
*
*-----Print
*
      iPL = iPL_espf()
*
* Some warnings
*
      Call Get_iScalar('nSym',nSym)
      nIrrep=nSym
      If (nSym.gt.1) Then
         Write(6,'(A)')' Symmetry cannot be used together with ESPF.'
         Call Quit_OnUserError()
       End If
*
* Set on the System Bit 11
*
      Call Get_iScalar('System Bitswitch',iOption)
      iOption=iOr(iOption,2**11)
      Call Put_iScalar('System Bitswitch',iOption)
*
* Some initializations
*
      Call espf_init(natom,nAtQM,ipCord,ipIsMM,ipExt)
      ipMltp = ip_Dummy
      nGrdPt = 0
      ipGrid = ip_Dummy
      ipDGrd = ip_Dummy
      isNAC_tmp = isNAC
*
* Read the input and compute the external potential
*
      Call StatusLine(' espf:',' Reading input')
      Call ReadIn_ESPF(natom,ipCord,ipExt,MltOrd,iRMax,DeltaR,Forces,
     &                Show_espf,ipIsMM,StandAlone,iGrdTyp,DoTinker,
     &                DoGromacs,DynExtPot,ipMltp,natMM,lMorok,DoDirect,
     &                ipGradCl,EnergyCl)
*
* If the present calculation does not use ESPF but the Direct scheme
*
      If (DoDirect) Then
         Call No_ESPF(natom,Forces,DoTinker)
         Goto 98
      End If
*
      nMult = MltOrd * nAtQM
      If (iPL.ge.2) Write(6,'(/,A,I2,A,i4,A,i6)')
     &   ' Number of ESPF operators (nMult=',MltOrd,' * nAtQM=',
     &     nAtQM,'): ',nMult
*
* Compute the grid around the molecule
*
      Call StatusLine(' espf:',' Making the grid')
      If (iGrdTyp.eq.1) Then
         If(nGrdPt.eq.0) Call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,
     &      DeltaR,Forces,ipIsMM,-iGrdTyp,ipDGrd,nAtQM)
         Call GetMem('ESPF_Grid','ALLO','REAL',ipGrid,3*nGrdPt)
         Call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,Forces,
     &     ipIsMM,iGrdTyp,ipDGrd,nAtQM)
         If (iPL.ge.2) Then
            Write(6,'(A)') ' PNT Grid (Warning: no grid derivatives)'
            Write(6,'(A)') ' (C. Chipot and J. Angyan,'//
     &                     ' Henri Poincare University, Nancy, France)'
            Write(6,'(5X,I5,A)') nGrdPt,' grid points'
         End If
      Else
         Call MkGrid(natom,ipCord,ipGrid,nGrdPt,iRMax,DeltaR,Forces,
     &      ipIsMM,iGrdTyp,ipDGrd,nAtQM)
         If (iPL.ge.2) Then
            Write(6,'(A)')' GEPOL Grid, using United Atoms radii'
            Write(6,'(5X,I5,A)') nGrdPt,' grid points'
         End If
      End If
*
*     If this is a standalone call to &ESPF, there are 2 options:
*        1) static external potential: compute here the ESPF contributions
*        2) dynamic external potential: nothing more to compute here
*
      If (StandAlone .and. DynExtPot) Goto 98
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
      Call GetMem('DerivB','Allo','Real',ipDB,nGrdPt*3*nAtQM)
      Call InitDB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,
     &           ipTTT,ipExt,ipDB,ipIsMM,iRMax,DeltaR,iGrdTyp,ipDGrd)
      If (iGrdTyp.eq.2 .and. ipDGrd.ne.ip_Dummy)
     &   Call GetMem('ESPF_DGrid','Free','Real',ipDGrd,3*nGrdPt*3*nAtQM)
*
* Here we must distinguish between an energy run and a gradient run
*
      If (.not.Forces) Then
         Call StatusLine(' espf:',' Computing energy components')
         Call Get_iArray('nBas',nBas,nSym)
         nBas0 = nBas(0)
         nSize = nBas0*(nBas0+1)/2 + 4
         Call Allocate_Work(ipH,nSize)
         iComp = 1
         iSyLbl = 1
         Label = 'OneHam  '
         iRc = -1
         Call iRdOne(iRc,1,Label,iComp,nInts,iSyLbl)
         If (iRc.ne.0) Then
            Write (6,'(A)')' ESPF: Error reading ONEINT'
            Write (6,'(A,A8)')' Label = ',Label
            Call QTrace()
            Call Abend()
         End If
         If (nInts+4.ne.nSize) Then
            Write (6,'(A,2I5)')' ESPF: nInts+4.ne.nSize',nInts+4,nSize
            Call QTrace
            Call Abend()
         End If
         iRc = -1
         Call RdOne(iRc,0,Label,iComp,Work(ipH),iSyLbl)
         Call Get_dScalar('PotNuc',RepNuc)
         Call espf_energy(nBas0,natom,nGrdPt,ipExt,ipGrid,ipB,
     &                    Work(ipH),nSize-4,RepNuc,EnergyCl,DoTinker,
     &                    DoGromacs,DynExtPot)
         Call Put_dScalar('PotNuc',RepNuc)
         Call WrOne(iRc,0,Label,iComp,Work(ipH),iSyLbl)
         If (iRC.ne.0) then
            Write (6,*)'ESPF: Error writing to ONEINT'
            Write (6,'(A,A8)')'Label=',Label
            Call QTrace
            Call Abend()
         End If
         Call Free_Work(ipH)
         If (iPL.ge.3) Write(6,*) 'The 1-e hamiltonian is now updated.'
         If (iPL.ge.2) Write(6,'(A,F16.10)')
     &                 ' Nuclear energy, including Ext Pot = ',RepNuc
      Else
         Call StatusLine(' espf:',' Computing gradient components')
         Call espf_grad(natom,nAtQM,nGrdPt,ipExt,ipGrid,ipB,ipDB,ipIsMM,
     &                  ipGradCl,DoTinker,DoGromacs)
         If (ipMltp.eq.ip_Dummy)
     &      Call GetMem('ESPFMltp','Allo','Real',ipMltp,nMult)
         Call espf_mltp(natom,MltOrd,nMult,nGrdPt,ipTTT,ipMltp,ipGrid,
     &                  ipIsMM,ipExt,iPL)
      End If
      Call ClsSew(1)
*
98    Continue
*
* Save data in the ESPF.DATA file
*
      Call espf_write(MltOrd,iRMax,DeltaR,iGrdTyp,nGrdPt,DoTinker,
     &                DoGromacs,lMorok,ipMltp,nMult,ipIsMM,natom,
     &                Show_espf,Forces,DoDirect)
*
* Exit
*
      isNAC = isNAC_tmp
      If (.not. (StandAlone .and. DynExtPot)) Then
         Call GetMem('DerivB','FREE','REAL',ipDB,nGrdPt*3*nAtQM)
         Call GetMem('CartTensor','FREE','REAL',ipT,iSize1)
         Call GetMem('TT','FREE','REAL',ipTT,iSize2)
         Call GetMem('TTT','FREE','REAL',ipTTT,iSize3)
         Call GetMem('ExtPot*TTT','FREE','REAL',ipB,nGrdPt)
      End If
      If (ipMltp.ne.ip_Dummy)
     &   Call GetMem('ESPFMltp','FREE','REAL',ipMltp,nMult)
      If (ipGrid.ne.ip_Dummy)
     &   Call GetMem('ESPF_Grid','FREE','REAL',ipGrid,3*nGrdPt)
      Call GetMem('AtomCoord','FREE','REAL',ipCord,3*natom)
      Call GetMem('ExtPot','FREE','REAL',ipExt,natom*MxExtPotComp)
      Call GetMem('IsMM for atoms','Free','Inte',ipIsMM,natom)
      If (DoGromacs.and.Forces) Then
         Call GetMem('GradCl','FREE','REAL',ipGradCl,3*natom)
      End If
*
*     Slapaf needs to know that the gradient is NOT translational
*     and rotational invariant.
*
      If (.not.Forces .and. natMM.gt.0) Then
         Call Get_iScalar('System BitSwitch',iOption)
         iOption=iOr(iOption,2**7)
         iOption=iOr(iOption,2**8)
         Call Put_iScalar('System BitSwitch',iOption)
      End If
*
      Call QExit('espf')
      iReturn=0
      Return
      End
