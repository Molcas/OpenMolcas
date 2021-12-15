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
      Subroutine espf_energy (nBas0,natom,nGrdPt,ipExt,ipGrid,ipB,h1,
     &                        nh1,RepNuc,EnergyCl,DoTinker,DoGromacs,
     &                        DynExtPot)
      Implicit Real*8 (A-H,O-Z)
*
*      Compute the integrals <mu|B/R_grid|nu>, where B weights every
*      point of the grid and R_grid is the distance to one grid point.
*
#include "espf.fh"
*
      Character*180 Line,Get_Ln
      External Get_Ln
      Character*8 Label
      Logical DoTinker,DoGromacs,DynExtPot
      Real*8 h1(nh1)
      Dimension opnuc(1),idum(1)
*
      iPL = iPL_espf()
*
* Read the MM contribution to the total energy and add it
* to the Nuclear Repulsion term
*
      If (DoTinker) Then
         ITkQMMM = IsFreeUnit(30)
         Call Molcas_Open (ITkQMMM,'QMMM')
         Line = ' '
         Do While (Index(Line,'TheEnd ') .eq. 0)
            Line=Get_Ln(ITkQMMM)
            If (Index(Line,'MMEnergy ').ne.0) Call Get_F1(2,TkE)
         End Do
         Close (ITkQMMM)
         TkE = TkE * ToHartree
         RepNuc_old = RepNuc
         RepNuc = RepNuc + TkE
         If (iPL.ge.3) Write(6,3000) RepNuc_old,TkE,RepNuc
      Else If (DoGromacs) Then
         RepNuc_old = RepNuc
         RepNuc = RepNuc + EnergyCl
         If (iPL.ge.3) Write(6,3000) RepNuc_old,EnergyCl,RepNuc
      End If
*
 3000 Format(/,' RepNuc + MM = ',F13.8,' + ',F13.8,' = ',F13.8)
*
*     Call to DrvPot to compute the integrals
*     Here we don't care about opnuc (nuclear potential)
*
      nSize = nBas0*(nBas0+1)/2 + 4
      If (nSize .ne. (nh1+4)) Then
         Write(6,*) 'In espf_energy, nSize ne nh1',nSize,nh1+4
         Call Abend()
      End If
      opnuc = Dum
*
      ncmp = 1
      iAddPot = 1
      If (iPL.ge.4) Then
         Do i = 1, NGrdPt
            Write(6,1234) i,(Work(ipGrid+(i-1)*3+j),j=0,2),
     &                       Work(ipB+(i-1))
 1234       Format('Grid point ',I4,/,4F12.6)
         End Do
      End If
      Call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipB),nGrdPt,iAddPot)
      Label = 'Pot     '
      iComp = 1
      iSyLbl = 1
      iRc = -1
      Call iRdOne(iRc,1,Label,iComp,idum,iSyLbl)
      nInts=idum(1)
      If (iRc.ne.0) Then
         Write (6,'(A)')' ESPF: Error reading ONEINT'
         Write (6,'(A,A8)')' Label = ',Label
         Call Abend()
      End If
      If (nInts+4.ne.nSize) Then
         Write (6,'(A,2I5)')' ESPF: nInts+4.ne.nSize',nInts+4,nSize
         Call Abend()
      End If
      Call GetMem('IntOnGrid','Allo','Real',ipInt,nSize)
      Call RdOne(iRc,0,Label,iComp,Work(ipInt),iSyLbl)
      If(iPL .ge. 4) Call TriPrt(Label,' ',Work(ipInt),nBas0)
*
*     The core Hamiltonian must be updated
*
      call daxpy_(nInts,One,Work(ipInt),1,h1,1)
      If (DynExtPot) Then
         iSyLbl=1
         iRc=-1
         iOpt=0
         iComp=1
         Label='OneHamRF'
         Call WrOne(iRc,iOpt,Label,iComp,Work(ipInt),iSyLbl)
      End If
      Call GetMem('IntOnGrid','Free','Real',ipInt,nSize)
*
*     The electrostatic energy between the external potential
*     and the nuclei is added to the nuclear energy
*
      EQC = ExtNuc(ipExt,natom)
      RepNuc = RepNuc + EQC
      If (IsStructure().eq.1) Then
        Call Add_Info('PotNuc',[RepNuc],1,6)
      Else
        Call Add_Info('PotNuc',[RepNuc],1,12)
      End If
*
      Return
      End
