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
      Subroutine AverMEP(Kword,Eint,Poli,ici,SumElcPot
     &                  ,NCountField,PertElcInt
     &                  ,iQ_Atoms,nBas,nOcc,natyp,nntyp)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "qminp.fh"
#include "qm1.fh"
#include "qm2.fh"
#include "qmcom.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"
#include "warnings.fh"
      Dimension Eint(MxQCen,10),Poli(MxQCen,10)
      Dimension SumElcPot(MxQCen,10)
      Dimension PertElcInt(MxBas*(MxBas+1)/2),SumOld(MxQCen,10)
      Dimension iCent(MxBas*MxBas)
      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6)
      Dimension nOcc(*),natyp(*),ForceNuc(MxAt,3)

      Character*4 Kword
      Character*20 MemLab,MemLab1
      Logical Exist
      Dimension iiDum(1)

      Call UpCase(Kword)
***************
* This subroutine include three different options. All have to do with the
* calculation of a Mean Electrostatic Potential, Field and Field gradients,
* and to evalute the perturbation of them in the One electron Hamiltoniam
* this perturbation is added to the SEWARD One Electron File.
* Option 1: Add the components of Potential, etc
* Option 2: Obtain the Average
* Option 3: Calculate the electrostatic perturbation energy integrals
*           and add them to the one-electron file.
* Calculation involve up to the field gradients because the charge density
* is expanded to the quadrupoles. If the expansion is bigger the number 10
* must be changed, but also all the eqscf and eqras subroutines shoud be
* change in the last option the array nMlt is used instead of the
* number so if a smaller expantion is used, non problem,
* since this index take care of that.
******************
*
*-- The keywords and their labels.
*
      If(Kword(1:4).eq.'ADD ') Go To 101
      If(Kword(1:4).eq.'AVER') Go To 102
      If(Kword(1:4).eq.'PERT') Go To 103

101   Continue
      Do 1001, i=1,iCi
        Do 1002, j=1,10 !Charges (1),Dipoles(3),Quadrupoles(6)
          SumOld(i,j)=SumElcPot(i,j)
          SumElcPot(i,j)=SumOld(i,j)+Eint(i,j)+Poli(i,j)
1002    Continue
1001  Continue
      If(iPrint.ge.9) then
        Write(6,*)'Total Sum Potential'
        Do 1010, i=1,iCi
           Write(6,*)(SumElcPot(i,j),j=1,10)
1010    Continue
      Endif
      Go to 9999

102   Continue
      Do 2001, i=1,iCi
        Do 2002, j=1,10 !Charges (1),Dipoles(3),Quadrupoles(6)
          AvElcPot(i,j)=SumElcPot(i,j)/Dble(NCountField)
2002    Continue
*
*-- The order of Field gradients is changed in order to follow
*-- the same order than Molcas
*
      AvTemp=AvElcPot(i,8)         ! This change is due to the
      AvElcPot(i,8)=AvElcPot(i,7)  ! different order of quadrupoles
      AvElcPot(i,7)=AvTemp         ! in QmStat and Molcas.
2001  Continue                     ! QmStat:xx,xy,yy,xz,yz,zz
                                   ! Molcas:xx,xy,xz,yy,yz,zz

********************************
* This multiplication comes because the off-diagonal
* quadrupoles must be multiply by two since we use
* a triangular form to compute the Interaction
* Energy with the Electric Field Gradient.
* Since it is easier multiply the Average potential
* than the quadrupole for each pair of basis, we perform
* the multiplication here
***********************
      AvElcPot(i,6)=2.0d0*AvElcPot(i,6)
      AvElcPot(i,7)=2.0d0*AvElcPot(i,7)
      AvElcPot(i,9)=2.0d0*AvElcPot(i,9)
***********************
      If(iPrint.ge.9) then
        Write(6,*)'Total Averg Potential'
        Do 1020, i=1,iCi
           Write(6,*)(AvElcPot(i,j),j=1,10)
1020    Continue
      Endif

      Go to 9999


103   Continue
*
*----First we read the multipoles expansion for each pair of basis.
*----The index iCent(i) will give us to which center belongs each pair of basis.
*
      Call GetMem('Dummy','Allo','Inte',iDum,nBas**2)
      Call MultiNew(iQ_Atoms,nBas,nOcc,natyp,nntyp,iMME
     &             ,iCent,iWork(iDum),nMlt,outxyz,SlExpQ,.false.)
      Call GetMem('Dummy','Free','Inte',iDum,nBas**2)


**********************
* Calculate the forces for the nuclei
* these forces will compensate parcially
* the forces due to the electrons
* They will be printed and added to the
* RUNFILE in the optimization procedure
* after Alaska module
*********************
* This model do not work
* To calculate the forces in the nuclei
* with a Slater representation since
* you have to calculate the field
* in a set of point charges and not in
* distributed charges as the field is calculated
* when used Slater representation
* also there are a more dark and complicated
* problem about the string interaction keeping
* together the distributed electronic charge
* and the point nuclear charge under different
* forces.
**********************
      Do 2030, i=1,iQ_Atoms
         Do 2032, j=1,3
           ForceNuc(i,j)=ChaNuc(i)*AvElcPot(i,j+1)
2032     Continue
2030  Continue
      iLuField=63
      iLuField=IsFreeUnit(iLuField)
      Call OpnFl(FieldNuc,iLuField,Exist)
         Write(6,*)'FieldNuc',FieldNuc
      Do 2036, i=1,iQ_Atoms
         Write(iLuField,*)(ForceNuc(i,j),j=1,3)
2036  Continue
      Close(iLuField)

      If(iPrint.ge.9) then
        Write(6,*)'Nuclei charge and Forces'
        Do 1030, i=1,iQ_Atoms
           Write(6,*)ChaNuc(i),(ForceNuc(i,j),j=1,3)
1030    Continue
      Endif
*********************

      nTyp=0
      Do 3000, i=1,nMlt
        nTyp=nTyp+i*(i+1)/2
3000  Continue
      Do 3001, i=1,(nBas*(nBas+1)/2)
          PertElcInt(i)=0.0d0
3001  Continue

*-- Put quadrupoles in Buckinghamform.
*
      Do 191, i1=1,nBas
        Do 192, i2=1,i1
          indMME=i2+i1*(i1-1)/2
           Do 194, j=5,10
             Work(iMME(j)+indMME-1)=
     &                Work(iMME(j)+indMME-1)*1.5
194        Continue
          Tra=Work(iMME(5)+indMME-1)
     &    +Work(iMME(8)+indMME-1)+Work(iMME(10)+indMME-1)
          Tra=Tra/3
          Work(iMME(5)+indMME-1)=Work(iMME(5)+indMME-1)-Tra
          Work(iMME(8)+indMME-1)=Work(iMME(8)+indMME-1)-Tra
          Work(iMME(10)+indMME-1)=Work(iMME(10)+indMME-1)-Tra
192     Continue
191   Continue



      irc=-1
      Lu_One=49
      Lu_One=IsFreeUnit(Lu_One)
      Call OpnOne(irc,0,'ONEINT',Lu_One)
      If(irc.ne.0) then
        Write(6,*)
        Write(6,*)'ERROR! Could not open one-electron integral file.'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif

*
*---We Read the size of the unperturbed Hamiltonian 'OneHam 0' in OneInt.
*
      irc=-1
      iOpt=1
      iSmLbl=1
      nSize=0
      Call iRdOne(irc,iOpt,'OneHam 0',1,iiDum,iSmLbl)
      nSize=iiDum(1)
      If(irc.ne.0) then
        Write(6,*)
        Write(6,*)'ERROR! Failed to read number of one-electron i'
     &//'ntegrals.'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif
      If(nSize.eq.0) then
        Write(6,*)
        Write(6,*)'ERROR! Problem reading size of unperturbed'
     &//' Hamiltonian in OneInt'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif

*---Memory allocation for the unperturbed Hamiltonian
      Write(MemLab,*)'MAver'
      Call GetMem(MemLab,'Allo','Real',iH0,nSize+4)
      irc=-1
      iOpt=0
      iSmLbl=0

*---Read the unperturbed Hamiltonian
      Call RdOne(irc,iOpt,'OneHam 0',1
     &          ,Work(iH0),iSmLbl) !Collect non perturbed integrals
      Write(MemLab1,*)'MAver1'
      Call GetMem(MemLab1,'Allo','Real',iH1,nSize+4)
      If(iPrint.ge.9) then
        Call TriPrt('Non Perturb One-e',' ',Work(iH0),nBas)
      Endif
*

*---We perform the multiplication for each pair of basis in a triangular form.
*---The perturbation is added to the unperturbed Hamiltonian 'iH0'.
*
      kaunta=0
      Do 3003, iB1=1,nBas
        Do 3004, iB2=1,iB1
          kaunta=kaunta+1
          indMME=iB2+iB1*(iB1-1)/2
          Do 3005, iTyp=1,nTyp
            PertElcInt(indMME)=PertElcInt(indMME)
     &        +AvElcPot(iCent(kaunta),iTyp)*Work(iMME(iTyp)+indMME-1)
3005      Continue
          Work(iH1+kaunta-1)=Work(iH0+kaunta-1)+PertElcInt(indMME)
3004    Continue
3003  Continue

      If(iPrint.ge.9) then
        Call TriPrt('H0+Elec One-e',' ',Work(iH1),nBas)
      Endif

*---The non-Electrostatic perturbation is added. The PertNElcInt array comes
*---throught the include file qminp.fh.

      If(iPrint.ge.10) then
        Call TriPrt('PertNElcInt-e',' ',PertNElcInt,nBas)
      Endif

      iTriBasQ=nBas*(nBas+1)/2
      Call DaxPy_(iTriBasQ,ONE,PertNElcInt,iONE,Work(iH1),iONE)

      If(iPrint.ge.9) then
        Call TriPrt('H0+Elec+nonEl One-e',' ',Work(iH1),nBas)
      Endif

*----The perturbed Hamiltonian 'H1' is writen in OneInt.
      irc=-1
      iOpt=0
      iSmLbl=1
      Call WrOne(irc,iOpt,'OneHam  ',1
     &          ,Work(iH1),iSmLbl) !Write  perturbed integrals
      If(iPrint.ge.9) then
        Call TriPrt('Perturb One-e',' ',Work(iH1),nBas)
      Endif

      If(iPrint.ge.10) then
        Call TriPrt('Non Perturb One-e AGAIN',' ',Work(iH0),nBas)
      Endif

      Call ClsOne(irc,Lu_One)

      Call GetMem(MemLab,'Free','Real',iH0,nSize+4)
      Call GetMem(MemLab1,'Free','Real',iH1,nSize+4)


9999  Continue

      Return
      End
