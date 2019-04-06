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
      Subroutine Print_MPPROP(rMP,xrMP,xnrMP,nij,nElem,lMax,EC,Polar,
     &                        Lbl,nAtoms,iANr,NoField,CoC,Coor,
     &                        nOcOb,Energy_Without_FFPT,ip_Ene_Occ,
     &                        MpProp_Level,Bond_Threshold,nReal_Centers)
*
*  rMP : Multipole moments moved to center of charge
* xrMP : Multipole moments kept on the expansion center
*
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
      Real*8 rMP(nij,nElem),xrMP(nij,nElem),xnrMP(nij,nElem),CoC(3),
     &       EC(3,nij),Polar(6,nij),Polar_M(6),Coor(3,nAtoms)
      Integer iANr(nAtoms),nBas(8), lHeader
      Parameter (lHeader = 144)
      Character*(lHeader) Header
      Character*(LENIN) Lbl(nAtoms)
      Character*8 Method,Label
      Character*6 fName
      Logical NoField, Exist, Text, Bond_OK, Check_Bond
      Logical Use_Two_Centers_For_Atoms, Found
      Parameter (Use_Two_Centers_For_Atoms = .False.)
*
      MxMP = lMax
      If (lMax .gt. MpProp_Level) MxMP = MpProp_Level
*
      Call Qpg_carray('Relax Method',Found,iSize)
      If(Found) Then
         Call Get_cArray('Relax Method',Method,8)
      Else
         Method = 'UNDEF   '
      EndIf
      Call Get_cArray('Seward Title',Header,lHeader)
      j = 0
      Text = .False.
      Last_NonBlank = 0
      Do i = 1, 72
         If (Text .or. Header(i:i) .ne. ' ') Then
            Text = .True.
            j = j + 1
            Header(j:j) = Header(i:i)
            If (Header(i:i) .ne. ' ') Last_NonBlank = j
         End If
      End Do
*
      fName = 'MPPROP'
      iUnit = 11
      Call OpnFl(fName,iUnit,Exist)

      nCenters = nReal_Centers
      If (Use_Two_Centers_For_Atoms) Then
         nCenters = nCenters + nAtoms
      End If

      Write(iUnit,'(A)')
     &   '**************************************************'
      Write(iUnit,'(A)') '* Molecule'
      Write(iUnit,'(A)') Header(1:Last_NonBlank)
      Write(iUnit,'(A)') '* Method'
      Write(iUnit,'(A)') Method
      Write(iUnit,'(A)') '* Level of Multipoles and Polarizabilities'
      Write(iUnit,'(2I5)') MxMP,1
      Write(iUnit,'(A)') '* All centers'
      Write(iUnit,'(I5)') nCenters

* Insert atoms
      Do iAtom = 1, nAtoms
         ii = iAtom*(iAtom+1)/2
         Call ReExpand(xnrMP,nij,nElem,CoC,EC(1,ii),ii,lMax)

* Informations on the atom
         Write(iUnit, '(2I5,4X,A)') 2,1,Lbl(iAtom)
* The expansion center
         Write(iUnit,'(3F20.10)') (EC(iElem,ii),iElem = 1,3)
* The multipole moments
         iStrt = 1
         Do l = 0, MxMP
            iEnd  = iStrt + (l+1)*(l+2)/2 - 1
            Write(iUnit,'(3F20.10)') (xrMP(ii,iElem)+xnrMP(ii,iElem),
     &                                            iElem = iStrt,iEnd)
            iStrt = iEnd + 1
         End Do
* The size parameters (never higher than quadrupole)
         lMax_for_size = min(MxMP,2)
         iStrt = 1
         Do l = 0, lMax_for_size
            iEnd  = iStrt + (l+1)*(l+2)/2 - 1
            Write(iUnit,'(3F20.10)') (Zero,iElem = iStrt,iEnd)
            iStrt = iEnd + 1
         End Do
* Polarizabilities
         If (NoField) Then
            Write(iUnit,'(3F20.10)') (Zero,iElem=1,6)
         Else
            Polar_M(1) = Polar(1,ii)
            Polar_M(2) = Polar(2,ii)
            Polar_M(3) = Polar(4,ii)
            Polar_M(4) = Polar(3,ii)
            Polar_M(5) = Polar(5,ii)
            Polar_M(6) = Polar(6,ii)
            Write(iUnit,'(3F20.10)') (Polar_M(iElem),iElem=1,6)
         End If

         Call ReExpand(xnrMP,nij,nElem,EC(1,ii),CoC,ii,lMax)
      End Do
      If (Use_Two_Centers_For_Atoms) Then
* Insert real atoms and fake values
         Do iAtom = 1, nAtoms

* Informations on the center
            Write(iUnit, '(2I5,4X,A)') iANr(iAtom),1,Lbl(iAtom)
* The expansion center
            Write(iUnit,'(3F20.10)') (Coor(iElem,iAtom),iElem = 1,3)
* The multipole moments
            iStrt = 1
            Do l = 0, MxMP
               iEnd  = iStrt + (l+1)*(l+2)/2 - 1
               Write(iUnit,'(3F20.10)') (Zero,iElem = iStrt,iEnd)
               iStrt = iEnd + 1
            End Do
* The size parameters (never higher than quadrupole)
            lMax_for_size = min(MxMP,2)
            iStrt = 1
            Do l = 0, lMax_for_size
               iEnd  = iStrt + (l+1)*(l+2)/2 - 1
               Write(iUnit,'(3F20.10)') (One,iElem = iStrt,iEnd)
               iStrt = iEnd + 1
            End Do
* Polarizabilities
            Write(iUnit,'(3F20.10)') (Zero,iElem=1,6)

         End Do
      End If
* Insert bonds
      Do jAtom = 1, nAtoms
         Do iAtom = jAtom+1, nAtoms
            ij = iAtom*(iAtom-1)/2 + jAtom

            Bond_OK = Check_Bond(Coor(1,iAtom),Coor(1,jAtom),
     &                           iANr(iAtom),iANr(jAtom),Bond_Threshold)
            If (Bond_OK) Then
               Label = '       '
               Write(Label,'(I3,A,I3)') jAtom,'-',iAtom
               j = 0
               Do i = 1,8
                  If (Label(i:i) .ne. ' ') Then
                     j = j + 1
                     Label(j:j) = Label(i:i)
                  End If
               End Do
* Informations on the center
               Write(iUnit, '(2I5,4X,A)') 2,1,Label(1:j)
* The expansion center
               Write(iUnit,'(3F20.10)') (EC(iElem,ij),iElem = 1,3)
* The multipole moments
               iStrt = 1
               Do l = 0, MxMP
                  iEnd  = iStrt + (l+1)*(l+2)/2 - 1
                  Write(iUnit,'(3F20.10)')
     &               (xrMP(ij,iElem),iElem = iStrt,iEnd)
                  iStrt = iEnd + 1
               End Do
* The size parameters (never higher than quadrupole)
               lMax_for_size = min(MxMP,2)
               iStrt = 1
               Do l = 0, lMax_for_size
                  iEnd  = iStrt + (l+1)*(l+2)/2 - 1
                  Write(iUnit,'(3F20.10)') (One,iElem = iStrt,iEnd)
                  iStrt = iEnd + 1
               End Do
* Polarizabilities
               If (NoField) Then
                  Write(iUnit,'(3F20.10)') (Zero,iElem=1,6)
               Else
                  Polar_M(1) = Polar(1,ij)
                  Polar_M(2) = Polar(2,ij)
                  Polar_M(3) = Polar(4,ij)
                  Polar_M(4) = Polar(3,ij)
                  Polar_M(5) = Polar(5,ij)
                  Polar_M(6) = Polar(6,ij)
                  Write(iUnit,'(3F20.10)') (Polar_M(iElem),iElem=1,6)
               End If
            End If

         End Do
      End Do
*
******** Molecular properties
*
      Write(iUnit,'(A)') '* Molecule properties'
      Call Allocate_Work(iScratch_ele,nElem)
      Call Allocate_Work(iScratch_nuc,nElem)
      Do iElem = 1, nElem
         Work(iScratch_ele+iElem-1)=DDot_(nij,[One],0,rMP(1,iElem),1)
         Work(iScratch_nuc+iElem-1)=DDot_(nij,[One],0,xnrMP(1,iElem),1)
      End Do
      iStrt = 1
      Do l = 0, MxMP
         iEnd  = iStrt + (l+1)*(l+2)/2 - 1
         Write(iUnit,'(3F20.10)')
     &      (Work(iScratch_ele+iElem-1)+Work(iScratch_nuc+iElem-1),
     &                                          iElem = iStrt,iEnd)
         iStrt = iEnd + 1
      End Do
      Call Free_Work(iScratch_ele)
      Call Free_Work(iScratch_nuc)

      If (NoField) Then
         Write(iUnit,'(3F20.10)') (Zero,iElem=1,6)
      Else
         Polar_M(1)=DDot_(nij,[One],0,Polar(1,1),6)
         Polar_M(2)=DDot_(nij,[One],0,Polar(2,1),6)
         Polar_M(3)=DDot_(nij,[One],0,Polar(4,1),6)
         Polar_M(4)=DDot_(nij,[One],0,Polar(3,1),6)
         Polar_M(5)=DDot_(nij,[One],0,Polar(5,1),6)
         Polar_M(6)=DDot_(nij,[One],0,Polar(6,1),6)
         Write(iUnit,'(3F20.10)') (Polar_M(iElem),iElem=1,6)
      End If
*
******** Orbital information
*
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
*
      Write(iUnit,'(A)') '* Orbital information'
      Write(iUnit,'(2I5)') nBas(1), nOcOb
      Write(iUnit,'(F20.10)') Energy_Without_FFPT
      Write(iUnit,'(3F20.10)') (Work(ip_Ene_Occ+i),i=0,nOcOb-1)
*
      Call Free_Work(ip_Ene_Occ)
*
      Close(iUnit)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
