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
      Subroutine isoloop(Double)
      Implicit Real*8(a-h,o-z)
#include "Molcas.fh"
#include "WrkSpc.fh"
      Character*2 Element(MxAtom)
      Logical Double
*                                                                      *
************************************************************************
*                                                                      *
      Write (6,*)
      Call CollapseOutput(1,'   Isotopic shifts:')
      Write (6,'(3X,A)')    '   ----------------'
      Write (6,*)
*
      Call Get_nAtoms_All(nAtoms_All)
      Call Allocate_Work(ipCoor,3*nAtoms_All)
      Call Get_Coord_All(Work(ipCoor),nAtoms_All)
      Call Get_Name_All(Element)
*
      Write (6,*)
      Write (6,*)
      Write (6,*)'****************************************'
      Write (6,*)'*                                      *'
      Write (6,*)'* Isotope shifted frequencies in cm-1  *'
      Write (6,*)'*                                      *'
      Write (6,*)'****************************************'
      Write (6,*)
      n=3*nAtoms_All
      nTemp=6*2*n**2
      Call Getmem('ISOLOOP','ALLO','REAL',ipT,nTemp)
      Call Isotop_i(n,Element,nAtoms_All,Work(ipT),nTemp,Work(ipCoor),
     &              Double)
      Call Getmem('ISOLOOP','FREE','REAL',ipT,nTemp)
      Call Free_Work(ipCoor)
*
      Call CollapseOutput(0,'   Isotopic shifts:')
      Write(6,*)
*
      Return
      End
      Subroutine Isotop_i(n,Element,mAtoms,Temp,nTemp,Coord,Double)
      Use Isotopes
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "periodic_table.fh"
#include "constants2.fh"
#include "WrkSpc.fh"
      Integer AtNum, IsoNum, AtNum2, Subs, Subs2
      Real*8 Temp(nTemp), Coord(3,mAtoms)
      Real*8 MassIso
      Real*8 Mass(MxAtom)
      Character*2 Element(mAtoms)
      Logical EQ, Double, lmass, Found, Changed
      Character*3 cmass(MxAtom)
      Real*8 umass(MxAtom)
*                                                                      *
************************************************************************
*                                                                      *
      ipH=1
*
      Call qpg_iscalar('iMass',lmass)
      If (lmass) Then
         Call Get_iScalar('iMass',iMass)
         Call Get_cArray('cmass',cmass,3*iMass)
         Call Get_dArray('umass',umass,iMass)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Read the unsymmetrized Hessian from run file
*
      Call Get_dArray('FC-Matrix',Temp(ipH),n**2)
*
*---- Put in the initial masses and take backup copy
*
      Call Get_Mass_All(Mass,mAtoms)
      Call GetMem('DefaultMass','ALLO','REAL',ipDefMass,mAtoms)
      Call dCopy_(mAtoms,Mass,1,Work(ipDefMass),1)
*
      Do i=1,mAtoms
         Coord(1,i)=Abs(Coord(1,i))
         Coord(2,i)=Abs(Coord(2,i))
         Coord(3,i)=Abs(Coord(3,i))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Loop over common isotops
*
      ipH2=ipH+2*n**2
      ip1=ipH2+2*n**2
      ip2=ip1+2*n**2
      ipVal=ip2+2*n**2
      ipVec=ipVal+2*n**2
      Call Initialize_Isotopes()
*                                                                      *
************************************************************************
*                                                                      *
*---- Single substitutions
*
      Write (6,*)
      Write (6,*) ' Single substitutions:'
      Write (6,*) ' -----------------------'
      Write (6,*)
      Do i = 1, mAtoms
*
         If (i.gt.1) Then
            If(EQ(Coord(1,i),Coord(1,i-1))) Go To 94
         End If
*
         AtNum=iNuclearChargeFromSymbol(Element(i))
         If (lmass) Then
*
*           Do according to user list of isotopes.
*
            Subs=0
            Do nAt=1,iMass
               If (Element(i).eq.cmass(nAt)) Then
                  Mass(i)=UtoAU*umass(nAt)
                  Subs=1
               End If
            End Do
         Else
*
*           Get list of natural isotopes
*
            Subs=ElementList(AtNum)%Nat
         End If
*
         dMass = Mass(i)
         Do k = 1, Subs
            If (.not.lmass) Then
               IsoNum=ElementList(AtNum)%Isotopes(k)%A
               If (IsoNum.eq.nInt(dMass/UtoAU)) Cycle
               Call Isotope(IsoNum,AtNum,Mass(i))
            End If
*
            call dcopy_(n**2,Temp(ipH),1,Temp(ipH2),1)
            Write(6,*) 'Masses:'
            Write(6,*) '======='
            Write(6,'(20I4)') (nInt(mass(l)/UtoAU),l=1,mAtoms)
            Write(6,*)
            Write(6,*)
            Write(6,*) 'Frequencies:'
            Write(6,*) '============'
            call freq_i(n,Temp(ipH2),mass,Temp(ip1),Temp(ip2),
     &               Temp(ipVec),Temp(ipVal),ineg)
            call GFPrnt_i(Temp(ipVal),n)
         End Do
*        Put back the original mass
         Mass(i) = dMass
 94      Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Full substitutions
*
      Write (6,*)
      Write (6,*) ' Full substitutions:'
      Write (6,*) ' -----------------------'
      Write (6,*)
      Do iElement=1,Num_Elem

         Found = .False.
         Do i = 1, mAtoms
            Found = Found .or. (PTab(iElement).eq.Element(i))
         End Do
*
*        Process if element found in the molecule.
*
         If (Found) Then
*
*           Loop over all natural isotopes of this element.
*
            Do k = 1, ElementList(iElement)%Nat
               IsoNum=ElementList(iElement)%Isotopes(k)%A
               Call Isotope(IsoNum,iElement,MassIso)
*
*              Substitute all instances.
*
               Changed = .False.
               Do i = 1, mAtoms
                  If (PTab(iElement).eq.Element(i)) Then
                     If (IsoNum.ne.nInt(Work(ipDefMass+i-1)/UtoAU))
     &                  Changed=.True.
                     Mass(i)=MassIso
                  End If
               End Do
*
               If (Changed) Then
                  call dcopy_(n**2,Temp(ipH),1,Temp(ipH2),1)
                  Write(6,*) 'Masses:'
                  Write(6,*) '======='
                  Write(6,'(20I4)') (nInt(mass(l)/UtoAU),l=1,mAtoms)
                  Write(6,*)
                  Write(6,*)
                  Write(6,*) 'Frequencies:'
                  Write(6,*) '============'
                  Call Freq_i(n,Temp(ipH2),mass,Temp(ip1),Temp(ip2),
     &                        Temp(ipVec),Temp(ipVal),ineg)
                  Call GFPrnt_i(Temp(ipVal),n)
               End If
*
            End Do
*
*           Restore the initial mass array
*
            Call dCopy_(mAtoms,Work(ipDefMass),1,Mass,1)
*
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Double substitutions
*
      If (Double) Then
      Write (6,*)
      Write (6,*) ' Double substitutions:'
      Write (6,*) ' -----------------------'
      Write (6,*)
*
      Do i = 1, mAtoms
         AtNum=iNuclearChargeFromSymbol(Element(i))
         If (lmass) Then
            Subs=0
            Do nAt=1,iMass
               If (Element(i).eq.cmass(nAt)) Then
                  Mass(i)=UtoAU*umass(nAt)
                  Subs=1
               End If
            End Do
         Else
            Subs=ElementList(AtNum)%Nat
         End If
*
*        All isotopes first atom
         dMass1 = Mass(i)
         Do k = 1, Subs
            If (.not.lmass) Then
               IsoNum=ElementList(AtNum)%Isotopes(k)%A
               If (IsoNum.eq.nInt(dMass1/UtoAU)) Cycle
               Call Isotope(IsoNum,AtNum,Mass(i))
            End If
*
*           Second atom
            Do j = i+1, mAtoms
               AtNum2=iNuclearChargeFromSymbol(Element(j))
               If (lmass) Then
                  Subs2=0
                  Do nAt=1,iMass
                     If (Element(j).eq.cmass(nAt)) Then
                        Mass(j)=UtoAU*umass(nAt)
                        Subs2=1
                     End If
                  End Do
               Else
                  Subs2=ElementList(AtNum2)%Nat
               End If
*
*              All isotopes for second atom
               dMass2 = Mass(j)
               Do l = 1, Subs2
                  If (.not.lmass) Then
                     IsoNum=ElementList(AtNum2)%Isotopes(l)%A
                     If (IsoNum.eq.nInt(dMass2/UtoAU)) Cycle
                     Call Isotope(IsoNum,AtNum2,Mass(j))
                  End If
*
                  call dcopy_(n**2,Temp(ipH),1,Temp(ipH2),1)
                  Write(6,*) 'Masses:'
                  Write(6,*) '======='
                  Write(6,'(20I4)')
     &                 (nInt(mass(m)/UtoAU),m=1,mAtoms)
                  Write(6,*)
                  Write(6,*)
                  Write(6,*) 'Frequencies:'
                  Write(6,*) '============'
                  call freq_i(n,Temp(ipH2),mass,Temp(ip1),
     &                Temp(ip2),Temp(ipVec),Temp(ipVal),ineg)
                  call GFPrnt_i(Temp(ipVal),n)

               End do
*              Put back the original mass
               Mass(j) = dMass2
            End do   ! End inner loop over atoms
*
         End Do
*        Put back the original mass
         Mass(i) = dMass1
      End Do ! End outer loop over atoms
*                                                                      *
************************************************************************
*                                                                      *
      End if
      Call GetMem('DefaultMass','FREE','REAL',ipDefMass,mAtoms)
*
*
      Return
      End
      Subroutine FREQ_i(nX,H,mass,Tmp1,Tmp2,EVec,EVal,iNeg)
      Implicit Real*8 (a-h,o-z)
      Real*8  Tmp1(nX,nX), Tmp2(nX,nX),H(nX,nX),
     &       EVec(2*nX,nX),
     &       EVal(2*nX),mass(*)

#include "constants2.fh"
*
      iprint=0
      call dcopy_(nX**2,0.0d0,0,Tmp1,1)
      Do i = 1, nX
       ii=(i-1)/3+1
       Do j=1,nX
            H(i,j) = H(i,j)/mass(ii)
       End Do
      End Do
*
*-----Compute eigenvectors and eigenfunctions
*
      nAux = 2 * nX
      iOpt=1
      islct=0
      If ( nX.gt.0 ) then
        Call Not_DGeEv(iOpt,H,nX,
     &             EVal,EVec,nX,
     &             iSlct,nX,Tmp2,nAux)
      End If
*
*-----Compute the harmonic frequencies
*
      iNeg=0
      Do 649 iHarm = 1, 2*nX, 2
         jHarm = (iHarm+1)/2
         temp = EVal(iHarm)
         If (temp.ge.0.0d0) Then
            EVal(jHarm) = Sqrt(temp)*autocm
         Else
            iNeg=iNeg+1
            EVal(jHarm) = -Sqrt(Abs(temp))*autocm
         End If
 649  Continue
      If (iPrint.ge.99) Call RecPrt('Converted EVal',' ',EVal,1,nX)
*
*-----Normalize
*
      Do iHarm = 1, nX
         r2=DDot_(nX,EVec(1,iHarm),2,EVec(1,iHarm),2)
         r2=1.0d0/Sqrt(r2)
         Call DScal_(nX,r2,EVec(1,iHarm),2)
      End Do
      If (iPrint.ge.99) Call RecPrt('Normalized EVec',' ',
     &                               EVec,nX*2,nX)
*
*-----Order, from low to high.
*
      Do iHarm = 1, nX-1
         Do jHarm = iHarm+1, nX
            If (EVal(jHarm).lt.EVal(iHarm)) Then
               rlow=EVal(iHarm)
               EVal(iHarm)=EVal(jHarm)
               EVal(jHarm)=rLow
               Call DSwap_(nX,EVec(1,iHarm),2,EVec(1,jHarm),2)
            End If
         End Do
      End Do
*
      Return
      End
      Subroutine GFPrnt_i(EVal,nDim)
      Implicit Real*8 (a-h,o-z)
      Real*8 EVal(nDim)
      Character*80 Format, Line*120
*
      LUt=6
      Inc = 6
      Do iHarm = 1, nDim, Inc
         Jnc=Min(Inc,nDim-iHarm+1)
         Write(Format,'(A,I3,A)') '(5X,A10,1x,',Jnc,'I10)'
         Write (LUt,Format) ' ',(i,i=iHarm,iHarm+Jnc-1)
         Write (LUt,*)
*
         Write(Format,'(A,I3,A)') '(A12,1x,',Jnc,'F10.2)'
         Line=' '
         Write (Line,Format) 'Freq.',(EVal(i),i=iHarm,iHarm+Jnc-1)
         Do i = 1, 120
            If (Line(i:i).eq.'-') Line(i:i)='i'
         End Do
         Write (LUt,'(5X,A)') Line
         Write (LUt,*)
*
         Write (LUt,*)
      End Do
      Return
      End
