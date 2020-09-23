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
      Subroutine Print_Local(rMP,nij,nElem,Coor,nAtoms,C_o_C,Q_Nuc,lMax,
     &                       Lbl_Center,rMPq,EC,Polar,NoField,Temp,
     &                       xrMP,xxRMP,xnrMP,iANr,nOcOb,
     &                       Energy_Without_FFPT,ip_Ene_Occ,
     &                       MpProp_Level,Bond_Threshold,XHole,
     &                       XHoleLoc,D2, ChPol, ChPolBB,LIonize)
      use Real_Spherical
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "status.fh"
      Real*8 rMP(nij,nElem), Coor(3,nAtoms), A(3), C_o_C(3),
     &       Q_Nuc(nAtoms), rMPq(nElem), EC(3,nij), Polar(6,nij),
     &       Polar_M(6), Temp(nij), CX(3), xrMP(nij,nElem),
     &       xxrMP(nij,nElem), xnrMP(nij,nElem), XHoleLoc(nij),
     &       ChPol(6,nij), ChPolBB(6,nij)
      Integer iANr(nAtoms)
      Character*(LENIN) Lbl_Center(nAtoms)
      Character*80 Banner_Line(2)
      Logical NoField, Center_OK, Check_Bond, get_BasisType, XHole
      Logical LIonize
*
      Real*8 Sites
      Real*8 TP
      Real*8 TP_Q
      Real*8 LI, LA
      Real*8 LI_TOT
      Real*8 LI_MIN
*
      Sites = 0.0
      LI_TOT = 0.0
      LI_MIN = 99.0
*                                                                      *
************************************************************************
*                                                                      *
      Center_OK = .True.
      Call Sphere(lMax)
*                                                                      *
************************************************************************
*                                                                      *
*     Make copy of rMP so we keep a copy of the
*     electronic contributions to the multipoles.
*
      iDim = nij*nElem
      Call dCopy_(iDim,rMP,1,xrMP,1)
*                                                                      *
************************************************************************
*                                                                      *
*     Put info to check file.
*
C     Call RecPrt('Print_Local: rMP',' ',rMP,nij,nElem)
      Do iElem = 1, nElem
         Call Add_Info('LoProp MP',rMP(1,iElem),nij,4)
      End Do
      If (.Not.NoField) Then
         Do iPol = 1, 6
            call dcopy_(nij,Polar(iPol,1),6,Temp,1)
            Call Add_Info('LoProp a',Temp,nAtoms*(nAtoms+1)/2,4)
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Binom()
*                                                                      *
************************************************************************
*                                                                      *
*     Set up array for nuclear contributions.
*
      Call dCopy_(nij*nElem,[Zero],0,xnrMP,1)
      ij = 0
      Do iAtom = 1, nAtoms
         Do jAtom = 1, iAtom
            ij = ij + 1
            If (iAtom.eq.jAtom) Then
               xnrMP(ij,1)=Q_nuc(iAtom)
               CX(1)=Coor(1,iAtom)
               CX(2)=Coor(2,iAtom)
               CX(3)=Coor(3,iAtom)
               Call ReExpand(xnrMP,nij,nElem,CX,EC(1,ij),ij,lMax)
            End If
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over all domains
*
*
      Write (6,*)
      Banner_Line(1)='The Localized properties'
      Call Banner (Banner_Line(1),1,80)
      Write(6,*)

      nReal_Centers = 0
      ij = 0
      Do iAtom = 1, nAtoms
         Do jAtom = 1, iAtom
            ij = ij + 1
            If (iAtom.eq.jAtom) Then
*--------------Atomic domain
               Charge_center=Q_nuc(iAtom)
               CX(1)=Coor(1,iAtom)
               CX(2)=Coor(2,iAtom)
               CX(3)=Coor(3,iAtom)
               Center_OK = .True.
            Else
*--------------Bond domain
               Charge_center=Zero
               CX(1)=(Coor(1,iAtom)+Coor(1,jAtom))*Half
               CX(2)=(Coor(2,iAtom)+Coor(2,jAtom))*Half
               CX(3)=(Coor(3,iAtom)+Coor(3,jAtom))*Half
               Center_OK = Check_Bond(Coor(1,iAtom),Coor(1,jAtom),
     &                           iANr(iAtom),iANr(jAtom),Bond_Threshold)
            End If
            call dcopy_(3,EC(1,ij),1,A,1)
*
*---------- Print out values for each domain.
*
            If (Center_OK) Then
            nReal_Centers = nReal_Centers + 1
            Write (6,*)
            Write (6,*)
            If (iAtom.eq.jAtom) Then
               Write (6,*)       '===================='
               Write (6,'(A,A)') ' ATOMIC DOMAIN: ', Lbl_Center(iAtom)
               Write (6,*)       '===================='
            Else
               Write (6,*)       '========================='
               Write (6,'(A,A,A,A)') ' BOND DOMAIN: ',
     &                         Lbl_Center(iAtom), ',', Lbl_Center(jAtom)
               Write (6,*)       '========================='
            End If
            Write (6,'(A,3F12.8,A)') ' Domain center:  :',CX,' / bohr'
            Write (6,'(A,3F12.8,A)') ' Expansion center:', A,' / bohr'
            Write (6,'(A, F12.8  )') ' Total charge    :',
     &                               Charge_center+rMP(ij,1)
**
            If (LIonize) Then
               If (.Not.NoField) Then
                  TP  = (Polar(1,ij)+Polar(3,ij)+Polar(6,ij))/Three
                  TP_Q = Charge_center+rMP(ij,1)
                  Sites = Sites + 1.0
                  If(TP.gt.0.0) Then
                     LI = 0.623*(TP_Q+1.742)*(1/TP)**(1.0/3.0)
                     LA = LI*(TP_Q+0.258)/(TP_Q+1.742)
                     Write (6,*)
                     Write (6,'(A)') 'Local Ionization Energy'
                     Write (6,'(10F12.8)') LI
                     Write (6,*)
                     If (LI.lt.LI_MIN) LI_MIN = LI
                     LI_TOT = LI_TOT + LI
                     Write (6,'(A)') 'Local Electron Affinity'
                     Write (6,'(10F12.8)') LA
                     Write (6,'(A)') 'Local Electronegativity'
                     Write (6,'(10F12.8)') (LI+LA)/2.0
                     Write (6,'(A)') 'Local Chemical Hardness'
                     Write (6,'(10F12.8)') (LI-LA)/2.0

                  Else
                     Write(6,'(A)') 'Negative isotropic polarizability!'
                     Write(6,'(A)') 'The local ionization energy will'
                     Write(6,'(A)') 'not be computed for this site'
                  EndIf
               EndIf
            End If
**
            Write (6,*)
            Write (6,'(A)') ' Electronic multipole moments:'
            iStrt = 1


            Do l = 0, lMax
               iEnd = iStrt + (l+1)*(l+2)/2 - 1
               If (l .eq. 0) Write (6,'(A)') 'Electronic Charge'
               If (l .eq. 1) Write (6,'(A)') 'Electronic Dipole'
               If (l .eq. 2) Write (6,'(A)') 'Electronic Quadrupole'
               If (l .eq. 3) Write (6,'(A)') 'Electronic Octupole'
               If (l .eq. 4) Write (6,'(A)') 'Electronic Hexadecapole'
               If (l .ge. 5) Write (6,'(A,I2)')
     &                  'Electronic multipole moment with l = ',l
               Write (6,'(10F12.8)') (rMP(ij,iElem),iElem=iStrt,iEnd)
c****************************************************
               If (iAtom .eq. jAtom .AND. l .ge. 1) Then
                  Write (6,'(A)') '... with nuclear contribution'
                  Write (6,'(10F12.8)')
     &               (rMP(ij,iElem)+xnrMP(ij,iElem),iElem=iStrt,iEnd)
               End If
c****************************************************
               iStrt = iEnd + 1
               Write (6,*)
            End Do
            If (lMax.ge.1) Then
               Dip_Tot=Sqrt(rMP(ij,2)**2+rMP(ij,3)**2+rMP(ij,4)**2)
               Write (6,*)
               Write (6,'(A,F12.8)')'Dipole magnitude:',Dip_Tot
               Write (6,*)
            End If
            If (.Not.NoField) Then
               Write (6,*)
               Call TriPrt(
     &         'Symmetrized Local Polarizability Tensor',' ',
     &         Polar(1,ij),3)
               Pol_Tot=(Polar(1,ij)+Polar(3,ij)+Polar(6,ij))/Three
               Write (6,*)
               Write (6,'(A,F12.8)')'Isotropic Polarizability:',Pol_Tot
               Write (6,*)
            End If
            If (XHole) then
              x2Dip=XHoleLoc(ij)
              Write(6,'(A,F12.8)')'Exchange hole second-moment:',x2Dip
              Write(6,*)
            Endif
            End If ! Center_OK
*
*debug
C           Call Allocate_Work(ip_EVec,9)
C           Call Allocate_Work(ip_EVal,6)
C           call dcopy_(9,[Zero],0,Work(ip_EVec),1)
C           call dcopy_(3,[One],0,Work(ip_EVec),4)
C           call dcopy_(6,Polar(1,ij),1,Work(ip_EVal),1)
C           Call Jacob(Work(ip_EVal),Work(ip_EVec),3,3)
C           Call TriPrt('EVal',' ',Work(ip_EVal),3)
C           Call RecPrt('EVec',' ',Work(ip_EVec),3,3)
C           Call Free_Work(ip_EVal)
C           Call Free_Work(ip_EVec)
*debug
*           Reexpand all multipole moments to the center of mass
*
            Call ReExpand(rMP,nij,nElem,A,C_o_C,ij,lMax)
C??            Call ReExpand(xrMP,nij,nElem,A,C_o_C,ij,lMax)
            Call ReExpand(xnrMP,nij,nElem,EC(1,ij),C_o_C,ij,lMax)
*
         End Do
      End Do
*     Write the charge capacitances for the bonds
      If (.Not.NoField) Then
         Write(6,*) '=== Charge capacitance for bonds ==='
         ij = 0
         Do jAtom = 1, nAtoms
            Do iAtom = jAtom+1, nAtoms
               ij = iAtom*(iAtom-1)/2 + jAtom
               If (iAtom.eq.jAtom) Then
               Else
                  CX(1)=(Coor(1,iAtom)-Coor(1,jAtom))
                  CX(2)=(Coor(2,iAtom)-Coor(2,jAtom))
                  CX(3)=(Coor(3,iAtom)-Coor(3,jAtom))
                  CRX = ChPolBB(1,ij)*CX(1)
     &                 +ChPolBB(2,ij)*CX(2)
     &                 +ChPolBB(4,ij)*CX(3)

                  CRY = ChPolBB(2,ij)*CX(1)
     &                 +ChPolBB(3,ij)*CX(2)
     &                 +ChPolBB(5,ij)*CX(3)

                  CRZ = ChPolBB(4,ij)*CX(1)
     &                 +ChPolBB(5,ij)*CX(2)
     &                 +ChPolBB(6,ij)*CX(3)

                  CRN = SQRT(CRX*CRX+
     &                 CRY*CRY+
     &                 CRZ*CRZ)

                  DR = SQRT(CX(1)*CX(1)+
     &                 CX(2)*CX(2)+
     &                 CX(3)*CX(3))
                  If ( 0.001.lt.(CRN/(DR*DR*DR))) Then
                     Write(6,' (A, A, F12.8)') Lbl_Center(iAtom),
     &                                        Lbl_Center(jAtom),
     &                                        CRN/(DR*DR*DR)
                  End If
               End If
            End Do
         End Do
         Write(6,*) '=== =========================== ==='
         Write(6,*)
      End If
C      *** Local ionzations are implemented by A. Holt
      If (LIonize) Then
         If (.Not.NoField) Then
            Write (6,'(A)') '===   Local Ionization energy   ==='
            Write (6,'(A)') 'The local ionization energies are'
            Write (6,'(A)') 'computed using the expression for atoms'
            Write (6,'(A)') 'found in; J.Phys.Chem. 1996, 100,4828'
            Write (6,'(A)')
            Write (6,'(A)') 'Average local ionization energy'
            Write (6,'(10F12.8)') LI_TOT/Sites*27.2114
            Write (6,'(A)') 'Lowest local ionization energy (eV)'
            Write (6,'(10F12.8)') LI_MIN*27.2114
            Write (6,'(A)') 'HOMO energy (absolute value, eV)'
            Write (6,'(10F12.8)') ABS(Work(ip_Ene_Occ+nOcOb-1))*27.2114
            Write(6,'(A)') '=== =========================== ==='
            Write (6,'(A)')
         EndIf
      EndIf

*
*     Write out the molecular properties
*
      Write (6,*)
      Write (6,*)
      Write (6,*)
      Banner_Line(1)='The Molecular Multipole Moments'
      Call Banner (Banner_Line(1),1,80)
      Write (6,'(A,3F12.8,A)') ' Expansion center:',C_o_C,' / bohr'
      Write (6,*)
      Write (6,*)
      iElem = 0
      Do l = 0, lMax
         Write (6,*)
         Write (6,'(A,I1)') 'l=',l
         Write (6,*)
         Write (6,'(A)')
     &         'xyz    Nuclear        Electronic     Molecular   '
         Write (6,*)
         Do ix = l, 0, -1
         Do iy = l-ix,0,-1
            iz = l-ix-iy
            iElem = iElem + 1
            rMP_Tot_Electronic=DDot_(nij,[One],0,rMP(1,iElem),1)
            rMP_Tot_Nuclear=rMPq(iElem)
            rMP_Tot = rMP_Tot_Nuclear + rMP_Tot_Electronic
            Write (6,'(3I1,3F16.8)') ix,iy,iz,rMP_Tot_Nuclear,
     &                               rMP_Tot_Electronic,rMP_Tot
         End Do
         End Do
      End Do
      If (.Not.NoField) Then
         Do iPol = 1, 6
            Polar_M(iPol)=DDot_(nij,[One],0,Polar(iPol,1),6)
         End Do
         Call TriPrt(
     &         'Molecular Polarizability Tensor',' ',Polar_M,3)
*debug
C           Call Allocate_Work(ip_EVec,9)
C           Call Allocate_Work(ip_EVal,6)
C           call dcopy_(9,Zero,0,Work(ip_EVec),1)
C           call dcopy_(3,One,0,Work(ip_EVec),4)
C           call dcopy_(6,Polar_M,1,Work(ip_EVal),1)
C           Call Jacob(Work(ip_EVal),Work(ip_EVec),3,3)
C           Call TriPrt('EVal',' ',Work(ip_EVal),3)
C           Call RecPrt('EVec',' ',Work(ip_EVec),3,3)
C           Call Free_Work(ip_EVal)
C           Call Free_Work(ip_EVec)
*debug
      End If
      If(XHole) then
        Write(6,*)
        Write(6,'(A,F12.8)')'Molecular exchange hole second moment:',D2
        Write(6,*)
      Endif
*
*---------- Generate mpprop file
*
      Call Print_MPPROP(rMP,xrMP,xnrMP,nij,nElem,lMax,EC,Polar,
     &                  Lbl_Center,nAtoms,iANr,NoField,C_o_C,Coor,
     &                  nOcOb,Energy_Without_FFPT,ip_Ene_Occ,
     &                  MpProp_Level,Bond_Threshold,nReal_Centers)
*
      Call Allocate_Work(iScratch_New,nij*(2*lMax+1))
      Call Allocate_Work(iScratch_Org,nij*(2*lMax+1))
*
*---------- Compare the molecular multipole moments to the ones arrising from truncation
*           xrMP contains the original multipole moments (electr. + nuclear)
*           xxrMP contains the approximate multipole moments.
*
      Call dCopy_(iDim,rMP,1,xrMP,1)
      call daxpy_(iDim,One,xnrMP,1,xrMP,1)
      Call dCopy_(iDim,xrMP,1,xxrMP,1)

      iPrint = 1
      Do l = lMax-1, 0, -1
         Call CutOff_Error(l,lMax,xrMP,xxrMP,nij,EC,C_o_C,nElem,
     &                     Work(iScratch_New),Work(iScratch_Org),
     &                     nAtoms,iPrint,rms)
      End Do

      Call Free_Work(iScratch_New)
      Call Free_Work(iScratch_Org)
      Call Sphere_Free()
*
*------ Print warning if non-ANO basis sets have been used
*
      If (.NOT. get_BasisType('ANO')) Then
         Write(6,*)
         Write(6,*) 'WARNING: The calculation were performed with at '
         Write(6,*) '         least one non-ANO basis! The results'
         Write(6,*) '         might therefore be erroneous.'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(ChPol)
      End
