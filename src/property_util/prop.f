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
      Subroutine Prop(Short,qplab,cen1,cen2,nIrrep,nBas,nTot,Occ,ThrSV,
     &                PrEl,PrNu,lpole,labs,tmat,temp)
************************************************************************
*                                                                      *
*     purpose: preprocessing of tables for different tensor            *
*              properties                                              *
*                                                                      *
*     Short           logical option for either Short (Short           *
*                     =.true., the total electronic contribution)      *
*                     long (Short=.false., orbital contributions)      *
*                     output                                           *
*     oplab           operator label as defined in SEWARD              *
*     cen1(1:3)       coordinates of centre no.1 (operator)            *
*     cen2(1:3)       coordinates of centre no.2 (gauge)               *
*     nIrrep          the number of irreducible representations        *
*     nBas            the number of functions in each representa-      *
*     (0:nIrrep-1)    tion                                             *
*     nTot            the total number of elements supplied for        *
*                     each component of the property tensor; equal     *
*                     either to 1 (total electronic and total nuc-     *
*                     lear contributions) or to the dimension of       *
*                     the basis set.                                   *
*     Occ(1:nTot)     occupation numbers for all eigenvectors,         *
*                     a dummy for Short outputs                        *
*     ThrSV           threshold for occupation numbers; If             *
*                     Occ(i).le.ThrSV the contribution will not        *
*                     be printed                                       *
*     PrEl(1:nTot,    matrix elements for all components 1,2,...,      *
*          1:maxlab)  maxlab, nTot entries for each component          *
*                     maxlab=(lpole+1)*(lpole+2)/2                     *
*     PrNu(1:maxlab)  nuclear contributions for each component         *
*     labs(1:maxlab)  labels for each component                        *
*     temp(1:maxlab)  auxiliary storage area                           *
*                                                                      *
*                                                                      *
*     lpole is the value for l in l-pole moments                       *
*     allocate integer storage area of the size appropriate            *
*     for the actual lpole value: size as below                        *
*                                                                      *
* 2000 Dept. of Chem. Phys., Univ. of Lund, Sweden                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "constants2.fh"
#include "real.fh"
#include "WrkSpc.fh"
      parameter (lmax=16)
      Character*2 lab2
      Character*3 lab3
      Character*4 lab4
      Character*5 lab5
      Character*8 qplab,oplab
      Character*(lmax) labs(1:(lpole+1)*(lpole+2)/2),lab
      Character*80 Line
      Logical Short
      Integer nBas(0:nIrrep-1)
      Integer Cho_X_GetTol
      External Cho_X_GetTol
      Real*8 cen1(1:3), cen2(1:3), Occ(1:nTot),
     &       PrEl(1:nTot,1:(lpole+1)*(lpole+2)/2),
     &       PrNu(1:(lpole+1)*(lpole+2)/2)
      Real*8 tmat(1:(lpole+1)*(lpole+2)/2,1:(lpole+1)*(lpole+2)/2),
     &       temp(1:(lpole+1)*(lpole+2)/2)
      Character*(lmax) labsAug(1:(lpole+1)*(lpole+2)/2+1)
      Real*8 PrElAug(1:nTot,1:(lpole+1)*(lpole+2)/2+1),
     &       PrNuAug(1:(lpole+1)*(lpole+2)/2+1)
      Real*8 Molecular_Charge
      Save Molecular_Charge
      Data Molecular_Charge/0.0D0/
      Logical StoreInfo, Reduce_Prt
      External Reduce_Prt
*                                                                      *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
*                                                                      *
************************************************************************
*                                                                      *
*     Some properties can be calculated at many, many centers,
*     and writing all values through Add_Info affects performance,
*     those will be deactivated and the summation is printed by the
*     calling function.
*     This only affects the call to Add_Info, otherwise the properties
*     are calculated and printed normally
*
      StoreInfo=.True.
*                                                                      *
************************************************************************
*                                                                      *
*
*     decipher the operator label:oplab
*
      oplab=qplab
      Call UpCase(oplab)
      Read (oplab,'(a4,a2,i2)') lab4,lab2,l
      If (lab4.eq.'MLTP') Then
         If (lpole.lt.0) Then
            Return
         Else
            Go To 200
         End If
      End If
      If (lab4(1:2).eq.'EF') Go To 300
      If (lab4.eq.'DMS ') Go To 500
      If (lab4(1:3).eq.'PAM') Go To 600
      If (lab4(1:3).eq.'CNT') Go To 700
*
*     invalid label supplied
*
      Write (6,'(//1x,a,a8,a//)')
     &' The label: *',oplab,'* is not a valid label ... stop'
      Call Quit_OnUserError()
*                                                                      *
************************************************************************
*                                                                      *
*     Multipole moment section ... generate labels for printing
*
  200 Continue
      If (lPole.gt.lMax) Then
         Write(6,*) 'Prop: lPole.gt.lMax'
         Write(6,*) 'lPole=',lPole
         Write(6,*) 'Increase lMax and recompile!'
         Call Abend()
      End If
      ilab=0
      Do ix=lpole,0,-1
         Do iy=lpole-ix,0,-1
            iz=lpole-ix-iy
            ilab=ilab+1
            lab='                '
            Do izz=lmax,lmax-iz+1,-1
               lab(izz:izz)='Z'
            End Do
            Do iyy=lmax-iz,lmax-iz-iy+1,-1
               lab(iyy:iyy)='Y'
            End Do
            Do ixx=lmax-iz-iy,lmax-lpole+1,-1
               lab(ixx:ixx)='X'
            End Do
            labs(ilab)=lab
         End Do
      End Do
*
      maxlab=ilab
      Call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)
*
*     Print cartesian moments
*
*----------------------------------------------------------------------*
      If (iPL.eq.2.and.Short) Then
*----------------------------------------------------------------------*

         Do i = 1, MaxLab
            Work(ipPrTot+i-1)=PrNu(i)-PrEl(1,i)
         End Do
*
*----    New style output
*
         Line=' '
         If (lPole.eq.0) Then
            Line='Charge (e):'
            Fact=1.0D0
         Else If (lPole.eq.1) Then
            Line='Dipole Moment (Debye):'
            Fact=Debye
         Else If (lPole.eq.2) Then
            Line='Quadrupole Moment (Debye*Ang):'
            Fact=Debye*Angstrom
         Else If (lPole.eq.3) Then
            Line='Octapole Moment (Debye*Ang**2):'
            Fact=Debye*Angstrom**2
         Else If (lPole.eq.4) Then
            Line='Hexadecapole Moment (Debye*Ang**3):'
            Fact=Debye*Angstrom**3
         Else
            Line=''
            If (lpole.le.9) Then
               Write (Line,'(I1)') lPole
               iSt=2
            Else
               Write (Line,'(I2)') lPole
               iSt=3
            End If
            Line(iSt:iSt+26)='th-pole Moment (Debye*Ang**'
            iSt=iSt+27
            If (lpole.le.10) Then
               Write(Line(iSt:iSt),'(I1)') lpole-1
               iSt=iSt+1
            Else
               Write(Line(iSt:iSt+1),'(I2)') lpole-1
               iSt=iSt+2
            End If
            Line(iSt:iSt+2)='):'
            Fact=Debye*Angstrom**(lPole-1)
         End If
         Write(6,'(6X,A)') Line(:mylen(Line))
         If (lpole.gt.0) Then
            Write (6,'(6X,A,3F10.4)') 'Origin of the operator (Ang)=',
     &            (cen1(i)*Angstrom,i=1,3)
         End If
         If (lPole.eq.0) Then
            tmp=Work(ipPrTot)
            Write (6,'(6X,A,A,F10.4)')
     &                       labs(1),'=',Work(ipPrTot  )*Fact
            Molecular_Charge=Work(ipPrTot  )*Fact
         Else If (lPole.eq.1) Then
            tmp=Sqrt(Work(ipPrTot  )**2
     &              +Work(ipPrTot+1)**2
     &              +Work(ipPrTot+2)**2)
            Write (6,'(4X,4(A,A,ES12.4))')
     &                       labs(1),'=',Work(ipPrTot  )*Fact,
     &                       labs(2),'=',Work(ipPrTot+1)*Fact,
     &                       labs(3),'=',Work(ipPrTot+2)*Fact,
     &            '           Total','=',tmp*Fact
            If (Abs(Molecular_Charge).gt.0.90D0) Then
               Write(6,'(6X,A)') 'Center of Charge (Ang)'
               X_Coor=Angstrom*(Work(ipPrTot  )/Molecular_Charge)
               Y_Coor=Angstrom*(Work(ipPrTot+1)/Molecular_Charge)
               Z_Coor=Angstrom*(Work(ipPrTot+2)/Molecular_Charge)
               Write (6,'(6X,3(A,A,F14.8))')
     &                       labs(1),'=',X_Coor,
     &                       labs(2),'=',Y_Coor,
     &                       labs(3),'=',Z_Coor
               Molecular_Charge=Zero
            End If
            Call Put_DArray('Dipole moment',Work(ipPrTot),3)
c            Call peek_iScalar('xml opened',isopen)
c            If(isopen.eq.1) Then
               Call xml_dDump('dipole','Dipole moment','Debye',1,
     &                        Work(ipPrTot),3,1)
c            End If
         Else If (lPole.ge.2) Then
            ip_=ipPrTot
            tmp=0.0D0
            Do i = 0, Maxlab-1
               tmp=Max(tmp,ABS(Work(ipPrTot+i)))
            End Do
            ip_=ipPrTot
            Do i = 1, maxlab, 4
               jMax=min(maxlab-i,3)
               Write (6,'(4X,4(A,A,ES12.4))')
     &          (labs(i+j),'=',Work(ip_+j)*Fact,j=0,jMax)
               ip_ = ip_ + 4
            End Do
         End If
         If (lpole.ge.2.and.lpole.le.4) Then
*
*           Transform cartesian moments to multipole moments
*
            inp=0
            Call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
            inp=1
            Call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
            Do i = 1, MaxLab
               Work(ipPrTot+i-1)=PrNu(i)-PrEl(1,i)
            End Do
*
            If (lPole.ge.3) Then
               Write(6,'(6X,A,I1,A)') 'In traceless form (Debye*Ang**',
     &              (lPole-1),')'
            Else
               Write(6,'(6X,A,I1,A)') 'In traceless form (Debye*Ang)'
            End If
            ip_=ipPrTot
            Do i = 1, maxlab, 4
               jMax=min(maxlab-i,3)
               Write (6,'(4X,4(A,A,ES12.4))')
     &          (labs(i+j),'=',Work(ip_+j)*Fact,j=0,jMax)
               ip_ = ip_ + 4
            End Do
*
         End If
*
*----------------------------------------------------------------------*
      Else If (iPL.ge.3.or.(.Not.Short.and.iPL.eq.2)) Then
*----------------------------------------------------------------------*
*
*----    Old style output
*
         If (lpole.eq.0) lab5=' 0-th'
         If (lpole.eq.1) lab5=' 1-st'
         If (lpole.eq.2) lab5=' 2-nd'
         If (lpole.eq.3) lab5=' 3-rd'
         If (lpole.gt.3) Then
           lab3='-th'
           Write (lab5,'(i2,a3)') lpole,lab3
         End If
         Write (6,'(//6x,a5,a,3(f12.8,a))') lab5,
     &         ' cartesian moments: origin at (',
     &         cen1(1),',',cen1(2),',',cen1(3),')'
         Write (6,'(6x,76(''-''))')
         sig=-One
         Call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,
     &              maxlab,labs,Work(ipPrTot),iPL,0)
         If (lpole.eq.1) Then
            Write (6,'(6x,76(''-''))')
            Write(6,'(6x,a,3f16.8,3x,a)')   'Total             ',
     &               (Work(ipPrTot+j)*Debye,j=0,2), 'Debye'
            Call Put_DArray('Dipole moment',Work(ipPrTot),3)
         End If
*
         If (lpole.gt.1.and.lpole.le.4) Then
           Write (6,'(//6x,a,i2,a,3(f12.8,a))') 'Cartesian ',lpole,
     &            '-pole moment: origin at (',
     &            cen1(1),',',cen1(2),',',cen1(3),')'
           Write (6,'(6x,76(''-''))')
*
*          Transform cartesian moments to multipole moments
*
           inp=0
           Call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
           inp=1
           Call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
*
*          Print l-pole cartesian moments
*
           sig=-One
           Call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,
     &                maxlab,labs,Work(ipPrTot),iPL,0)
         End If
*
*----------------------------------------------------------------------*
      Else
*----------------------------------------------------------------------*
         Do i = 1, MaxLab
            Work(ipPrTot+i-1)=PrNu(i)-PrEl(1,i)
         End Do
         If (lpole.eq.1)
     &      Call Put_DArray('Dipole moment',Work(ipPrTot),3)
         If (lpole.ge.2.and.lpole.le.4) Then
*
*           Transform cartesian moments to multipole moments
*
            inp=0
            Call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
            inp=1
            Call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
            Do i = 1, MaxLab
               Work(ipPrTot+i-1)=PrNu(i)-PrEl(1,i)
            End Do
         End If
*----------------------------------------------------------------------*
      End If ! iPL
*----------------------------------------------------------------------*
      Go To 999
*                                                                      *
************************************************************************
*                                                                      *
*     electric field section
*
  300 Continue
      Read (oplab,'(a3,i5)') lab3,icen
*
      If (lab3.eq.'EF0') Then
         MaxLab=1
*        set labels
         labs(1)='                '
*        Print cartesian components of the electric field gradient
*        tensor at the given centre
*
         If (iPL.ge.3.or.(.Not.Short.and.iPL.eq.2)) Then
           Write (6,'(//6x,a,i5,1x,a,3(f12.8,a))')
     &          ' Electric potential:  centre no.',icen,'(',
     &           cen1(1),',',cen1(2),',',cen1(3),')'
           Write (6,'(6x,72(''-''))')
         Else
           If (icen.eq.1) Then
             Write(6,'(//6X,A)') 'Electric potential:'
           End If
         End If
         sig=-One
      End If
      If (lab3.eq.'EF1') Then
         MaxLab=3
*        set labels
         labs(1)='               X'
         labs(2)='               Y'
         labs(3)='               Z'
*        Print cartesian components of the electric field vector
*        at the given centre
*
         If (iPL.ge.3.or.(.Not.Short.and.iPL.eq.2)) Then
           Write (6,'(//6x,a,i5,1x,a,3(f12.8,a))')
     &            ' Electric field:  centre no.',icen,'(',
     &             cen1(1),',',cen1(2),',',cen1(3),')'
           Write (6,'(6x,72(''-''))')
         Else
           If (icen.eq.1) Then
             Write(6,'(//6X,A)') 'Electric field:'
             Write(6,'(5X,6A16)') (labs(i),i=1,MaxLab)
           End If
         End If
         sig=+One
      End If
      If (lab3.eq.'EF2') Then
*        Seven components because we add r*r to the "normal" 6
         MaxLab=7
*        set labels
         labs(1)='  (2*XX-YY-ZZ)/2'
         labs(2)='          1.5*XY'
         labs(3)='          1.5*XZ'
         labs(4)='  (2*YY-ZZ-XX)/2'
         labs(5)='          1.5*YZ'
         labs(6)='  (2*ZZ-XX-YY)/2'
*        Actually, r*r is already stored as the 6th component
*        we now move it to the 7th using the "augmented" arrays
         Do j=1,5
           Do iOcc=1,nTot
             PrElAug(iOcc,j)=PrEl(iOcc,j)
           End Do
           PrNuAug(j)=PrNu(j)
           labsAug(j)=labs(j)
         End Do
         Do iOcc=1,nTot
*          Generate the actual value of the 6th component
           PrElAug(iOcc,6)=-PrEl(iOcc,1)-PrEl(iOcc,4)
           PrElAug(iOcc,7)=PrEl(iOcc,6)
         End Do
         PrNuAug(6)=PrNu(6)
         PrNuAug(7)=Zero
         labsAug(6)=labs(6)
         labsAug(7)='     RR=XX+YY+ZZ'
*
*        Print cartesian components of the electric field gradient
*        tensor at the given centre
*
         If (iPL.ge.3.or.(.Not.Short.and.iPL.eq.2)) Then
           Write (6,'(//6x,a,i5,1x,a,3(f12.8,a))')
     &          ' Electric field gradient:  centre no.',icen,'(',
     &           cen1(1),',',cen1(2),',',cen1(3),')'
           Write (6,'(6x,78(''-''))')
         Else
           If (icen.eq.1) Then
             Write(6,'(//6X,A)') 'Electric field gradient:'
             Write(6,'(5X,6A16)') (labsAug(i),i=1,MaxLab)
           End If
         End If
         sig=+One
      End If
      Call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)
*
*     Print the values using "augmented" arrays if needed
*
      If (Maxlab.eq.7) Then
        Call Prout(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrElAug,PrNuAug,
     &             maxlab,labsAug,Work(ipPrTot),iPL,icen)
        MaxLab=6 ! Reset so call to Add_Info is correct!
*
      Else
        Call Prout(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,
     &             maxlab,labs,Work(ipPrTot),iPL,icen)
      End If
*
      If (lab3.eq.'EF2') Then
         Tmp = Work(ipPrTot+2)
         Work(ipPrTot+2)=Work(ipPrTot+3)
         Work(ipPrTot+3)=Tmp
         If (iPL.ge.3.or.(.Not.Short.and.iPL.eq.2))
     &      Call Print_EigenValues(Work(ipPrTot),3)
         Tmp = Work(ipPrTot+2)
         Work(ipPrTot+2)=Work(ipPrTot+3)
         Work(ipPrTot+3)=Tmp
      End If
*
*     do not write the different electric field components through Add_Info
      StoreInfo=.False.
      Go To 999
*                                                                      *
************************************************************************
*                                                                      *
*     diamagnetic shielding section
*
  500 Continue
      Read (oplab,'(a4,i2,i2)') lab4,icen2,icen1
*
      maxlab=9
      Call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)
*
*     set labels
      labs(1)='   (YO*YG+ZO*ZG)'
      labs(2)='         (XO*YG)'
      labs(3)='         (XO*ZG)'
      labs(4)='         (YO*XG)'
      labs(5)='   (XO*XG+ZO*ZG)'
      labs(6)='         (YO*ZG)'
      labs(7)='         (ZO*XG)'
      labs(8)='         (ZO*YG)'
      labs(9)='   (XO*XG+YO*YG)'
*
*     Print cartesian components of the diamagnetic shielding
*     tensor at the given centre (O=cen1) with the gauge ori-
*     gin a centre G=cen2
*
      Write (6,'(//6x,a,i3,1x,a,3(f12.8,a)/1x,a,3(f12.8,a))')
     &       ' Diamagnetic shielding:   centre no.',icen1,'(',
     &       cen1(1),',',cen1(2),',',cen1(3),')',
     &       '                       gauge origin at (',
     &       cen2(1),',',cen2(2),',',cen2(3),')'
      Write (6,'(6x,78(''-''))')
      sig=+One
      Call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,
     &           labs,Work(ipPrTot),iPL,0)
*                                                                      *
************************************************************************
*                                                                      *
*     PAM section ... generate labels for printing
*
  600 Continue
      ilab=0
      Do ix=lpole,0,-1
         Do iy=lpole-ix,0,-1
            iz=lpole-ix-iy
            ilab=ilab+1
            lab='                '
            Do izz=lmax,lmax-iz+1,-1
               lab(izz:izz)='Z'
            End Do
            Do iyy=lmax-iz,lmax-iz-iy+1,-1
               lab(iyy:iyy)='Y'
            End Do
            Do ixx=lmax-iz-iy,lmax-lpole+1,-1
               lab(ixx:ixx)='X'
            End Do
            labs(ilab)=lab
         End Do
      End Do
*
      maxlab=ilab
      Call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)
*
*     Print cartesian moments
*
      If (lpole.eq.0) lab5=' 0-th'
      If (lpole.eq.1) lab5=' 1-st'
      If (lpole.eq.2) lab5=' 2-nd'
      If (lpole.eq.3) lab5=' 3-rd'
      If (lpole.gt.3) Then
        lab3='-th'
        Write (lab5,'(i2,a3)') lpole,lab3
      End If
      Write (6,'(//6x,a5,a,3(f12.8,a))') lab5,
     &      ' cartesian moments: origin at (',
     &      cen1(1),',',cen1(2),',',cen1(3),')'
      Write (6,'(6x,76(''-''))')
      sig=One
      Call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,
     &           labs,Work(ipPrTot),iPL,0)
*
      If (lpole.gt.1.and.lpole.le.4) Then
        Write (6,'(//6x,a,i2,a,3(f12.8,a))') 'Cartesian ',lpole,
     &         '-pole moment: origin at (',
     &         cen1(1),',',cen1(2),',',cen1(3),')'
        Write (6,'(6x,76(''-''))')
*
*
*
        inp=0
        Call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
        inp=1
        Call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
*
*       Print 0-pole PAM integrals (sig=1 in opposite multipole moments)
*
        sig=One
        Call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,
     &             maxlab,labs,Work(ipPrTot),iPL,0)
      End If
*      Go To 999
*                                                                      *
************************************************************************
*                                                                      *
*     Contact term section
*
  700 Continue
      Read (oplab,'(a3,i5)') lab3,icen
*
      maxlab=1
      Call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)
*
*     set labels
      labs(1)='      Delta(R-C)'
*
      IF (iPL.ge.3.or.(.Not.Short.and.iPL.eq.2)) Then
        Write (6,'(//6x,a,i5,1x,a,3(f12.8,a))')
     &       ' Contact term:  centre no.',icen,'(',
     &        cen1(1),',',cen1(2),',',cen1(3),')'
        Write (6,'(6x,78(''-''))')
      Else
        If (icen.eq.1) Then
          Write(6,'(//6X,A)') 'Contact term:'
          Write(6,'(5X,6A16)') (labs(i),i=1,MaxLab)
        End If
      End If
      sig=+One
      Call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,
     &           labs,Work(ipPrTot),iPL,icen)
*
*     do not write the contact term through Add_Info
      StoreInfo=.False.
*                                                                      *
************************************************************************
*                                                                      *
  999 Continue
      iTol=5
      iTol_E0=8
      iTol_E1=Cho_X_GetTol(iTol_E0)
      iTol = INT( DBLE(iTol) * DBLE(iTol_E1)/DBLE(iTol_E0))
      If (StoreInfo) Call Add_Info(OpLab,Work(ipPrTot),maxlab,iTol)
      Call GetMem('PrTot','Free','Real',ipPrTot,maxlab)
      Return
      End
