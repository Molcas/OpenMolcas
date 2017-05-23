************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine RotateOrbB(CMO,Col,ipLbl,nComp,nBas,nOrb2Loc,
     &                      Maximisation,ThrRot,PctSkp,Debug)
C
C     Author: T.B. Pedersen
C
C     Purpose: rotate orbitals (Jacobi Sweeps) for Boys localisation.
C
      Implicit Real*8 (a-h,o-z)
      Real*8  CMO(nBas,*), Col(nOrb2Loc,2)
      Integer ipLbl(nComp)
      Logical Maximisation, Debug
#include "WrkSpc.fh"
#include "real.fh"

      Character*10 SecNam
      Parameter (SecNam = 'RotateOrbB')

      Character*80 Txt

      xDone = 0.0d0
      If (Debug) iCouple=0
      Do iMO1=1,nOrb2Loc-1
        Do iMO2=iMO1+1,nOrb2Loc

          If (Debug) Then
            iCouple = iCouple + 1
            Write(6,'(A9,I5)') 'Couple n:',iCouple
            Write(6,'(A9,I5)') '    MO1 :',iMO1
            Write(6,'(A9,I5)') '    MO2 :',iMO2
          End If

          iMO_s = iMO1
          iMO_t = iMO2

          Ast = 0.0d0
          Bst = 0.0d0
          Do iComp = 1,nComp
             ip0 = ipLbl(iComp) - 1
             iss = ip0 + nOrb2Loc*(iMO_s-1) + iMO_s
             itt = ip0 + nOrb2Loc*(iMO_t-1) + iMO_t
             ist = ip0 + nOrb2Loc*(iMO_t-1) + iMO_s
             Ast = Ast
     &           + Work(ist)**2
     &           - 2.5d-1*(Work(iss)-Work(itt))**2
             Bst = Bst
     &           + Work(ist)*(Work(iss)-Work(itt))
          End Do

          If ((Ast.eq.0.0d0).and.(Bst.eq.0.0d0)) Then
            cos4alpha=-1.0d0
            sin4alpha=0.0d0
          Else
            cos4alpha=-Ast/sqrt(Ast**2+Bst**2)
            sin4alpha= Bst/sqrt(Ast**2+Bst**2)
          End If
          Tst=abs(cos4alpha)-1.0d0
          If (Tst .gt. 0.0d0) Then
            If (Tst .gt. 1.0d-10) Then
              Write(Txt,'(A,D18.10)') 'Actual: cos4alpha = ',cos4alpha
              Call SysAbendMsg(SecNam,'-1.0d0 < cos4alpha < 1.0d0',
     &                         Txt)
            Else
               If (cos4alpha .lt. 0.0d0) Then
                  cos4alpha=-1.0d0
               Else
                  cos4alpha=1.0d0
               End If
            End If
          End If

          Alpha1=acos(cos4alpha)/4.0d0
          Alpha2=asin(sin4alpha)/4.0d0
          If (Alpha2.lt.0.0d0) Alpha1=Alpha2+PI
          Alpha=Alpha1
          If (.Not.Maximisation) Then
            Gamma_rot=Alpha-PI/4.0d0
          Else
            Gamma_rot=Alpha
          End If
          If (Debug) Then
            Write(6,'(A9,F10.5)') '   Ast :',Ast
            Write(6,'(A9,F10.5)') '   Bst :',Bst
            Write(6,'(A9,F10.5)') 'Alpha1 :',Alpha1
            Write(6,'(A9,F10.5)') 'Alpha2 :',Alpha2
            Write(6,'(A9,F10.5)') ' Gamma :',Gamma_rot
          End If

          Tsts = sin(Gamma_rot)
          Tstc = 1.0d0 - cos(Gamma_rot)
          If (abs(Tsts).gt.ThrRot .or. abs(Tstc).gt.ThrRot) Then
             Call Rot_st(CMO(1,iMO_s),CMO(1,iMO_t),nBas,Gamma_rot,Debug)
             Call UpdateB(Col,nOrb2Loc,ipLbl,nComp,Gamma_rot,
     &                    iMO_s,iMO_t,Debug)
             xDone = xDone + 1.0d0
          End If

         End Do
      End Do

      If (nOrb2Loc .gt. 1) Then
         xOrb2Loc = dble(nOrb2Loc)
         xTotal = xOrb2Loc*(xOrb2Loc-1.0d0)/2.0d0
         PctSkp = 1.0d2*(xTotal-xDone)/xTotal
      Else
         PctSkp = 0.0d0
      End If

      End
