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
* Copyright (C) Yannick Carissan                                       *
*               2005, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RotateOrb(cMO,PACol,nBasis,nAtoms,PA,
     &                     Maximisation,nOrb2loc,Name,
     &                     nBas_per_Atom,nBas_Start,ThrRot,PctSkp,
     &                     Debug)
c
c     Author: Yannick Carissan.
c
c     Modifications:
c        - October 6, 2005 (Thomas Bondo Pedersen):
c          Array PACol introduced in argument list.
c
      Implicit Real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
#include "Molcas.fh"
      Real*8 cMO(nBasis,*), PACol(nOrb2Loc,2)
      Integer nBas_per_Atom(*),nBas_Start(*)
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Logical Maximisation, Debug
      Character*80 Txt
      Character*(LENIN8) Name(*),PALbl
c
      xDone = 0.0d0
      If (Debug) Then
         Write(6,*) 'RotateOrb[Debug]: nBas_per_Atom: ',
     &              (nBas_per_Atom(i),i=1,nAtoms)
         iCouple=0
      End If
      Do iMO1=1,nOrb2Loc-1
        Do iMO2=iMO1+1,nOrb2Loc
c
          If (Debug) Then
            iCouple = iCouple + 1
            Write(6,'(a9,i5)') 'Couple n:',iCouple
            Write(6,'(a9,i5)') '    MO1 :',iMO1
            Write(6,'(a9,i5)') '    MO2 :',iMO2
          End If
c
          iMO_s=iMO1
          iMO_t=iMO2
          SumA=Zero
          SumB=Zero
          Do iAt=1,nAtoms
            PAst=PA(iMO_t,iMO_s,iAt)
            PAss=PA(iMO_s,iMO_s,iAt)
            PAtt=PA(iMO_t,iMO_t,iAt)
            If (Debug) Then
              Write(6,*) 'In RotateOrb'
              Write(6,*) '------------'
              PALbl='PA__'//Name(nBas_Start(iAt))(1:LENIN)
              Call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
              Write(6,*) '**************************'
              Write(6,*) 'A :',iAt
              Write(6,*) '<',iMO_s,'|PA|',iMO_t,'> = ',PAst
              Write(6,*) '<',iMO_s,'|PA|',iMO_s,'> = ',PAss
              Write(6,*) '<',iMO_t,'|PA|',iMO_t,'> = ',PAtt
              Write(6,*) '**************************'
            End If
            SumA=SumA+PAst**2-0.25d0*(PAss-PAtt)**2
            SumB=SumB+PAst*(PAss-PAtt)
          End Do
          Ast=SumA
          Bst=SumB
c
          If ((Ast.eq.Zero).and.(Bst.eq.Zero)) Then
            cos4alpha=-One
            sin4alpha=Zero
          Else
            cos4alpha=-Ast/sqrt(Ast**2+Bst**2)
            sin4alpha= Bst/sqrt(Ast**2+Bst**2)
          End If
          Tst=abs(cos4alpha)-1.0d0
          If (Tst .gt. 0.0d0) Then
            If (Tst .gt. 1.0d-10) Then
              Write(Txt,'(A,D18.10)') 'Actual: cos4alpha = ',cos4alpha
              Call SysAbendMsg('RotateOrb','-1.0d0 < cos4alpha < 1.0d0',
     &                         Txt)
            Else
               If (cos4alpha .lt. 0.0d0) Then
                  cos4alpha=-1.0d0
               Else
                  cos4alpha=1.0d0
               End If
            End If
          End If
c
c-------- On choisit le cos car Alpha IN [0;PI/2]
c
          Alpha1=acos(cos4alpha)/4.0d0
          Alpha2=asin(sin4alpha)/4.0d0
          If (Alpha2.lt.Zero) Alpha1=Alpha2+PI
          Alpha=Alpha1
          If (.Not.Maximisation) Then
            Gamma_rot=Alpha-PI/4.0d0
          Else
            Gamma_rot=Alpha
          End If
          If (Debug) Then
            Write(6,'(a9,f10.5)') '   Ast :',Ast
            Write(6,'(a9,f10.5)') '   Bst :',Bst
            Write(6,'(a9,f10.5)') 'Alpha1 :',Alpha1
            Write(6,'(a9,f10.5)') 'Alpha2 :',Alpha2
            Write(6,'(a9,f10.5)') ' Gamma :',Gamma_rot
          End If
c
          Tsts = sin(Gamma_rot)
          Tstc = 1.0d0 - cos(Gamma_rot)
          If (abs(Tsts).gt.ThrRot .or. abs(Tstc).gt.ThrRot) Then
             Call Rot_st(cMO(1,iMO_s),cMO(1,iMO_t),nBasis,Gamma_rot,
     &                   Debug)
             Call UpdateP(PACol,Name,nBas_Start,
     &                    nOrb2Loc,nAtoms,PA,Gamma_rot,
     &                    iMO_s,iMO_t,Debug)
             xDone = xDone + 1.0d0
          End If
c
          If (Debug) Then
            Call RecPrt('MO after rotation',' ',cMO,nBasis,nBasis)
          End If
c
        End Do !iMO2
      End Do !iMO1
c
      If (nOrb2Loc .gt. 1) Then
         xOrb2Loc = dble(nOrb2Loc)
         xTotal = xOrb2Loc*(xOrb2Loc-1.0d0)/2.0d0
         PctSkp = 1.0d2*(xTotal-xDone)/xTotal
      Else
         PctSkp = 0.0d0
      End If
c
      Return
      End
