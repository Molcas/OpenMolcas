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
      Subroutine LA_Morok(nAtom,CorG,iMode)
      Implicit Real*8 (a-h,o-z)
      Real*8 CorG(3,nAtom)
c
c     Morokuma's scaling scheme:
c       k = (q_LA - q_QM)/(q_MM - q_QM), k = constant
c
c     iMode = 1 => the LA gradient is distributed
c     on the frontier QM and MM atoms according to:
c     dE/dq_QM = dE/dq_QM + dE/dq_LA * dq_LA/dq_QM
c     dE/dq_MM = dE/dq_MM + dE/dq_LA * dq_LA/dq_MM
c
c     iMode = 2 => the LA position is updated
c     q_LA = q_QM + k * (q_MM - q_QM)
c
#include "espf.fh"
      Include 'stdalloc.fh'
*
      Logical Exist,Exist2,isOkLA,isOkMM,isOkQM,lMorok
      Logical DoTinker,DoGromacs
      Character*10 ESPFKey
      Character*180 Line
      Character*180 Get_Ln
      Character*256 Message
      External Get_Ln
*
      Integer, Dimension(:), Allocatable :: AT,GroToMol
      Integer, Dimension(:,:), Allocatable :: DefLA
      Real*8, Dimension(:), Allocatable :: FactLA
*
*
*define _DEBUG_
*
      iPL = iPL_espf()
      lMorok = .False.
      DoTinker = .False.
      DoGromacs = .False.
      Call F_Inquire('ESPF.DATA',Exist)
      If (Exist) Then
         IPotFl = IsFreeUnit(1)
         Call Molcas_Open(IPotFl,'ESPF.DATA')
10       Line = Get_Ln(IPotFl)
         ESPFKey = Line(1:10)
         If (ESPFKey.eq.'LA_MOROK  ') Then
            lMorok = .True.
         Else If (ESPFKey.eq.'TINKER    ') Then
            DoTinker = .True.
         Else If (ESPFKey.eq.'GROMACS   ') Then
            DoGromacs = .True.
         Else If (ESPFKey.eq.'ENDOFESPF ') Then
            Goto 11
         EndIf
         Goto 10
11       Close (IPotFl)
      End If
      If (.not.lMorok) Return
#ifdef _DEBUG_
      iPL = 4
      Call RecPrt('LA_Morok: coord or grad:',' ',CorG,3,nAtom)
#endif
c
c Tinker part
c
      Exist = .False.
      If (DoTinker) Then
         Call F_Inquire('QMMM',Exist)
      End If
      If(Exist) Then
         ITkQMMM=IsFreeUnit(25)
         Call Molcas_Open(ITkQMMM,'QMMM')
         Line = ' '
         Do While (Index(Line,'TheEnd ') .eq. 0)
            Line=Get_Ln(ITkQMMM)
            If (Index(Line,'LAH').ne.0) Then
               Call Get_I1(2,iLA)
               Call Get_I1(3,iMM)
               Call Get_I1(4,iQM)
               Call Get_F1(5,Fact)
               If (iMM.lt.1 .or. iQM.lt.1) Then
                  Write (6,*) 'LA_Morok: link atoms badly defined'
                  Write (6,*) '          check each LA connectivity'
                  Call Quit_OnUserError()
               End If
#ifdef _DEBUG_
               Write(6,*)
               Write(6,*) 'LA_Morok: LAH ',iLA,' between ',iQM,
     &                                             ' and ',iMM
               Write(6,*) '          scaling factor : ',Fact
#endif
               If (iMode .eq. 1) Then
                  If (iPL.ge.2) Write(6,*) 'LA_Morok: scaling gradients'
                  Do iXYZ = 1, 3
                     CorG(iXYZ,iQM)=CorG(iXYZ,iQM)+
     &                              CorG(iXYZ,iLA)*(One-Fact)
                     CorG(iXYZ,iMM)=CorG(iXYZ,iMM)+
     &                              CorG(iXYZ,iLA)*Fact
                     CorG(iXYZ,iLA)=Zero
                  End Do
               Else If (iMode .eq. 2) Then
                  If (iPL.ge.2) Write(6,*)'LA_Morok: updating positions'
                  Do iXYZ = 1, 3
                     CorG(iXYZ,iLA) = CorG(iXYZ,iQM) + (CorG(iXYZ,iMM)
     &                              - Corg(iXYZ,iQM))*Fact
                  End Do
               Else
                  Write (6,*) 'LA_Morok: wrong iMode'
                  Call Quit_OnUserError()
               End If
            End If
         End Do
         Close (ITkQMMM)
      End If
*
* Gromacs part
*
      Exist = .False.
      If (DoGromacs) Then
         Call Qpg_iArray('LA Def',Exist,nLink)
      End If
      If (Exist) Then
         nLink = nLink/3
         Call mma_allocate(DefLA,3,nLink)
         Call mma_allocate(FactLA,nLink)
         Call Get_iArray('LA Def',DefLA,3*nLink)
         Call Get_dArray('LA Fact',FactLA,nLink)
* Check for consistency
         Call Qpg_iArray('Atom Types',Exist2,nTot)
         If (.Not.Exist2) Then
            Message = 'LA_Morok: no atom type info on runfile'
            Call WarningMessage(2,Message)
            Call Abend
         End If
         Call mma_allocate(AT,nTot)
         Call Get_iArray('Atom Types',AT,nTot)
         Do iLink = 1,nLink
            isOkLA = AT(DefLA(1,iLink)).Eq.QM
            isOkQM = AT(DefLA(2,iLink)).Eq.QM
            isOkMM = AT(DefLA(3,iLink)).Eq.MMI
            If (.Not.(isOkLa.And.isOkQM.And.isOkMM)) Then
               Message = 'Link atoms badly defined. Check input!'
               Call WarningMessage(2,Message)
               Call Abend()
            End If
         End Do
* Generate vector for translating from Gromacs to Molcas numbering
         Call mma_allocate(GroToMol,nTot)
         iAtIn = 1
         iAtOut = 1
         Do iAt = 1,nTot
            If (AT(iAt).Eq.QM.Or.AT(iAt).Eq.MMI) Then
               GroToMol(iAt) = iAtIn
               iAtIn = iAtIn+1
            Else If (AT(iAt).Eq.MMO) Then
               GroToMol(iAt) = iAtOut
               iAtOut = iAtOut+1
            Else
               Message = 'LA_Morok: unknown atom type'
               Call WarningMessage(2,Message)
               Call Abend()
            End If
         End Do
* Apply Morokuma scheme to gradient...
         If (iMode.Eq.1) Then
            If (iPL.GE.2) Then
               Write(6,*) 'Applying Morokuma scheme to gradient'
            End If
            Do iLink = 1,nLink
               iLA = GroToMol(DefLA(1,iLink))
               iQM = GroToMol(DefLA(2,iLink))
               iMM = GroToMol(DefLA(3,iLink))
               Fact = FactLA(iLink)
               Do ixyz = 1,3
                  CorG(ixyz,iQM) = CorG(ixyz,iQM)
     &                           + CorG(ixyz,iLA)*(1-Fact)
                  CorG(ixyz,iMM) = CorG(ixyz,iMM)
     &                           + CorG(ixyz,iLA)*Fact
                  CorG(ixyz,iLA) = Zero
               End Do
            End Do
* ...or to position
         Else If (iMode.Eq.2) Then
            If (iPL.GE.2) Then
               Write(6,*) 'Applying Morokuma scheme to positions'
            End If
            Do iLink = 1,nLink
               iLA = GroToMol(DefLA(1,iLink))
               iQM = GroToMol(DefLA(2,iLink))
               iMM = GroToMol(DefLA(3,iLink))
               Fact = FactLA(iLink)
               Do ixyz = 1,3
                  CorG(ixyz,iLA) = CorG(ixyz,iQM)
     &                           +(CorG(ixyz,iMM) - CorG(ixyz,iQM))*Fact
               End Do
            End Do
         Else
            Message = 'LA_Morok: wrong iMode'
            Call WarningMessage(2,Message)
            Call Abend()
         End If
         Call mma_deallocate(DefLA)
         Call mma_deallocate(FactLA)
         Call mma_deallocate(AT)
         Call mma_deallocate(GroToMol)
      End If
c
#ifdef _DEBUG_
      Call RecPrt('LA_Morok: coord or grad:',' ',CorG,3,nAtom)
#endif
      Return
      End
