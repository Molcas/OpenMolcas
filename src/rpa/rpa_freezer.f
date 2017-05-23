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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine RPA_Freezer()
C
C     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
C
C     Figure out symmetry distribution of frozen occupied orbitals.
C
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Integer  RPA_iUHF
      External RPA_iUHF

      Logical Freeze
      Logical Prnt

      Integer iUHF
      Integer iSym
      Integer iSpin
      Integer ip_Fro, l_Fro

      ! set restricted(1)/unrestricted(2)
      iUHF=RPA_iUHF()

      ! freeze orbitals (if requested)
      iSpin=1
      Freeze=nFreeze(iSpin).gt.0
      Do While (.not.Freeze .and. iSpin.lt.iUHF)
         iSpin=iSpin+1
         Freeze=Freeze.or.nFreeze(iSpin).gt.0
      End Do
      If (Freeze) Then
         Prnt=iPrint.ge.0
         l_Fro=nSym
         Call GetMem('OccFrz','Allo','Inte',ip_Fro,l_Fro)
         Do iSpin=1,iUHF
            If (nFreeze(iSpin).gt.0) Then
               Call RPA_Frz(nFreeze(iSpin),Prnt,nSym,
     *                      Work(ip_OccEn(iSpin)),
     *                      nFro(1,iSpin),nOcc(1,iSpin),
     *                      iWork(ip_Fro))
               Do iSym=1,nSym
                  nFro(iSym,iSpin)=nFro(iSym,iSpin)+iWork(ip_Fro-1+iSym)
               End Do
            End If
         End Do
         Call GetMem('OccFrz','Free','Inte',ip_Fro,l_Fro)
      End If

      ! correct number of active occupied orbitals
      Do iSpin=1,iUHF
         Do iSym=1,nSym
            nOcc(iSym,iSpin)=nOcc(iSym,iSpin)-nFro(iSym,iSpin)
         End Do
      End Do

      End
************************************************************************
      Subroutine RPA_Frz(nFre,Prnt,nSym,EOcc,nFro,nOcc,nFro1)
      Implicit None
      Integer nFre
      Logical Prnt
      Integer nSym
      Real*8  EOcc(*)
      Integer nFro(nSym)
      Integer nOcc(nSym)
      Integer nFro1(nSym)
#include "WrkSpc.fh"

      Character*7 SecNam
      Parameter (SecNam='RPA_Frz')

      Integer  Cho_iRange
      External Cho_iRange

      Integer ip_Point, l_Point
      Integer ip_E, l_E
      Integer ip_Pivot, l_Pivot
      Integer ip_iOcc, l_iOcc
      Integer iCount
      Integer NumFre
      Integer iFre
      Integer i
      Integer iSym
      Integer ip1

      Real*8 xMin

      Integer j
      Integer iOcc
      iOcc(j)=iWork(ip_iOcc-1+j)

      If (nSym.lt.1 .or. nSym.gt.8) Then
         Write(6,'(A,I6)') 'nSym=',nSym
         Call RPA_Warn(3,SecNam//': illegal nSym')
      Else If (nSym .eq. 1) Then
         nFro1(1)=max(nFre,0)
         Return
      Else
         Call iZero(nFro1,nSym)
      End If
      If (nFre.lt.1) Return

      l_Point=nFre
      l_iOcc=nSym
      l_E=nOcc(1)
      Do iSym=2,nSym
         l_E=l_E+nOcc(iSym)
      End Do
      l_Pivot=l_E
      If (nFre.gt.l_E) Then
         Call RPA_Warn(4,SecNam//': too many orbitals to freeze')
      End If
      Call GetMem('ScrPnt','Allo','Inte',ip_Point,l_Point)
      Call GetMem('iOcc','Allo','Inte',ip_iOcc,l_iOcc)
      Call GetMem('ScrOccE','Allo','Real',ip_E,l_E)
      Call GetMem('Pivot','Allo','Inte',ip_Pivot,l_Pivot)

      iCount=0
      ip1=ip_iOcc-1
      Do iSym=1,nSym
         iWork(ip1+iSym)=iCount
         iCount=iCount+nOcc(iSym)
      End Do

      iCount=1
      Do iSym=1,nSym
         Call dCopy_(nOcc(iSym),EOcc(iCount+nFro(iSym)),1,
     *                         Work(ip_E+iOcc(iSym)),1)
         iCount=iCount+nFro(iSym)+nOcc(iSym)
      End Do

      xMin=-1.0d15
      NumFre=nFre
      Call dScal_(l_E,-1.0d0,Work(ip_E),1)
      Call CD_DiaMax(Work(ip_E),l_E,iWork(ip_Pivot),iWork(ip_Point),
     *               NumFre,xMin)
      If (NumFre.ne.nFre) Then
         Write(6,'(2(A,I12))') 'NumFre=',NumFre,'  nFre=',nFre
         Call RPA_Warn(3,SecNam//': NumFre != nFre')
      End If

      Do iFre=1,nFre
         iSym=Cho_iRange(iWork(ip_Point-1+iFre),iWork(ip_iOcc),nSym,
     *                   .false.)
         nFro1(iSym)=nFro1(iSym)+1
      End Do

      If (Prnt) Then
         Write(6,'(/,3X,A,A,A)') 'Output from ',SecNam,':'
         Write(6,'(A,I5,A)')
     &   'The',nFre,' lowest occupied orbitals have been frozen.'
         Write(6,'(A)') 'List of frozen occupied orbitals:'
         Do iFre=1,nFre
            iCount=iWork(ip_Point-1+iFre)
            iSym=Cho_iRange(iCount,iWork(ip_iOcc),nSym,.false.)
            i=iCount-iOcc(iSym)
            Write(6,'(1X,A,I5,A,I1,A,F15.8)')
     &      'Occupied orbital',i,' of symmetry ',iSym,
     &      ' and energy ',-Work(ip_E-1+iCount)
         End Do
         Call xFlush(6)
      End If

      Call GetMem('Pivot','Free','Inte',ip_Pivot,l_Pivot)
      Call GetMem('OccE','Free','Real',ip_E,l_E)
      Call GetMem('iOcc','Free','Inte',ip_iOcc,l_iOcc)
      Call GetMem('Point','Free','Inte',ip_Point,l_Point)

      End
