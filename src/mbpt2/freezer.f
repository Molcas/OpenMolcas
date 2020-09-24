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
      SubRoutine Freezer(EAll,nFre,nFro,nFro1,nOcc,nBas,nSym,LocPrt)

      Implicit Real*8 (a-h,o-z)
      Real*8  EAll(*)
      Integer nFro(nSym), nFro1(nSym), nOcc(nSym), nBas(nSym)
      Logical LocPrt
#include "WrkSpc.fh"

      Character*7 SecNam
      Parameter (SecNam = 'Freezer')

      Integer  Cho_iRange
      External Cho_iRange

      Integer iOcc(8)

C     For nSym=1, simply transfer nFre to nFro1.
C     Else initialize nFro1 array.
C     ------------------------------------------

      If (nSym.lt.1 .or. nSym.gt.8) Then
         Write(6,*) SecNam,': illegal nSym = ',nSym
         Call qTrace()
         Call SysAbendMsg(SecNam,'illegal nSym',' ')
      Else If (nSym .eq. 1) Then
         nFro1(1) = nFre
         Return
      Else
         Call Cho_iZero(nFro1,nSym)
      End If

C     Set up array of active occupied orbital energies.
C     -------------------------------------------------

      lPoint = nFre

      iOcc(1) = 0
      lEOcc   = nOcc(1)
      Do iSym = 2,nSym
         iOcc(iSym) = lEOcc
         lEOcc = lEOcc + nOcc(iSym)
      End Do

      Call GetMem('ScrOcc','Allo','Real',ipEOcc,lEOcc)
      Call GetMem('Pivot','Allo','Inte',ipPivot,lEOcc)
      Call GetMem('Point','Allo','Inte',ipPoint,lPoint)

      iCount = 1
      Do iSym = 1,nSym
         kAll = iCount + nFro(iSym)
         kOcc = ipEOcc + iOcc(iSym)
         Call dCopy_(nOcc(iSym),EAll(kAll),1,Work(kOcc),1)
         iCount = iCount + nBas(iSym)
      End Do

C     Find pointers to lowest nFre occupied orbital energies.
C     -------------------------------------------------------

      xMin   = -1.0D15
      NumFre = nFre
      Call dScal_(lEOcc,-1.0D0,Work(ipEOcc),1) ! DiaMax finds MAX values
      Call CD_DiaMax(Work(ipEOcc),lEOcc,iWork(ipPivot),iWork(ipPoint),
     &               NumFre,xMin)
      If (NumFre .ne. nFre) Then
         Write(6,*) SecNam,': an error occurred in CD_DiaMax!'
         Write(6,*) 'NumFre = ',NumFre,' != ',nFre,' = nFre'
         Call qTrace()
         Call SysAbendMsg(SecNam,'CD_DiaMax failure',' ')
      End If

C     Set up nFro1 array.
C     -------------------

      Do iFre = 1,nFre
         iSym = Cho_iRange(iWork(ipPoint-1+iFre),iOcc,nSym,.false.)
         nFro1(iSym) = nFro1(iSym) + 1
      End Do

C     If requested, print.
C     --------------------

      If (LocPrt) Then
         Write(6,'(/,3X,A,A,A)') 'Output from ',SecNam,':'
         Write(6,'(1X,A,I5,A)')
     &   'The',nFre,' lowest occupied orbitals have been frozen.'
         Write(6,'(1X,A)') 'List of frozen occupied orbitals:'
         Do iFre = 1,nFre
            kOcc = iWork(ipPoint-1+iFre)
            jSym = Cho_iRange(kOcc,iOcc,nSym,.false.)
            jOcc = kOcc - iOcc(jSym)
            Write(6,'(1X,A,I5,A,I1,A,F15.8)')
     &      'Occupied orbital',jOcc,' of symmetry ',jSym,
     &      ' and energy ',-Work(ipEOcc-1+kOcc)
         End Do
      End If

C     Free memory.
C     ------------

      Call GetMem('Point','Free','Inte',ipPoint,lPoint)
      Call GetMem('Pivot','Free','Inte',ipPivot,lEOcc)
      Call GetMem('ScrOcc','Free','Real',ipEOcc,lEOcc)

      End
