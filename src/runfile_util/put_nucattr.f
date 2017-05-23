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
      Subroutine Put_NucAttr()

      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Character*8 Label
      Integer nSym, nBas(8)
*
      Logical DoEMPC
      Common /EmbPCharg/ DoEMPC
*
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)

      nLT=nBas(1)*(nBas(1)+1)/2
      Do i=2,nSym
         nLT=nLT+nBas(i)*(nBas(i)+1)/2
      End Do
      nLT_=nLT
      If (DoEMPC) nLT=2*nLT
      Call Getmem('tempAtr','Allo','Real',ipAttr,nLT)
      ipXFdInt=ipAttr+nLT_

      irc    = -1
      iOpt   = 6
      iComp  = 1
      iSyLbl = 1
      Label  = 'Attract '
      Call RdOne(irc,iOpt,Label,iComp,Work(ipAttr),iSyLbl)
      If (irc .ne. 0) Then
         Write(6,*) 'Put_NucAttr: RdOne returned ',irc
         Write(6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
         Call SysAbendMsg('Put_NucAttr','I/O error in RdOne',' ')
      End If

      If (DoEMPC) Then
         irc    = -1
         iOpt   = 2
         iComp  = 1
         iSyLbl = 1
         Label  = 'XFdInt  '
         Call RdOne(irc,iOpt,Label,iComp,Work(ipXFdInt),iSyLbl)
         If (irc .ne. 0) Then
            Write(6,*) 'Put_NucAttr: RdOne returned ',irc
            Write(6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
            Call SysAbendMsg('Put_NucAttr','I/O error in RdOne',' ')
         End If
         Call daxpy_(nLT_,1.0d0,Work(ipXFdInt),1,Work(ipAttr),1)
      End If
*
      Call Put_dArray('Nuc Potential',Work(ipAttr),nLT_)
*
#ifdef _DEBUG_
      iAttr=ipAttr
      Do i=1,nSym
         Call TriPrt('Attr Inte','',Work(iAttr),nBas(i))
         iAttr=iAttr+nBas(i)*(nBas(i)+1)/2
      End Do
#endif
*
      Call Getmem('tempAtr','Free','Real',ipAttr,nLT)
*
      Return
      End
