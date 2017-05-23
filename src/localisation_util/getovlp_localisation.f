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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine GetOvlp_Localisation(S,Storage,nBas,nSym)
C
C     Thomas Bondo Pedersen, Dec. 2005.
C
C     Purpose: read the overlap matrix and return in S in lower
C              triangular storage if Storage="Tri" else in full square
C              storage.
C
      Implicit None
      Real*8      S(*)
      Character*3 Storage
      Integer     nSym
      Integer     nBas(nSym)
#include "WrkSpc.fh"

      Character*20 SecNam
      Parameter (SecNam = 'GetOvlp_Localisation')
      Logical Debug
      Parameter (Debug = .False.)

      Character*3 myStorage
      Character*8 Label
      Integer     irc, iOpt, iComp, iSyLbl
      Integer     iSym, l_Scr, ip_Scr, kTri, kSq, l_Tri

      l_Tri = nBas(1)*(nBas(1)+1)/2
      Do iSym = 2,nSym
         l_Tri = l_Tri + nBas(iSym)*(nBas(iSym)+1)/2
      End Do
      l_Scr = l_Tri + 4
      Call GetMem('OvlpScr','Allo','Real',ip_Scr,l_Scr)

      irc    = -1
      iOpt   = 2
      iComp  = 1
      iSyLbl = 1
      Label  = 'Mltpl  0'
      Call RdOne(irc,iOpt,Label,iComp,Work(ip_Scr),iSyLbl)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': RdOne returned ',irc
         Write(6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
         Call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
      End If

      myStorage = Storage
      Call UpCase(myStorage)
      If (myStorage .eq. 'TRI') Then
         Call dCopy_(l_Tri,Work(ip_Scr),1,S,1)
      Else
         kTri = ip_Scr
         kSq  = 1
         Do iSym = 1,nSym
            Call Tri2Rec(Work(kTri),S(kSq),nBas(iSym),Debug)
            kTri = kTri + nBas(iSym)*(nBas(iSym)+1)/2
            kSq  = kSq  + nBas(iSym)**2
         End Do
      End If

      Call GetMem('OvlpScr','Free','Real',ip_Scr,l_Scr)

      End
