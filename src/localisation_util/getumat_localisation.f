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
      SubRoutine GetUmat_Localisation(U,C,S,X,Scr,lScr,nBas,nOrb)
C
C     Author: T.B. Pedersen
C
C     Purpose: compute transformation matrix U=C^TSX.
C
      Implicit None
      Real*8  U(*), C(*), S(*), X(*)
      Integer lScr
      Real*8  Scr(lScr)
      Integer nBas, nOrb

      Character*80 Txt
      Character*20 SecNam
      Parameter (SecNam = 'GetUmat_Localisation')

      Real*8 d0, d1
      Parameter (d0 = 0.0d0, d1 = 1.0d0)

      Integer Need

      If (nOrb.lt.1 .or. nBas.lt.1) Return

      Need = nBas*nOrb
      If (lScr .lt. Need) Then
         Write(Txt,'(A,I9,A,I9)')
     &   'lScr =',lScr,'     Need =',Need
         Call SysAbendMsg(SecNam,
     &                   'Insufficient dimension of scratch array!',Txt)
      End If

      Call DGEMM_('N','N',nBas,nOrb,nBas,d1,S,nBas,X,nBas,d0,Scr,nBas)
      Call DGEMM_('T','N',nOrb,nOrb,nBas,d1,C,nBas,Scr,nBas,d0,U,nOrb)

      End
