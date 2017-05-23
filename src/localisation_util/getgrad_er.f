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
      SubRoutine GetGrad_ER(Functional,GradNorm,R,CMO,nBasis,nOrb2Loc,
     &                      Timing)
C
C     Thomas Bondo Pedersen, November 2005.
C
C     Purpose: compute ER Functional and its gradient norm.
C              The R matrix is computed as R(i,j) = (ij|jj).
C
C     Note: symmetry is NOT allowed (but is not tested!).
C
      Implicit Real*8 (a-h,o-z)
      Real*8  R(nOrb2Loc,nOrb2Loc), CMO(nBasis,nOrb2Loc)
      Logical Timing
#include "WrkSpc.fh"

      Character*10 SecNam
      Parameter (SecNam = 'GetGrad_ER')

      Character*80 Txt

      Parameter (nSym = 1)
      Integer nOcc(nSym), ipCMO(nSym)

C     Initialization.
C     ---------------

      Functional = 0.0d0
      GradNorm = 0.0D0
      If (nOrb2Loc.lt.1 .or. nBasis.lt.1) Return

C     Transpose CMO (only the part to be localised).
C     ----------------------------------------------

      lCMO = nOrb2Loc*nBasis
      Call GetMem('CMO_T','Allo','Real',ipCMO(1),lCMO)
      ip0 = ipCMO(1) - 1
      Do i = 1,nOrb2Loc
         Call dCopy_(nBasis,CMO(1,i),1,Work(ip0+i),nOrb2Loc)
      End Do

C     Compute R.
C     ----------

      nOcc(1) = nOrb2Loc
      irc = -1
      Call Cho_Get_Rij(irc,ipCMO,nOcc,R,Timing)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I6)') 'Cho_Get_Rij returned',irc
         Call SysAbendMsg(SecNam,'Calculation of ER gradient failed:',
     &                    Txt)
      End If
      Call GetMem('CMO_T','Free','Real',ipCMO(1),lCMO)

C     Compute gradient norm and functional.
C     -------------------------------------

      Do i = 1,nOrb2Loc-1
         Functional = Functional + R(i,i)
         Do j = i+1,nOrb2Loc
            GradNorm = GradNorm + (R(i,j)-R(j,i))**2
         End Do
      End Do
      Functional = Functional + R(nOrb2Loc,nOrb2Loc)
      GradNorm = 4.0d0*sqrt(GradNorm)

      End
