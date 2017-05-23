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
      SubRoutine GenBMp_Localisation(D,C,X,nShell,iSym,ColD,ColC,ColX,
     &                               PreFix)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Real*8 D(nShell,nShell), C(nShell,*), X(nShell,*)
      Character*1 ColD, ColC, ColX
      Character*2 PreFix
#include "inflocal.fh"

      Character*12 FilNam

C     Generate bitmap for density.
C     ----------------------------

      Write(FilNam,'(A2,A5,I1,A4)') PreFix,'Dnsty',iSym,'.bmp'
      Call GenBMp_Loc(D,nShell,nShell,FilNam,ColD)

C     Generate bitmap for original MOs.
C     ---------------------------------

      Write(FilNam,'(A2,A5,I1,A4)') PreFix,'MOrig',iSym,'.bmp'
      Call GenBMp_Loc(C,nShell,nOrb2Loc(iSym),FilNam,ColC)

C     Generate bitmap for localised MOs.
C     ----------------------------------

      Write(FilNam,'(A2,A5,I1,A4)') PreFix,'MOloc',iSym,'.bmp'
      Call GenBMp_Loc(X,nShell,nOrb2Loc(iSym),FilNam,ColX)

      End
