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
      Subroutine Print_Input(Title,nSym,nBas,wSet,nSet)
      Implicit Real*8 (a-h,o-z)

#include "mxdm.fh"
#include "mxave.fh"

      Dimension wSet(nSet),nBas(MxSym)
      Character*72 Title
      Parameter (lPaper=132)

      Call Banner(Title,1,lPaper-7)
      Write(6,*)
      Write(6,*)
      Write(6,*)'Number of symmetries:',nSym
      Write(6,*)'Basis functions:',(nBas(iSym),iSym=1,nSym)
      Write(6,*)'Normalized weights:',(wSet(iS),iS=1,nSet)

      Return
      End
