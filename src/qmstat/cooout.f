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
      Subroutine Cooout(Head,Cordst,nPart,nCent)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"

      Dimension Cordst(MxCen*MxPut,3)
      Character Head*200

      Write(6,*)
      Write(6,*)
      Write(6,'(A)')Head
      kaunter=0
      Do 1, i=1,nPart
        Write(6,*)'Molecule ',i
        Do 2, j=1,nCent
          kaunter=kaunter+1
          Write(6,*)(Cordst(kaunter,ii),ii=1,3)
2       Continue
1     Continue

      Return
      End
