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
      SubRoutine Trnab(Win,Wout,nvec,na,nb)
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
      Real*8 Win(na,nb,nVec), Wout(na,nb,nvec)
*
      Do iVec = 1, nVec
        Call DGeTmo(Win(1,1,ivec),na,na,nb,Wout(1,1,ivec),nb)
      End Do
      Return
      End
