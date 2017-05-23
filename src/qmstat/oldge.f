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
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      Subroutine Oldge(iAcc,Etot,Eold,Ract,Rold)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"

      iAcc=iAcc-1
      Etot=Eold
      Ract=Rold
      icCom=0
      Do 100, i=1,nPart
        Do 101, j=1,nCent
          icCom=icCom+1
          Do 102, k=1,3
            Cordst(icCom,k)=OldGeo(icCom,k)
102       Continue
101     Continue
100   Continue
      Return
      End
