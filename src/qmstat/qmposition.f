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
      Subroutine QMPosition(EHam,Cordst,Coord,Forcek
     &                     ,dLJrep,Ract,iQ_Atoms)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"

      Dimension Cordst(MxCen*MxPut,3),Coord(MxAt*3)

*
*-- First the harmonic potential that keeps QM close to centre.
*
      dDepart=(Cordst(1,1)-Coord(1))**2
     &       +(Cordst(1,2)-Coord(2))**2
     &       +(Cordst(1,3)-Coord(3))**2
      EHam=Forcek*0.5d0*dDepart

*
*-- Second the repulsion with boundary that keeps QM away from
*   boundary.
*
      RMax=0.0d0
      Do 901, iAt=1,iQ_Atoms
        R=Cordst(iAt,1)**2+Cordst(iAt,2)**2+Cordst(iAt,3)**2
        R=sqrt(R)
        Diff=Ract-R
        EHam=EHam+(dLJRep/Diff)**12
901   Continue

      Return
      End
