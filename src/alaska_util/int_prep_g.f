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
      Subroutine Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)
      Implicit Real*8 (a-h,o-z)
*
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
*
      Integer iSD4(0:nSD,4)
*
      Real*8  Coor(3,4)
      Integer iAOV(4), iStabs(4)
      Logical  Shijij
*
      iCnttp=iSD4(13,1)
      kCnttp=iSD4(13,3)
*
      If (AuxCnttp(iCnttp)) Then
         call dcopy_(3,Work(iSD4(8,2)),1,Coor(1,1),1)
      Else
         call dcopy_(3,Work(iSD4(8,1)),1,Coor(1,1),1)
      End If
      call dcopy_(3,Work(iSD4(8,2)),1,Coor(1,2),1)
*
      If (AuxCnttp(kCnttp)) Then
         call dcopy_(3,Work(iSD4(8,4)),1,Coor(1,3),1)
      Else
         call dcopy_(3,Work(iSD4(8,3)),1,Coor(1,3),1)
      End If
      call dcopy_(3,Work(iSD4(8,4)),1,Coor(1,4),1)
*
      Shijij= iSD4(11,1).eq.iSD4(11,3).and.
     &        iSD4(11,2).eq.iSD4(11,4)
*
      Do iQuad = 1, 4
         iAOV(iQuad)   = iSD4( 7,iQuad)
         iStabs(iQuad) = iSD4(10,iQuad)
      End Do
*
      Return
      End
