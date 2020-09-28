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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      Integer Function iDeg(Coor)
************************************************************************
*                                                                      *
* Object: to compute the degeneracy of a coordinate.                   *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             March '91                                                *
************************************************************************
      use Symmetry_Info, only: nIrrep, iOper
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Coor(3), Cx(3,8), r(3)
      Logical New
*
      iDeg = 1
      Cx(:,1) = Coor(:)
      Do i = 1, nIrrep-1
         r(1)=One
         If (iAnd(iOper(i),1).ne.0) r(1)=-One
         r(2)=One
         If (iAnd(iOper(i),2).ne.0) r(2)=-One
         r(3)=One
         If (iAnd(iOper(i),4).ne.0) r(3)=-One
         x=r(1)*Coor(1)
         y=r(2)*Coor(2)
         z=r(3)*Coor(3)
         New=.True.
         Do j = 1, iDeg
            If(New .and. x.eq.Cx(1,j)
     &             .and. y.eq.Cx(2,j)
     &             .and. z.eq.Cx(3,j)) New=.False.
         End Do
         If (New) Then
            iDeg=iDeg+1
            Cx(1,iDeg)=x
            Cx(2,iDeg)=y
            Cx(3,iDeg)=z
         End If
      End Do
*
      Return
      End
