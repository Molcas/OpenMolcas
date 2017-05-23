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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
      Subroutine CheckQuater(Q)
************************************************************
*
*   <DOC>
*     <Name>CheckQuater</Name>
*     <Syntax>Call CheckQuater(Q)</Syntax>
*     <Arguments>
*       \Argument{Q}{The quaternion to be checked, Dimension(4)}{Real*8}{in}
*     </Arguments>
*     <Purpose>Check whether the quaternion represents a rotation</Purpose>
*     <Dependencies>blas and Call SysAbendMsg</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Check whether the quaternion represents a rotation
*     </Description>
*    </DOC>
*
************************************************************
      Implicit none
#include "WrkSpc.fh"
#include "debug.fh"
#include "real.fh"
      Real*8 Q(0:3)
      Real*8 res
      Real*8 thrs,ddot_
      Real*8 angle
      Real*8 axis(3)
      real*8 modangle
      Parameter(thrs=1d-6)

      res=ddot_(4,Q,1,Q,1)

      if ( abs(res-One).gt.thrs ) then
        Call RecPrt("Quaternion tested",' ',Q,4,1)
        Call SysAbendMsg("CheckQuater",
     &      "Quaternion does not represent a rotation","")
      end if

      angle=modangle(Two * acos(Q(0)),Two*pi)

      if ( debug ) then
        Call RecPrt("Quaternion",' ',Q(0),4,1)
        Write(6,'(a8,f10.6,a3,f10.2,a3)') "Angle = ",
     &       angle,"Rad",180*angle/pi,"Deg"
      end if
      call dcopy_(3,Q(1),1,axis,1)
      Call normalizeVec(axis)
      if ( debug ) then
        Call RecPrt("Axis",' ',axis,3,1)
      end if

      Return
      End
