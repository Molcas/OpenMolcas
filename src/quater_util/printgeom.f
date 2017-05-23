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
      subroutine PrintGeom(iLU,Nat,title,Geom,lbl)
************************************************************
*
*   <DOC>
*     <Name>PrintGeom</Name>
*     <Syntax>PrintGeom(Nat,title,Geom,lbl)</Syntax>
*     <Arguments>
*       \Argument{iLU}{Logic unit number}{Integer}{in}
*       \Argument{Nat}{number of atoms}{Integer}{in}
*       \Argument{title}{title to be printed}{Character*20}{in}
*       \Argument{Geom}{XYZ coordinates, Dimension(3,nat)}{Real*8}{in}
*       \Argument{lbl}{Atom labels, Dimension(nat)}{Character*20}{in}
*     </Arguments>
*     <Purpose>Print the geometry</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Print the geometry.
*     </Description>
*    </DOC>
*
************************************************************
      implicit none
#include "debug.fh"
      integer nat,iat,icomp,iLU
      Real*8 Geom(3,nat)
      character*20 title
      character*20 lbl(nat)

      Write(6,'(a8i1)') '--- GEOM'
      Write(6,'(i4)') nat
      Write(6,*) title
      Do iat=1,nat
        if (debug) then
*        lbl is not defined. I'll comment out this statement
*        until some one needs it and replace it with the
*        the standard print out.
*        Write(6,'(i3,5x,a10,8x,3f16.8)') iat,lbl(iat),
*    &             (Geom(icomp,iat),icomp=1,3)
         Write(6,'(i3,5x,3f16.8)') iat,
     &             (Geom(icomp,iat),icomp=1,3)
        else
*        Write(6,   '(5x,a10,8x,3f16.8)')   lbl(iat),
*    &             (Geom(icomp,iat),icomp=1,3)
         Write(6,'(i3,5x,3f16.8)') iat,
     &             (Geom(icomp,iat),icomp=1,3)
        end if
      End Do
      Write(6,'(a8i1)') '--- GEOM'

      if (iLU.ne.-1) then
        Write(iLU,'(i4)') nat
        Write(iLU,*) title
        Do iat=1,nat
*          Write(iLU,   '(a10,8x,3f16.8)')   lbl(iat),
*    &               (Geom(icomp,iat),icomp=1,3)
         Write(6,'(i3,5x,3f16.8)') iat,
     &             (Geom(icomp,iat),icomp=1,3)
        End Do
      End If
      return
c Avoid unused argument warnings
      If (.False.) Call Unused_character(lbl)
      end
