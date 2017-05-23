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
      subroutine setvect(natoms,Geom,XYZ,V1,V2)
************************************************************
*
*   <DOC>
*     <Name>setvect</Name>
*     <Syntax>Call (natoms,Geom,XYZ,V1,V2)</Syntax>
*     <Arguments>
*       \Argument{natoms}{Number of atoms}{Integer}{in}
*       \Argument{Geom}{Geometry to be considered, xyz coordinates, Dimension(3,nat}{Real*8}{in}
*       \Argument{XYZ}{atom index array, Dimension(3)}{Integer}{in}
*       \Argument{V1}{output vector, Dimension(3)}{Real*8}{out}
*       \Argument{V2}{output vector, Dimension(3)}{Real*8}{out}
*     </Arguments>
*     <Purpose>Computes the vectors V1=XYZ(2)-XYZ(1) and V2=XYZ(3)-XYZ(1)
*       using the geometry geom</Purpose>
*     <Dependencies>blas</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*       Computes the vectors V1=XYZ(2)-XYZ(1) and V2=XYZ(3)-XYZ(1)
*       using the geometry geom
*     </Description>
*    </DOC>
*
************************************************************

      implicit none
#include "real.fh"
      Real*8 V1(3),V2(3)
      Real*8 O(3),A(3),B(3)
      Integer natoms
      Real*8 Geom(3,natoms)
      Integer XYZ(3)
c
      call dcopy_(3,Geom(1,XYZ(1)),1,O,1)
      call dcopy_(3,Geom(1,XYZ(2)),1,A,1)
      call dcopy_(3,Geom(1,XYZ(3)),1,B,1)
c
c daxpy(n,da,dx,incx,dy,incy)
      call dcopy_(3,A,1,V1,1)
      call daxpy_(3,-One,O,1,V1,1) ! V1=OA
      call dcopy_(3,B,1,V2,1)
      call daxpy_(3,-One,O,1,V2,1) ! V2=OB
c
      end
