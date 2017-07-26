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
*  setvect
*
*> @brief
*>   Computes the vectors \f$ V_1 = \mathit{XYZ}(2)-\mathit{XYZ}(1) \f$ and
*>   \f$ V_2 = \mathit{XYZ}(3)-\mathit{XYZ}(1) \f$ using the geometry \p geom
*> @author Y. Carissan
*>
*> @details
*> Computes the vectors \f$ V_1 = \mathit{XYZ}(2)-\mathit{XYZ}(1) \f$ and
*> \f$ V_2 = \mathit{XYZ}(3)-\mathit{XYZ}(1) \f$ using the geometry \p geom.
*>
*> @param[in]  natoms Number of atoms
*> @param[in]  Geom   Geometry to be considered, xyz coordinates
*> @param[in]  XYZ    atom index array
*> @param[out] V1     output vector
*> @param[out] V2     output vector
************************************************************************
      subroutine setvect(natoms,Geom,XYZ,V1,V2)
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
