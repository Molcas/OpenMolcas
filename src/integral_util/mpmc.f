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
      Module MpmC
*     Coordinate used to compute the multipole moment intergrals.
*     Each order of the multipole moments are associated with an
*     individual origin stored in this arrary.
*     Note that the array (0:n) is mapped onto (1,n+1)
      Real*8, Dimension(:,:), Allocatable :: Coor_MPM
      End Module MpmC

