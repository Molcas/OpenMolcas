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
      Subroutine Size_SO_block_g(iSD4,nSD,Petite,nSO,No_batch)
      Implicit Real*8 (a-h,o-z)
      Integer iSD4(0:nSD,4)
      Logical No_batch, Petite
*
      nSO = MemSO2_P(iSD4( 2,1),iSD4( 2,2),iSD4( 2,3),iSD4( 2,4),
     &               iSD4( 8,1),iSD4( 8,2),iSD4( 8,3),iSD4( 8,4))
      No_batch=nSO.eq.0
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(Petite)
      End
