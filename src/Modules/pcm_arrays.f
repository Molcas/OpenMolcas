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
      Module PCM_arrays
      Integer nTiles
      Real*8, Allocatable:: C_Tessera(:,:), Q_Tessera(:)
      Real*8, Allocatable:: dTes(:,:,:), dPnt(:,:,:,:), dRad(:,:,:),
     &                      dCntr(:,:,:,:)
      Real*8, Allocatable:: PCM_SQ(:,:)   ! PCM solvation charges
      Real*8, Allocatable:: PCMSph(:,:), PCMTess(:,:), Vert(:,:,:),
     &                      Centr(:,:,:)
      End Module PCM_arrays
