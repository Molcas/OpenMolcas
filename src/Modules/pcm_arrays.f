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
*     PCMSph: coordinates and radii of spheres                         *
*     PCMTess: coordinates and area of tesserae                        *
*     Vert: coordinates of tesserae vertices (3,*)                     *
*     Centr: coordinates of tesserae centers (3,*)                     *
*     SSph: exposed surface of each sphere                             *
*     PCMDM: PCM matrix                                                *
*     PCM_N: atoms where the spheres are centered                      *
*     PCMiSph: sphere to which each tessera belongs                    *
*     NVert: number of vertices of each tessera                        *
*     IntSph: sphere intersecated by each tessera edge (MxVert,*)      *
*     NewSph: parent spheres for each added sphere (2,*)               *
************************************************************************
      Module PCM_arrays
      Integer nTiles
      Real*8, Allocatable:: C_Tessera(:,:), Q_Tessera(:)
      Real*8, Allocatable:: dTes(:,:,:), dPnt(:,:,:,:), dRad(:,:,:),
     &                      dCntr(:,:,:,:)
      Real*8, Allocatable:: PCM_SQ(:,:)   ! PCM solvation charges
      Real*8, Allocatable:: PCMSph(:,:), PCMTess(:,:), Vert(:,:,:),
     &                      Centr(:,:,:), SSph(:), PCMDM(:,:)
      Integer, Allocatable:: PCM_N(:), PCMiSph(:), nVert(:),
     &                       IntSph(:,:), NewSph(:,:)
      End Module PCM_arrays
