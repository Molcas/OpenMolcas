!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
! MxSph: maximum number of spheres allowed                             *
! MxTs: maximum number of tesserae allowed                             *
! MxVert: maximum number of vertices per tessera                       *
! PCMSph: coordinates and radii of spheres                             *
! PCMTess: coordinates and area of tesserae                            *
! Vert: coordinates of tesserae vertices (3,:)                         *
! Centr: coordinates of tesserae centers (3,:)                         *
! SSph: exposed surface of each sphere                                 *
! PCMDM: PCM matrix                                                    *
! PCM_N: atoms where the spheres are centered                          *
! PCMiSph: sphere to which each tessera belongs                        *
! NVert: number of vertices of each tessera                            *
! IntSph: sphere intersected by each tessera edge (MxVert,:)           *
! NewSph: parent spheres for each added sphere (2,:)                   *
! PCM_SQ: PCM solvation charges                                        *
!***********************************************************************

module PCM_arrays

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MxSph = 1000, MxTs = 6000, MxVert = 20
integer(kind=iwp) :: nTiles
real(kind=wp), parameter :: DiagScale = 1.0694_wp
integer(kind=iwp), allocatable :: IntSph(:,:), NewSph(:,:), nVert(:), PCM_N(:), PCMiSph(:)
real(kind=wp), allocatable :: C_Tessera(:,:), Centr(:,:,:), dCntr(:,:,:,:), dPnt(:,:,:,:), dRad(:,:,:), dTes(:,:,:), MM(:,:), &
                              PCM_SQ(:,:), PCMDM(:,:), PCMSph(:,:), PCMTess(:,:), Q_Tessera(:), SSph(:), Vert(:,:,:)

public :: C_Tessera, Centr, dCntr, DiagScale, dPnt, dRad, dTes, IntSph, MM, MxSph, MxTs, MxVert, NewSph, nTiles, NVert, PCM_N, &
          PCM_SQ, PCMDM, PCMiSph, PCMSph, PCMTess, Q_Tessera, SSph, Vert

end module PCM_arrays
