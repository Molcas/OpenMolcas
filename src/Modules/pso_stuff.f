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
      Module PSO_Stuff
      Logical lPSO,lsa, Case_3C, Case_2C, Case_mp2
      Integer nnP(0:7), iOff_ij2K(8),npos(0:7,3),ipAOrb(0:7,2)
      Real*8, Allocatable:: DMdiag(:,:), Thpkl(:), G2(:,:), CMO(:,:),
     &                      Txy(:,:), V_k(:,:), U_k(:), Z_p_k(:,:),
     &                      G1(:,:), D0(:,:), DVar(:,:), DS(:), DSVar(:)
      Integer, Allocatable:: ij2K(:)
      Integer nG2, mG2
      Integer nG1, mG1
      Integer mCMO, kCMO
      Integer nDens, mDens
      Integer n_Txy, m_Txy
      Integer n_ij2K
      Integer nZ_p_k
      Integer nV_K, nSOs1
      Integer iD0Lbl
      End Module PSO_Stuff
