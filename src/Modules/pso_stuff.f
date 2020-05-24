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
     &                      G1(:,:), DS(:), DSVar(:)
      Integer, Allocatable:: ij2K(:)
      Integer nG2, mG2
      Integer nG1, mG1
      Integer mCMO, kCMO
      Integer nDens, nV_K, nZ_p_k, n_Txy, n_ij2K, nSOs1
      Integer ipD0, ipDVar,iD0Lbl,ip_V_K,ip_U_K,ip_Z_p_k,ip_Txy
      End Module PSO_Stuff
