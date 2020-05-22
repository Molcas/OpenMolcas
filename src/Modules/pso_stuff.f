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
      Integer ipD0, ipDVar,ipDS,ipDSVar,nDens,iD0Lbl,ipG2,nG2,
     &        ipCMO,MCMO,ipG1,nG1,ipCMOa,ipCMOb,ip_V_K,ip_U_K,nV_K,
     &        ip_Z_p_k,nZ_p_k,ip_Txy,n_Txy,ip_Thpkl,ip_ij2K,n_ij2K,
     &        ipDMdiag,nSOs1
      End Module PSO_Stuff
