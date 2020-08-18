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
      Subroutine Dens_Info(ijS,ipDij,ipDSij,mDCRij,ipDDij,ipTmp,
     &                     nr_of_Densities)
      use k2_arrays
      Implicit Real*8 (a-h,o-z)
*
      ipDij = ipOffD(1,ijS)
      mDCRij= ipOffD(2,ijS)
      nDij  = ipOffD(3,ijS)
      If (nr_of_Densities.eq.2) Then
         ipDSij=ipOffD(4,ijS)
      Else
         ipDSij= 1
      End If
*
      If (mDCRij*nDij.ne.0) Then
         ipDDij = ipTmp
         ipTmp = ipTmp + nDij*mDCRij
      Else
         ipDDij = 1
      End If
*
      Return
      End
