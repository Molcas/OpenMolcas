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
#include "k2.fh"
*
      ipDij = ipOffD(1,ijS)
      mDCRij= ipOffD(2,ijS)
      nDij  = ipOffD(3,ijS)
      ipDSij= ip_Dummy
      If (nr_of_Densities.eq.2) ipDSij=ipOffD(4,ijS)
      If (mDCRij*nDij.ne.0) Then
         ipDDij = ipTmp
         ipTmp = ipTmp + nDij*mDCRij
      Else
         ipDDij = ip_Dummy
      End If
*
      Return
      End
