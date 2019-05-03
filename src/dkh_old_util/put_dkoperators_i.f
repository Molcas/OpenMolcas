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
      subroutine put_dkoperators_i(i,string,iarray)
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
#include "dkhparameters.fh"
      character*(*) string
      dimension iarray(*)
*
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      koff=(i-1)*nwop
*
      Call ByteCopy(String,iArray(koff+1),MaxLength)
      return
      end
