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
      subroutine copy_dkoperators(i,iarray1,j,iarray2)
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
#include "dkhparameters.fh"
      parameter(maxscr=(maxlength-1)/4+1)
      dimension iarray1(*),iarray2(*)
      dimension iscr(maxscr)
      lwop=8/intrea()
      nwop=(maxlength-1)/lwop+1
      koff=(i-1)*nwop
      loff=(j-1)*nwop
      do k=1,nwop
        iscr(k)=iarray1(koff+k)
      end do
      do k=1,nwop
       iarray2(loff+k)=iscr(k)
      end do
      Return
      End
