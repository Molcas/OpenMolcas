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
      subroutine inidf
cbs   initializes the df on common block  with double facultatives
      implicit real*8 (a-h,o-z)
#include "para.fh"
#include "param.fh"
#include "dofuc.fh"
      df(0)=1.d0
      df(1)=1.d0
      do irun=2,ndfmx
      df(irun)=DBLE(irun)*df(irun-2)
      enddo
      do jbm=0,ndfmx-1
      do ibm=jbm,ndfmx
      dffrac(ibm,jbm)=df(ibm)/df(jbm)
      enddo
      enddo
      do jbm=1,ndfmx
      do ibm=0,jbm-1
      dffrac(ibm,jbm)=1d0/dffrac(jbm,ibm)
      enddo
      enddo
      return
      end
