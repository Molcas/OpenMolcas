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
      Subroutine TransMu(SqMu,nDim,TTot,Temp)
*                                                                      *
************************************************************************
*                                                                      *
*     Transform with TTOT the multipole moment integral matrix
*
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 SqMu(nDim*nDim),TTot(nDim*nDim),Temp(nDim*nDim)
*                                                                      *
************************************************************************
*                                                                      *
      Call DGEMM_('N','N',
     &            nDim,nDim,nDim,
     &            1.0d0,SqMu,nDim,
     &            TTot,nDim,
     &            0.0d0,Temp,nDim)
      Call DGEMM_('T','N',
     &            nDim,nDim,nDim,
     &            1.0d0,Ttot,nDim,
     &            Temp,nDim,
     &            0.0d0,SqMu,nDim)
C     Call RecPrt('Overout',' ',SqMu,nDim,nDim)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
