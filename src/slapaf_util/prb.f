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
      Subroutine PrB(dB,nq,nDim)
      Implicit Real*8 (a-h,o-z)
      Real*8 dB(nq,nDim,nDim)
*
      Do iq = 1, nQ
         Write (6,*) ' iq=',iq
         Do iDim = 1, nDim
            Write (6,'(9F10.6)')
     &            (dB(iq,iDim,jDim),jDim=1,nDim)
         End Do
      End Do
*
      Return
      End
