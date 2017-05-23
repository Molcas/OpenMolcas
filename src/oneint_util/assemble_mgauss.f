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
      Subroutine Assemble_mGauss(As,Ad,nAs)
      Implicit Real*8(a-h,o-z)
      Real*8 As(nAs), Ad(nAs,6)
*
      Call DaXpY_(nAs,1.0D0,Ad(1,1),1,As,1)
      Call DaXpY_(nAs,1.0D0,Ad(1,4),1,As,1)
      Call DaXpY_(nAs,1.0D0,Ad(1,6),1,As,1)
*
      Return
      End
