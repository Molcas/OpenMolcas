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
      Real*8 Function CovRadT(i)
      Implicit Real*8 (a-h,o-z)
#include "covradt_data.fh"
*
      If (i.gt.92) Then
*        Write (*,*) 'CovRadT: i.gt.92'
         CovRadT=1.50d0
      Else
         CovRadT=CovRadT_(i)
      End If
*
      Return
      End
