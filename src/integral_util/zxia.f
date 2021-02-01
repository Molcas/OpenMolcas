************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1990,2020, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine ZXia(Zeta,ZInv,N,M,Alpha,Beta)
************************************************************************
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             January 1990                                             *
************************************************************************
      Implicit None
#include "real.fh"
      Integer, Intent(In):: N,M
      Real*8, Intent(In):: Alpha(N), Beta(M)
      Real*8, Intent(InOut):: Zeta(N,M), ZInv(N,M)
      Integer j
*
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt(' In ZXia: Alpha',' ',Alpha,N,1)
      Call RecPrt(' In ZXia: Beta ',' ',Beta ,M,1)
#endif
*
      Do j = 1, M
         Zeta(:,j) = (Alpha(:)+Beta(j))
         ZInv(:,j) = One/Zeta(:,j)
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt( ' In ZXia: Zeta',' ',Zeta,N,M)
#endif
*
      Return
      End SubRoutine ZXia
