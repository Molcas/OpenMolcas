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
      Subroutine tetin(lmax)
      Implicit Real*8(A-H,O-Z)
#include "real.fh"
#include "welcom.fh"
*
      Do 100 k=0,lmax
         lm2=k/2
         Do 110 l=0,lm2
            tetint(k,l)=Zero
            m=k-l*2
            Do 120 i=0,l
               tetint(k,l)=tetint(k,l)+binom(l,i)*(-One)**i
     &                    / DBLE(m+i*2+1)
120         Continue
110      Continue
100   Continue
*
      Return
      End
