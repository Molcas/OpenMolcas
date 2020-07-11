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
      Subroutine bino(lmax)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "welcom.fh"
*
      Do 2 i=0,10
         Do 3 j=-1,10
            binom(i,j)=Zero
   3     Continue
   2  Continue
      binom(0,0)=One
      if(lmax.eq.0) go to 100
      Do 10 i=1,lmax
         Do 20 j=0,i
            binom(i,j)=binom(i-1,j-1)+binom(i-1,j)
20       Continue
10    Continue
*
100   Continue
      Return
      End
