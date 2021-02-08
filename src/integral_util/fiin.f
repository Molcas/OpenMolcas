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
      Subroutine fiin(lmax)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "welcom.fh"
*
      fiint(0,0)=Pi*Two
      Do 10 i=0,lmax
         Do 11 j=0,lmax-i
            fiint(i,j)=Zero
            Do 12 k=0,j
               a=binom(j,k)
               iexp=i+k
               tal=pi*Two*a*(-One)**k
               if(iexp.ne.0)then
                  Do 13 l=1,iexp
                     al=Two*DBLE(l)
                     tal=tal*(al-One)/al
 13               Continue
                Endif
                fiint(i,j)=fiint(i,j)+tal
 12         Continue
 11      Continue
 10   Continue
*
      Return
      End
