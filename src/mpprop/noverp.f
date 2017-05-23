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
      Subroutine NoverP(n,i,x)
      Implicit Real*8 (a-h,o-z)

      rn=1.0D0
      rp=1.0D0
      If((i.eq.0).or.(i.eq.n)) Then
         x=1
      Else
         Do j=1,i
            rn=rn*(n-j+1)
            rp=rp*j
         EndDo
         x=rn/rp
      EndIf
      Return
      End
