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
      Subroutine ylmnor(lmax)
      Implicit Real*8(A-H,O-Z)
#include "real.fh"
#include "welcom.fh"
*
      Do 10 i=0,lmax
         lm2=i/2
         Do 11 j=0,lm2
            Do 12 k=0,j
               anorm(i,j,k)=fiint(j-k,k)*tetint(i,j)
 12         Continue
 11      Continue
 10   Continue
      Do 20 i=0,lmax
         tal=One/anorm(i,0,0)
         lm2=i/2
         Do 21 j=0,lm2
            Do 22 k=0,j
               anorm(i,j,k)=anorm(i,j,k)*tal
 22         Continue
 21      Continue
 20   Continue
*
      Return
      End
