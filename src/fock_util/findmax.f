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

       SUBROUTINE FindMax(ipLaJ,transa,numA,numJ,ipMaxJ)

       Implicit Real*8 (a-h,o-z)
       Integer  ipLaJ,numA,numJ,ipMaxJ
       character*1 transA
       Real*8 LaJ, LJa
#include "WrkSpc.fh"

       LaJ(i,j) = Work(ipLaJ-1 + i + NumA*(j-1))
       LJa(i,j) = Work(ipLaJ-1 + i + NumJ*(j-1))

       If (transA.eq.'N') Then

          Do jvc=1,numJ

             XMax = abs(LaJ(1,jvc))

             Do ia=2,numA

                XMax = Max(XMax,abs(LaJ(ia,jvc)))

             End Do

             Work(ipMaxJ-1+jvc) = XMax

          End Do

       ElseIf (transA.eq.'T') Then

          Do jvc=1,numJ

             iLaJ = ipLaJ + jvc - 1

             XMax = abs(LJa(jvc,1))

             Do ia=2,numA

                XMax = Max(XMax,abs(LJa(jvc,ia)))

             End Do

             Work(ipMaxJ-1+jvc) = XMax

          End Do

       Else

          write(6,*)'FindMax: wrong input argument, transA= ',transA
          Call Abend

       EndIf


       Return
       End
