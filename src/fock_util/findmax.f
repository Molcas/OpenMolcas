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

#include "WrkSpc.fh"


       If (transA.eq.'N') Then

          Do jvc=1,numJ

             iLaJ = ipLaJ + numA*(jvc-1)

             XMax = abs(Work(iLaJ))

             Do ia=2,numA

                XMax = Max(XMax,abs(Work(iLaJ+ia-1)))

             End Do

             Work(ipMaxJ+jvc-1) = XMax

          End Do

       ElseIf (transA.eq.'T') Then

          Do jvc=1,numJ

             iLaJ = ipLaJ + jvc - 1

             XMax = abs(Work(iLaJ))

             Do ia=2,numA

                XMax = Max(XMax,abs(Work(iLaJ+numJ*(ia-1))))

             End Do

             Work(ipMaxJ+jvc-1) = XMax

          End Do

       Else

          write(6,*)'FindMax: wrong input argument, transA= ',transA
          Call Abend

       EndIf


       Return
       End
