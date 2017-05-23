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
      Subroutine IsItValid(Coo,CooRef,ValidOrNot)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"

      Parameter (dTroskel=1d-4)

      Dimension Coo(MxCen,3),CooRef(MxCen,3)

      Logical ValidOrNot

      ValidOrNot=.true.
*-- Lengths.
      Do 101, i=1,4
        Do 102, j=i+1,5
          dL_test=0.0d0
          dL_ref=0.0d0
          Do 103, k=1,3
            dL_test=dL_test+(Coo(i,k)-Coo(j,k))**2
            dL_ref=dL_ref+(CooRef(i,k)-CooRef(j,k))**2
103       Continue
          If(abs(dL_test-dL_ref).gt.dTroskel) then
            ValidOrNot=.false.
            Go To 999
          Endif
102     Continue
101   Continue

999   Continue

      Return
      End
