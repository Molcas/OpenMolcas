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
      Subroutine WrtLbl(Label1,n1,Label2,n2)
      Character*1 Label1(n1), Label2(n2)
      Do 10 i = 1, n2
         Label2(i) = ' '
 10   Continue
      Do 20 i = 1, Min(n1,n2)
         Label2(i) = Label1(i)
 20   Continue
      Return
      End
