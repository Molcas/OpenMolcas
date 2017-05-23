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
      Subroutine WrCheck(Label,Arr,n)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Arr(n)
      Character*(*) Label
      Write (6,*) Label, DDot_(n,Arr,1,Arr,1),
     &                   DDot_(n,Arr,1,One,0),n
      Return
      End
