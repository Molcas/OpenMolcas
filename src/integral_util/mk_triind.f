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
       Subroutine Mk_TriInd()
*
#include "TriInd.fh"
*
       ij = 0
       Do k = 0, I_Max-1
          Do i = 0, k
             j = k - i
             ij = ij + 1
             iTriInd(1,ij)=i
             iTriInd(2,ij)=j
          End Do
       End Do
*
       Return
       End
