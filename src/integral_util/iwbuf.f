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
      Subroutine iWBuf(Array,nArray)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
      Integer Array(nArray)
*
      Call idWBuf(Array,nArray/RtoI)
*
      Return
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine idWBuf(Array,nArray)
      Use Iso_C_Binding
      Integer :: nArray
      Integer, Target :: Array(nArray)
      Real*8, Pointer :: dArray(:)
      Call C_F_Pointer(C_Loc(Array(1)),dArray,[nArray])
      Call dWBuf(dArray,nArray)
      Nullify(dArray)
      End Subroutine idWBuf
*
      End
