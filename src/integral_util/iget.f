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
      Integer Function iGet(A,n)
      Use Iso_C_Binding
      Real*8, Target :: A(*)
      Integer, Pointer :: iA(:)
      Call C_F_Pointer(C_Loc(A),iA,[n])
      iGet=iA(n)
      Nullify(iA)
      Return
      End
