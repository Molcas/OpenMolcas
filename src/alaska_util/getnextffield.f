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
      Subroutine GetNextFfield(nextfld,fldname,nOrdOpf,ncmp,
     *                         force,lforce)
      Implicit Real*8 (A-H,O-Z)
      Character*(*) fldname
      Dimension force(lforce)
      nextfld=0
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_character(fldname)
         Call Unused_integer(nOrdOpf)
         Call Unused_integer(ncmp)
         Call Unused_real_array(force)
      End If
      End
