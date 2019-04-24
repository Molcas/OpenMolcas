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
#ifdef _HAVE_GRID_IT_
*     This should have never been used outside casvb_util
      Subroutine fmove(ia,ib,n)
      Real*8 ia(*),ib(*)
      Integer n
      Call fmove_cvb(ia,ib,n)
      End
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_fmove
      End
#endif
