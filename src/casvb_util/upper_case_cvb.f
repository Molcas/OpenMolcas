************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine upper_case_cvb(line,lenline)
      implicit real*8 (a-h,o-z)
      character*(*) line

      ichar_Aa=ichar('A')-ichar('a')
      do 100 ich=1,lenline
      if(line(ich:ich).ge.'a'.and.line(ich:ich).le.'z')
     >  line(ich:ich)=char(ichar(line(ich:ich))+ichar_Aa)
100   continue
      return
      end
