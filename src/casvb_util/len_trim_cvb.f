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
      function len_trim_cvb(a)
c  Length of string excluding trailing blanks
      implicit real*8(a-h,o-z)
      character*(*) a

      do 100 i=len(a),1,-1
      if(a(i:i).ne.' ')goto 200
100   continue
      i=0
200   len_trim_cvb=i
      return
      end
