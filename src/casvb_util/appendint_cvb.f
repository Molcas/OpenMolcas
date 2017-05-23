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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine appendint_cvb(c,number,iskip)
      implicit real*8 (a-h,o-z)
      character*(*)c
      character*10 format

      ibegin=len_trim_cvb(c)+1+iskip
      iend=len(c)

      format=' '
      if(number.ge.0)then
        limit=0
        do 100 itens=0,99
        limit=limit+9*(10**itens)
        if(number.le.limit)then
          write(format,'(a,i4.4,a)')'(i',itens+1,')'
          write(c(ibegin:iend),format)number
          return
        endif
100     continue
      else
        mnumber=-number
        limit=0
        do 200 itens=0,99
        limit=limit+9*(10**itens)
        if(mnumber.le.limit)then
          write(format,'(a,i4.4,a)')'(a,i',itens+1,')'
          write(c(ibegin:iend),format)'-',mnumber
          return
        endif
200     continue
      endif
      write(6,*)' Number too large in appendint :',number
      call abend_cvb()
      end
