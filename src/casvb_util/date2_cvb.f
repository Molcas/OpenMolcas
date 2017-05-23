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
      subroutine date2_cvb(delcpu)
      implicit real*8(a-h,o-z)
      character*120 line

      line=' '
      call datimx(line)
      write(6,'(6a,f10.3,a)')' CASVB completed on ',
     >  line(1:10),line(20:24),
     >  ' at ',line(12:19),' after',delcpu,' CPU seconds'
      return
      entry date1_cvb()
      line=' '
      call datimx(line)
      write(6,'(5a/)')' CASVB started on ',line(1:10),line(20:24),
     >  ' at ',line(12:19)
      return
      end
