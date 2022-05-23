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
      subroutine istkpush_cvb(iarr,ival)
      implicit real*8(a-h,o-z)
      dimension iarr(*)

      iarr(2)=iarr(2)+1
      if(iarr(2).gt.iarr(1))then
        write(6,*)' Stack dimension too small :',iarr(1)
        write(6,*)' Tried push of :',ival
        call abend_cvb()
      endif
      iarr(iarr(2))=ival
      return
      end
