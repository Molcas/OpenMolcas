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
      subroutine charinsert_cvb(cinsert,linsert,c,lc,ipos,idel)
      implicit real*8 (a-h,o-z)
      character*(*) cinsert,c
      character*300 buff

      buff(1:lc-(ipos+idel)+1)=c(ipos+idel:lc)
      c(ipos:ipos+linsert-1)=cinsert(1:linsert)
      c(ipos+linsert:lc+linsert)=buff(1:lc-(ipos+idel)+1)
      lc=lc+linsert-idel
      return
      end
