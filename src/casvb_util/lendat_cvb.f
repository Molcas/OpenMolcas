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
c  *******************
c  ** File handling **
c  *******************
      subroutine lendat_cvb(fileid,len)
c
c        return in len the length of record name in file ifil
c            len = -1 if record does not exist
c
      implicit real*8 (a-h,o-z)
c
      len=-1
      return
c Avoid unused argument warnings
      if (.false.) call Unused_real(fileid)
      end
