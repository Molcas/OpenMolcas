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
c  CASVB IO routines:
c
c  RDRS --> [RD/WR] [R/I]     [S]         [I]
c
c           Read/   Real/     Increment   Include also
c           write   integer   ioffset     integer offset
c
      subroutine rdr_cvb(vec,n,file_id,ioffset)
      implicit real*8 (a-h,o-z)
      dimension vec(n)

      call rdlow_cvb(vec,n,file_id,ioffset)
      return
      end
