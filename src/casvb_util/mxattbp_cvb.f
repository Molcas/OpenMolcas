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
      subroutine mxattbp_cvb(a,b,n1,n2,n3,c)
      implicit real*8 (a-h,o-z)
      dimension a(n2,n1),b(n2,n3),c(n1,n3)

      Call DGEMM_('T','N',n1,n3,n2,1.0d0,a,n2,b,n2,1.0d0,c,n1)

      return
      end
