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
      subroutine prgrad_cvb(grad,n)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension grad(n)

      if(ip(3).lt.2)return
      i1 = mstackr_cvb(norb*norb)
      call mxunfold_cvb(grad,w(i1),norb)
      write(6,'(/,a)')' Orbital gradient :'
      call mxprint_cvb(w(i1),norb,norb,0)
      if(n-nprorb.gt.0)then
        write(6,'(a)')' Structure coefficient gradient :'
        call mxprint_cvb(grad(nprorb+1),1,n-nprorb,0)
      endif
      call mfreer_cvb(i1)
      return
      end
