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
      subroutine cvbmn_cvb(icode)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"

#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
      dimension Dummy(1)

c  ICODE=0 standard casvb calculation
c  ICODE=1 variational calculation
c  ICODE=2 end of variational calculation (print summary)

      call cvbstart_cvb_lt9(icode)
      call main_cvb()
      call setretvals_cvb(esym,n_iter)
      call cvbfinish_cvb(icode)
      Call Lucia_Util('CLOSE',iDummy,iDummy,Dummy)
      return
      end
