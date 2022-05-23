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
      subroutine creatci_cvb(inumber,xident_ci,iaddr,iform,fileid)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      logical debug
      data debug/.false./

      ident_ci=inumber
      xident_ci=dble(ident_ci)
      ident_ci=inumber
      iaddr_ci(ident_ci)=iaddr
      iform_ci(ident_ci)=iform
      fileid_ci(ident_ci)=fileid
      if(debug)then
        write(6,*)' Creating CI vector :',inumber
        write(6,*)' Address            :',iaddr
        write(6,*)' Format             :',iform
        write(6,*)' File identifier    :',fileid
      endif
      return
      end
