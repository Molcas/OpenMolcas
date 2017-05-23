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
      subroutine wrlow_cvb(vec,n,fileid,ioffset)
      implicit real*8(a-h,o-z)
#include "io_cvb.fh"
#include "idbl_cvb.fh"
      dimension vec(n)
      logical newfile,debug
      data debug/.false./

      if(debug)write(6,*)' wrlow :',n,fileid,ioffset
      call mkfn_cvb(fileid,ibf)
      call ibf2unit_cvb(ibf,lu,newfile)

      if(newfile)call ioopn_cvb(filename(ibf),lu)

      ioffs=ioffset
      call dafupd_cvb(lu,ioffs)
      call dDaFile(lu,1,vec,n,ioffs)
      return
      end
