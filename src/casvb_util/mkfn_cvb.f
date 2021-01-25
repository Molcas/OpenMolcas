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
      subroutine mkfn_cvb(fileid,ibf)
      implicit real*8(a-h,o-z)
#include "io_cvb.fh"
      character*20 fn_tmp
      logical debug
      data debug/.false./

      do 100 i=1,nrec
      if(abs(fileid-fileids(i)).lt.thresh_io)then
        ibf=i
        goto 200
      endif
100   continue
      nrec=nrec+1
      if(nrec.gt.max_rec)then
        write(6,*)' nrec > max_rec in mkfn :',nrec,max_rec
        call abend_cvb()
      endif
      ibf=nrec
c  generate new file name
c  -> must be at most 8 characters to use daname
      fn_tmp=' '
      irec=int(fileid)
      ifile=nint(10*(fileid-irec))
      call appendint_cvb(fn_tmp,irec,0)
      call appendint_cvb(fn_tmp,ifile,0)
      filename(ibf)=fn_tmp(1:len_trim_cvb(fn_tmp))
      fileids(ibf)=fileid
      ifilio(ibf)=0
200   continue
      if(debug)then
        write(6,*)' IO information for identifier :',fileid
        write(6,*)' IBF is :',ibf
        write(6,*)' File name is :',filename(ibf)
      endif
      return
      end
