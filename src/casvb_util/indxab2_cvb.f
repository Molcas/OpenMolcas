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
      subroutine indxab2_cvb(indxa,indxb,nstra,nstrb,
     >  iocc,nsa,nsb)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension indxa(nsa),indxb(nsb)
      dimension nstra(mxirrep),nstrb(mxirrep)
      dimension iocc(norb+1)

      call izero(nstra,mxirrep)
      call izero(nstrb,mxirrep)
      inda=0
      indb=0
      do 300 iisym=1,mxirrep

      call loopstr0_cvb(iocc,index,nalf,norb)
400   continue
      irp=1
      do 500 ia=1,nalf
      irp=md2h(irp,ityp(iocc(ia)))
500   continue
      if(irp.ne.iisym)goto 600
      inda=inda+1
      nstra(iisym)=nstra(iisym)+1
      indxa(inda)=index
600   call loopstr_cvb(iocc,index,nalf,norb)
      if(index.ne.1)goto 400

      call loopstr0_cvb(iocc,index,nbet,norb)
700   continue
      irp=1
      do 800 ib=1,nbet
      irp=md2h(irp,ityp(iocc(ib)))
800   continue
      if(irp.ne.iisym)goto 900
      indb=indb+1
      nstrb(iisym)=nstrb(iisym)+1
      indxb(indb)=index
900   call loopstr_cvb(iocc,index,nbet,norb)
      if(index.ne.1)goto 700

300   continue
      return
      end
c  *************************
c  ** MOs, read and write **
c  *************************
