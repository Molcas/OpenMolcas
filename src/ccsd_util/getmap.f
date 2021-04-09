************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine getmap (lun,poss0,length,mapd,mapi,rc)
c
c     this routine reads mapd and mapi of the given mediade
c     from lun and reconstruct mapd to actual possitions poss0
c
c     lun   - Logical unit number of file, where mediate is stored (Input)
c     poss0 - initial possition in WRK, where mediate will be stored (Input)
c     length- overall length of mediate (Output)
c     mapd  - direct map matrix corresponding to given mediate (Output)
c     mapi  - inverse map matrix corresponding to given mediate (Output)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "filemgr.fh"
#include "ccsd1.fh"

#include "SysDef.fh"
c
       integer lun,poss0,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer poss,im,length
c
       rc=0
c
c1    read mapd
c
       if (iokey.eq.1) then
c      Fortran IO
       read (lun) ((mapd(i,j),i=0,512),j=1,6),
     & (((mapi(l,m,n),l=1,8),m=1,8),n=1,8)
c workaround for a bug in some Intel versions
#ifdef __INTEL_COMPILER
       else if (iokey.lt.0) then
       write(6,*) "this should never happen"
#endif
c
       else
c      MOLCAS IO
       call idafile (lun,2,mapd,513*6,daddr(lun))
       call idafile (lun,2,mapi,8*8*8,daddr(lun))
       end if
c
c2    change possitions in mapd to proper one and calculate overall length
c
       poss=poss0
       length=0
c
       do 100 im=1,mapd(0,5)
c
       mapd(im,1)=poss
       poss=poss+mapd(im,2)
       length=length+mapd(im,2)
c
 100    continue
c
       return
       end
