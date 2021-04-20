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
       subroutine wrtmediate (wrk,wrksize,
     & lun,mapd,mapi,rc)
c
c     this routine write required mediate to opened unformatted file
c     with number lun
c     it also store mapd and mapi of the given mediade
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "wrk.fh"
       integer lun,rc,rc1
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer im,length,poss0
c
       rc=0
c
c1    write mapd
c
      call wrtmap (lun,mapd,mapi,rc1)
c
c2    calculate overall length
c
       length=0
c
       do 100 im=1,mapd(0,5)
       length=length+mapd(im,2)
 100    continue
c
c     write mediate in one block
c
       if (length.eq.0) then
c     RC=1 : there is nothing to write, length of mediate is 0
       rc=1
       return
       end if
c
       poss0=mapd(1,1)
       call wri (lun,length,wrk(poss0))
c
       return
       end
