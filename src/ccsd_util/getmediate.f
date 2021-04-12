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
       subroutine getmediate (wrk,wrksize,
     & lun,poss0,mapd,mapi,rc)
c
c     this routine read required mediate from opened unformatted file
c     with number lun, and place it statring with the poss0
c     it also reads mapd and mapi of the given mediade, and reconstruct
c     mapd to actual possitions
c
c     lun   - Logical unit number of file, where mediate is stored (Input)
c     poss0 - initial possition in WRK, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Output)
c     mapi  - inverse map matrix corresponding to given mediate (Output)
c     rc    - return (error) code (Output)
c
c     N.B.
c     all mediates are storred as follows
c     1 - mapd, mapi
c     2 - one record with complete mediate
c
#include "wrk.fh"
       integer lun,poss0,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c     help variables
c
       integer length,rc1
c
       rc=0
c
c1    read mapd
c
      call getmap (lun,poss0,length,mapd,mapi,rc1)
c
c2    read mediate in one block
c
       if (length.eq.0) then
c     RC=1 : there is nothing to read, length of mediate is 0
       rc=1
       return
       end if
c
       call rea (lun,length,wrk(poss0))
c
       return
       end
