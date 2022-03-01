!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine ccsort_wrtmediate (wrk,wrksize,                       &
     & lun,mapd,mapi,rc)
!
!     this routine write required mediate to opened unformatted file
!     with number lun
!     it also store mapd and mapi of the given mediade
!
!     lun   - Logical unit number of file, where mediate will be stored (Input)
!     mapd  - direct map matrix corresponding to given mediate (Input)
!     mapi  - inverse map matrix corresponding to given mediate (Input)
!     rc    - return (error) code (Output)
!
!     N.B.
!     all mediates are storred as follows
!     1 - mapd, mapi
!     2 - one record with complete mediate
!
#include "wrk.fh"
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
!
!     help variables
!
       integer im,length,poss0
!
       rc=0
!
!1    write mapd
!
       write (lun) mapd,mapi
!
!2    calculate overall length
!
       length=0
!
       do 100 im=1,mapd(0,5)
       length=length+mapd(im,2)
 100    continue
!
!     write mediate in one block
!
       if (length.eq.0) then
!     RC=1 : there is nothing to write, length of mediate is 0
       rc=1
       return
       end if
!
       poss0=mapd(1,1)
       call ccsort_wri (lun,length,wrk(poss0))
!
       return
       end
