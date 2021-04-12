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
       subroutine ccsort_wrtmap (lun,mapd,mapi,rc)
c
c     this routine write required mapd and mapi to opened unformatted file
c     with number lun
c
c     lun   - Logical unit number of file, where mediate will be stored (Input)
c     mapd  - direct map matrix corresponding to given mediate (Input)
c     mapi  - inverse map matrix corresponding to given mediate (Input)
c     rc    - return (error) code (Output)
c
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
       rc=0
c
c1    write mapd
c
       write (lun) mapd,mapi
c
       return
       end
