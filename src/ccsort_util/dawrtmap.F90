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
       subroutine dawrtmap (lun,mapd,mapi,rc)
!
!     this routine write required mapd and mapi to opened unformatted file
!     with number lun
!
!     lun   - Logical unit number of file, where mediate will be stored (Input)
!     mapd  - direct map matrix corresponding to given mediate (Input)
!     mapi  - inverse map matrix corresponding to given mediate (Input)
!     rc    - return (error) code (Output)
!
#include "files_ccsd.fh"
#include "reorg.fh"

#include "SysDef.fh"
!
       integer lun,rc
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
!
       rc=0
!
!1    write mapd
!
       if (iokey.eq.1) then
!      Fortran IO
       write (lun) mapd,mapi
!
       else
!      MOLCAS IO
       call idafile (lun,1,mapd,3078,daddr(lun))
       call idafile (lun,1,mapi,512,daddr(lun))
       end if
!
       return
       end
