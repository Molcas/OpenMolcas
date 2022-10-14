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
       subroutine cct3_wri (lun,length,vector)
!
!     this routine write length-R8 numbers to opened unformatted file
!     with number lun at the given possition as one record
!
!     lun    - Logical unit number of file, where mediate will be stored (Input)
!     length - # of R8 numbers to be written  (Input)
!     vector - space, where numbers are stored  (Input)

!
#include "t31.fh"
!
       integer lun,length
       real*8 vector(1:length)
!
       if (iokey.eq.1) then
!      Fortran IO
       write (lun) vector
!
       else
!      MOLCAS IO
       call ddafile (lun,1,vector,length,daddr(lun))
       end if
!
       return
       end
