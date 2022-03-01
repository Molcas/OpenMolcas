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
       subroutine daopen (name,lun,recl,nrec)
!
!     this routine open direct acces file
!
!     name  - name of the file A8 (I)
!     lun   - logical unit number (I)
!     recl - record length in R8 (I)
!     nrec  - number of records (if needed) (I)
!
       integer lun,recl,nrec
       character*8 name
!
!     help variables
!
       integer recln,f_iostat
       logical is_error
!
#ifdef _DECAXP_
       recln=recl*2
#else
       recln=recl*8
#endif
!
       call molcas_open_ext2(lun,name,'direct','unformatted',           &
     &  f_iostat,.true.,recln,'unknown',is_error)
!       open (unit=lun,file=name,form='unformatted',access='direct',
!     & recl=recln)
!
       return
! Avoid unused argument warnings
      if (.false.) call Unused_integer(nrec)
       end
