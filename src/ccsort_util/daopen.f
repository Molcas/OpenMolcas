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
       subroutine daopen (name,lun,recl,nrec)
c
c     this routine open direct acces file
c
c     name  - name of the file A8 (I)
c     lun   - logical unit number (I)
c     recl - record length in R8 (I)
c     nrec  - number of records (if needed) (I)
c
       integer lun,recl,nrec
       character*8 name
c
c     help variables
c
       integer recln,f_iostat
       logical is_error
c
#ifdef _DECAXP_
       recln=recl*2
#else
       recln=recl*8
#endif
c
       call molcas_open_ext2(lun,name,'direct','unformatted',
     &  f_iostat,.true.,recln,'unknown',is_error)
c       open (unit=lun,file=name,form='unformatted',access='direct',
c     & recl=recln)
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nrec)
       end
