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
       subroutine dawrite (lun,irec0,vector,length,recl)
c
c     this routine write vector with required length to
c     opened direct access file lun starting from record number
c     irec0
c
c     lun   - logical unit of direct access file (I)
c     irec0 - initial recored number (I)
c     vector- vector (I)
c     length- number of R8 data to be readed (I)
c     recl  - length of one record in lun in R8 (I)
c
       real*8 vector(1:length)
       integer lun,irec0,length,recl
c
c     help variables
c
       integer ilow,iup,need,irec,i
c
       if (length.eq.0) then
       return
       end if
c
c*    def need,ilow,iup,irec
c
       need=length
       ilow=1
       irec=irec0
       iup=0
c
 1     if (recl.ge.need) then
       iup=iup+need
       else
       iup=iup+recl
       end if
c
       write (lun,rec=irec) (vector(i),i=ilow,iup)
c
       need=need-(iup-ilow+1)
       irec=irec+1
       ilow=ilow+recl
c
       if (need.gt.0) goto 1
c
       return
       end
