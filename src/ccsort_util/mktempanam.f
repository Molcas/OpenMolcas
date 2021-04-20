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
       subroutine mktempanam
c
c     this routine prepare names for TEMP and files as
c     TEMP001 - TEMPmbas and store them into
c     tmpnam and tmanam arrays (mbas-maximum number of basis functions)
c
c     variables used:
c     tmp-anam - array of TEMP file names (Transported through common /tmnames/)
c     this routine (I)
c
       implicit real*8 (a-h,o-z)
#include "reorg.fh"
       integer lun,itemp,k1
c
       lun=lunpublic
       call molcas_open(lun,'TEMP000')
c       open (unit=lun,file='TEMP000')
c
       itemp=0
       do 100 k1=1,9
       itemp=itemp+1
       if (itemp.gt.mbas) goto 500
       write (lun,99) k1
 99     format (6hTEMP00,i1)
 100    continue
c
       do 200 k1=10,99
       itemp=itemp+1
       if (itemp.gt.mbas) goto 500
       write (lun,199) k1
 199    format (5hTEMP0,i2)
 200    continue
c
       do 300 k1=100,mbas
       itemp=itemp+1
       if (itemp.gt.mbas) goto 500
       write (lun,299) k1
 299    format (4hTEMP,i3)
 300    continue
c
 500    rewind (lun)
c
       do 600 itemp=1,mbas
       read (lun,599) tmpnam(itemp)
 599    format (a7)
 600    continue
c
       rewind (lun)
       write (lun,*) ' File scratched'
       close (lun)
c
       return
       end
