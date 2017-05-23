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
      subroutine multi_opendir(FNAM, iunit)
      implicit none
c
c Direct fortran I/O with irregular data records
c
c Assume here RECL in byte units  (Same assumption in t3smat.f).
c
c PV/LAOG, 22 may 2003.
c
c
      character FNAM*(*)
      integer iunit, iost
#include "ndisk.fh"
      Logical is_error
*     open(unit=iunit, file=FNAM, access='direct',
*    $     form='unformatted', status='unknown', recl=nblock*8)
      Call MOLCAS_Open_Ext2(iUnit,FNam,
     &                      'direct','unformatted',
     &                      iost,.TRUE.,
     &                      nblock*8,'unknown',is_error)
      If (iost.gt.0 .or. is_error) Then
         Write (6,*) 'Multi_OpenDir: Error opening file!'
      End If
      return
      end



      subroutine multi_readir(G,lg,ifile,ias)
c
c  Direct fortran I/O with irregular data records
c
c  Each data record is assumed to begin aligned with a disk record
c  but may span several disk records.
c
c  The price to pay for this flexibility is
c  1) remembering the correspondance between data record numbers
c     and corresponding initial disk records.
c  2) wasting some disk space, in average about nu*nblock/2
c     where nu is the number of data records and nblock the disk
c     record size.
c
c  Arguments
c     G        Buffer (real*8 words)
c     lg       Buffer length
c     ifile    file unit
c     ias      direct access record to start with
c     (nblock  direct access record length, defined in include file)
c
c  PV/LAOG, 22 may 2003.
c
      implicit none
#include "ndisk.fh"
      integer lg, ifile, ias, iloc,irest,kas,k, last,iopt
      common/ioind/iopt(96)
      real*8 G(lg)
c
      iloc=1
      irest=lg
      kas=ias
c
      do while(irest.gt.0)
      k=min(irest,nblock)
      IF(kas.le.iopt(27))then
      call readir(G(iloc),k,ifile,kas)
      else
      call readir(G(iloc),k,ifile+1,kas-iopt(27))
      endif
      iloc=iloc+k
      irest=irest-k
      kas=kas+1
      enddo
      return
c
c
c
      entry multi_wridir(G,lg,ifile,ias,last)
c
      iloc=1
      irest=lg
      kas=ias

c
      do while(irest.gt.0)
      k=min(irest,nblock)
      IF(kas.le.iopt(27))then
      call wridir(G(iloc),k,ifile,kas)
      else
      call wridir(G(iloc),k,ifile+1,kas-iopt(27))
      endif
      iloc=iloc+k
      irest=irest-k
      kas=kas+1
      enddo
      last=kas-1
      return
      end
