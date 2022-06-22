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
