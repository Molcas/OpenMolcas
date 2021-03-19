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
      Subroutine gridExpandSelect(SelectStr)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "grid.fh"
      Character SelectStr*120, tmp*120
       ifirst=0
       ilen=0
       tmp=' '
       do i=1,120
         if(ifirst.eq.0.and.SelectStr(i:i).ne.' ') then
           ilen=ilen+1
           tmp(ilen:ilen)=SelectStr(i:i)
           ifirst=1
           goto 10
         endif
         if(ifirst.eq.1.and.SelectStr(i:i).ne.' ') then
           ilen=ilen+1
           tmp(ilen:ilen)=SelectStr(i:i)
         endif
         if(ifirst.eq.1.and.SelectStr(i:i).eq.' ') then
           ilen=ilen+1
           tmp(ilen:ilen)=SelectStr(i:i)
           ifirst=0
         endif
10     continue
       enddo
      if(ilen.lt.2) then
         write(6,*) 'SELEct section is incomplete'
         call Quit_OnUserError()
      endif
!      print *,'current',tmp
      nReq=0
1     istart=1
      iend=index(tmp(istart:),' ')
       ibr=index(tmp(istart:iend),':')
       if(ibr.eq.0) then
         write(6,*) 'Wrong format in SELEct section'
         write(6,*) 'Expecting : sign in >',tmp(istart:iend),'<'
         call Quit_OnUserError()
       endif
!       print *,'v01 >',tmp(istart:istart+ibr-2),'<'
       read(tmp(istart:istart+ibr-2),*,err=20,end=20) isymm
       ibrm=index(tmp(istart+ibr+1:iend),'-')
       if(ibrm.eq.0) then
! the only number
!       print *,'v02 >',tmp(istart+ibr:iend),'<'
         read(tmp(istart+ibr:iend),*,err=20,end=20) iibeg
         iReq(2*nReq+1)=isymm
         iReq(2*nReq+2)=iibeg
         nReq=nReq+1
         if(nReq.gt.MAXGRID) goto 30
       else
!         print *,'v03 >',tmp(istart+ibr:istart+ibr+ibrm-1),'<'
         read(tmp(istart+ibr:istart+ibr+ibrm-1),*,err=20,end=20) iibeg
!        print *,'v04 >',tmp(istart+ibr+ibrm+1:iend),'<'
         read(tmp(istart+ibr+ibrm+1:iend),*,err=20,end=20) iiend
         if(iiend.lt.iibeg) then
          write(6,*) 'Wrong data in SELEct section'
          call Quit_OnUserError()
         endif
        do i=iibeg,iiend
         iReq(2*nReq+1)=isymm
         iReq(2*nReq+2)=i
         nReq=nReq+1
         if(nReq.gt.MAXGRID) goto 30
        enddo
       endif
!      print *,'current',tmp(iend+1:)
      tmp=tmp(iend+1:)
      if(tmp.ne.' ') goto 1
      Return
20    write(6,*) 'Error in analyzing SELECT section'
      call Quit_OnUserError()
30    write(6,*) 'Too many Grids requested'
      call Quit_OnUserError()
      end
