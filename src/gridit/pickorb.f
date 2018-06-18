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

       Subroutine PickOrb_nosupport(ipNz,ipSort,ipGref,ipSort_ab,
     &  ipGref_ab,ipVol,ipE,ipOcc,ipE_ab,ipOcc_ab,
     &  nShowMOs,nShowMOs_ab,isener,nMOs,myTitle,ipType)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Include 'real.fh'
      Include 'WrkSpc.fh'
#include <SysDef.fh>
      Include 'grid.nosupport.fh'
      character myTitle*(*)


       ispin=index(myTitle,' spin')

       ik_ab=0
       ii_ab=0
       il_ab=0
       nShowMOs_ab=0

       ishift=0
       do i=0, nIrrep-1
        if(nBas(i).gt.0) Then
         do j=1,nBas(i)
          iWork(ipNZ+j-1+ishift)=i+1
          iWork(ipNZ+j-1+ishift+nMOs)=j
         enddo
        endif
        ishift=ishift+nBas(i)
       enddo

       eps=1d-6
c if no input, but TypeIndex contains 123
      if(isAuMO.eq.-1.and.isAll.ne.1) then
       iActive=0
       do i=0,nMOs-1
        if(iWork(ipType+i).ge.3
     &     .and.iWork(ipType+i).le.5) iActive=iActive+1
       enddo
       if(iActive.gt.0) then
        ii=0
        do i=0,nMOs
        if(iWork(ipType+i).ge.3
     &     .and.iWork(ipType+i).le.5) then
          iWork(ipGref+ii)=i+1
          ii=ii+1
         endif
        enddo
        nShowMOs=iActive
       goto 555

       endif
       endif

       if(isAuMo.eq.-1.and.ispin.gt.0) then
               isAuMO=1
               Region(1)=-2+eps
               Region(2)=-eps
               itRange=0
               isEner=0
       endif

        if(isAll.eq.0) then
          if(isAuMO.eq.1.and.itRange.eq.1) then
               s=Region(1)
               Region(1)=MIN(-Region(2),-Region(1))
               Region(2)=MAX(-Region(2),-s)
               s=Region(1)
               Region(1)=-Region(2)
               Region(2)=-s
          endif
          if(isAuMO.eq.-1 .and. itRange.eq.0) then
               isAuMO=1
               Region(1)=-2+eps
               Region(2)=-eps
          endif
          if(isAuMO.eq.-1 .and. isEner.eq.1) then
           Region(1)=-1000
           Region(2)=1000
          endif
        endif
        if(isAll.eq.1) then
         isAuto=1
         Region(1)=-1000
         Region(2)=1000
        endif
c

*
* 1. user defined number of orbitals. No auto function at all.
*
      if(isAuMO.eq.0) then
       ishift=0
       do i=0, nIrrep-1
         if(nBas(i).gt.0) Then
          iWork(ipSort+i)=ishift
* use Sort as temp
         endif
        ishift=ishift+nBas(i)
        enddo
c
      do i=1,nReq
        iia=iReq(i*2-1)
        iib=iReq(i*2)
        if(iia.le.0.or.iia.gt.nIrrep
     &     .or.iib.lt.0.or.iib.gt.nBas(iia-1)) then
        write (6,'(a)') 'Requested orbital does not exist'
        Call Quit_OnUserError()
c
        endif
        iWork(ipGref+i-1)=iWork(ipSort+iia-1)+iib
        if(isUHF.eq.1)
     *        iWork(ipGref_ab+i-1)=iWork(ipGref+i-1)
      enddo

      nShowMOs=nReq
        if(isUHF.eq.1) nShowMOs_ab=nReq
      goto 555
      endif
***************************************************************
* Well. The user didn't make an exact request. we need to choose orbitals.
       if(itRange.eq.0) then
           isEner=0
           R=Region(2)
           Region(2)=-Region(1)
           Region(1)=-R
       endif

       do i=0,nMOs-1
        Work(ipVol+i)=0.0
        iWork(ipSort+i)=0
        if(isEner.eq.0)  then
           Work(ipE+i)=-Work(ipOcc+i)
           if(isUHF.eq.1) then
           Work(ipE_ab+i)=-Work(ipOcc_ab+i)
           endif
        endif
       enddo


*  Well, now we need to choose rest (nGrid-1) grids.
*
*  Make stupid sorting...
*
       if(NoSort.eq.1) then
        ik=0
        do i=0,nMOs-1
          if(Work(ipE+i).gt.Region(1)
     &       .and.Work(ipE+i).lt.Region(2)) then
c       print *,'EE',Work(ipE+i), Region(1),Region(2)
           ik=ik+1
        iWork(ipSort+i)=ik
          endif
        enddo
        ik_ab=ik
       else

        ik=0
        do i=0,nMOs-1
          if(Work(ipE+i).gt.Region(1)
     &       .and.Work(ipE+i).lt.Region(2)) then
            do j=0,nMOs-1
              if(Work(ipE+j).ge.Work(ipE+i).and.
     &                 Work(ipE+j).ge.Region(1).and.
     &           Work(ipE+j).le.Region(2)) then
                 if(Work(ipE+j).eq.Work(ipE+i)) then
                   if(Work(ipOcc+j).le.Work(ipOcc+i)) then
                     iWork(ipSort+i)=iWork(ipSort+i)+1
                       if(ik.lt.iWork(ipSort+i)) ik=iWork(ipSort+i)
                   endif
                 else
                   iWork(ipSort+i)=iWork(ipSort+i)+1
                    if(ik.lt.iWork(ipSort+i)) ik=iWork(ipSort+i)
                 endif
               endif
             enddo
           endif
        enddo

        if(isUHF.eq.1) then
        ik_ab=0
          do i=0,nMOs-1
            if(Work(ipE_ab+i).gt.Region(1).and.
     &         Work(ipE_ab+i).lt.Region(2)) then
              do j=0,nMOs-1
                if(Work(ipE_ab+j).ge.Work(ipE_ab+i).and.
     &                Work(ipE_ab+j).ge.Region(1).and.
     &          Work(ipE_ab+j).le.Region(2)) then
                 if(Work(ipE_ab+j).eq.Work(ipE_ab+i)) then
                   if(Work(ipOcc_ab+j).le.Work(ipOcc_ab+i)) then
                     iWork(ipSort_ab+i)=iWork(ipSort_ab+i)+1
c        print *,'Here', ik_ab,iWork(ipSort_ab+i)
                       if(ik_ab.lt.iWork(ipSort_ab+i))
     &                           ik_ab=iWork(ipSort_ab+i)
                   endif
                 else
                   iWork(ipSort_ab+i)=iWork(ipSort_ab+i)+1
                    if(ik_ab.lt.iWork(ipSort_ab+i))
     &                          ik_ab=iWork(ipSort_ab+i)
                 endif
               endif
             enddo
           endif
        enddo
       endif
       endif
        if(isAuMO.eq.-1 .and. isEner.ne.0.and.isAll.eq.0) then
         ef=-1000.
         ief=0
         ef_ab=ef
         ief_ab=ief
        do i=0,nMOs-1
          if(Work(ipE+i).gt.ef.and.Work(ipOcc+i).gt.eps) then
            ef=Work(ipE+i)
            ief=i
          endif
          if(isUHF.eq.1) then
          if(Work(ipE_ab+i).gt.ef_ab.and.Work(ipOcc_ab+i).gt.eps) then
            ef_ab=Work(ipE_ab+i)
            ief_ab=i
          endif
          endif
        enddo
c       print *,'ef=',ef
c        print *,'ief=',ief, ief_ab
          ii=iWork(ipSort+ief)
          if(isUHF.eq.1) ii_ab=iWork(ipSort_ab+ief_ab)
        do i=0,nMOs-1
          if(iWork(ipSort+i).gt.ii+iMaxUp.or.
     *        iWork(ipSort+i).lt.ii-iMaxDown) iWork(ipSort+i)=0
          if(isUHF.eq.1) then
            if(iWork(ipSort_ab+i).gt.ii_ab+iMaxUp.or.
     *      iWork(ipSort_ab+i).lt.ii_ab-iMaxDown) iWork(ipSort_ab+i)=0
          endif
        enddo
        endif

c666   continue

      il=0
      do j=1,ik
      do i=0,nMOs-1
        if (iWork(ipSort+i).eq.j) then
          iWork(ipGRef+il)= i+1
          il=il+1
        endif
      enddo
      enddo
      if(isUHF.eq.1) then
      il_ab=0
      do j=1,ik_ab
      do i=0,nMOs-1
        if (iWork(ipSort_ab+i).eq.j) then
          iWork(ipGRef_ab+il_ab)= i+1
          il_ab=il_ab+1
        endif
      enddo
      enddo
      endif
      nShowMOs=il
      if(isUHF.eq.1) nShowMOs_ab=il_ab
555   return
      end
