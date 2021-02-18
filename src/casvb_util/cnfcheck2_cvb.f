************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine cnfcheck2_cvb(iconfs,nconf1,nel1,iocc)
      implicit real*8 (a-h,o-z)
      logical locc,lorbs,locc_only,lorbs_only
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


      dimension iconfs(noe,nconf1),iocc(noe)

      if(nconf1.eq.0)then
c  Special case -- iconfs will be ok and adhere to occ no definition.
        nconf1=1
        return
      endif
c  Perform basic checks of configurations :
c  First determine if a consistent definition (orb list or
c  occ numbers) has been used for the configurations.
c
c  We do this in two parses. locc_only/lorbs_only are set in
c  first run through if *any* iconf unamiguously adheres to
c  one of the definitions.
c
c  Second run through if necessary we check again in case iconfs
c  contains *both* types of definitions.
c
      locc_only=.false.
      lorbs_only=.false.
      do 100 iconf=1,nconf1
c  Consistency with occ no definition ?
      locc=.true.
      do 200 i=norb+1,noe
      if(iconfs(i,iconf).ne.0)locc=.false.
200   continue
      nsum=0
      do 300 i=1,norb
      if(iconfs(i,iconf).lt.0.or.iconfs(i,iconf).gt.2)locc=.false.
      nsum=nsum+iconfs(i,iconf)
300   continue
      if(nsum.ne.nel1)locc=.false.
c  Consistency with orb list definition ?
      lorbs=.true.
      do 400 i=nel1+1,noe
      if(iconfs(i,iconf).ne.0)lorbs=.false.
400   continue
      call izero(iocc,norb)
      do 500 i=1,nel1
      if(iconfs(i,iconf).ge.1.and.
     >   iconfs(i,iconf).le.norb)then
        iocc(iconfs(i,iconf))=iocc(iconfs(i,iconf))+1
      else
        lorbs=.false.
      endif
500   continue
      do 600 i=1,norb
      if(iocc(i).gt.2)lorbs=.false.
600   continue
c
      if(locc.and..not.lorbs)then
        locc_only=.true.
      elseif(lorbs.and..not.locc)then
        lorbs_only=.true.
      elseif((.not.lorbs).and..not.locc)then
        write(6,*)' Illegal configuration read ',iconf
        write(6,*)(iconfs(ii,iconf),ii=1,noe)
        call abend_cvb()
      endif
100   continue

      locc=locc_only
      lorbs=lorbs_only
      do 1100 iconf=1,nconf1
      if(locc_only.and.lorbs_only)then
c  Check again ...
c  Consistency with occ no definition ?
        locc=.true.
        do 1200 i=norb+1,noe
        if(iconfs(i,iconf).ne.0)locc=.false.
1200    continue
        nsum=0
        do 1300 i=1,norb
        if(iconfs(i,iconf).lt.0.or.iconfs(i,iconf).gt.2)locc=.false.
        nsum=nsum+iconfs(i,iconf)
1300    continue
        if(nsum.ne.nel1)locc=.false.
c  Consistency with orb list definition ?
        lorbs=.true.
        do 1400 i=nel1+1,noe
        if(iconfs(i,iconf).ne.0)lorbs=.false.
1400    continue
        call izero(iocc,norb)
        do 1500 i=1,nel1
        if(iconfs(i,iconf).ge.1.and.
     >     iconfs(i,iconf).le.norb)then
          iocc(iconfs(i,iconf))=iocc(iconfs(i,iconf))+1
        else
          lorbs=.false.
        endif
1500    continue
        do 1600 i=1,norb
        if(iocc(i).gt.2)lorbs=.false.
1600    continue
      endif
      if(locc.and.lorbs)then
c  Comment out following 5 lines if default should be occ no definition:
        call izero(iocc,norb)
        do 1700 i=1,nel1
        iocc(iconfs(i,iconf))=iocc(iconfs(i,iconf))+1
1700    continue
        call imove_cvb(iocc,iconfs(1,iconf),norb)
        if(noe-norb.gt.0) call izero(iconfs(norb+1,iconf),noe-norb)
      elseif(lorbs)then
        call izero(iocc,norb)
        do 1800 i=1,nel1
        iocc(iconfs(i,iconf))=iocc(iconfs(i,iconf))+1
1800    continue
        call imove_cvb(iocc,iconfs(1,iconf),norb)
        if(noe-norb.gt.0) call izero(iconfs(norb+1,iconf),noe-norb)
      endif
      if(iconf.le.500)then
c Test for repeated configurations :
        do 1900 jconf=1,iconf-1
        do 2000 iorb=1,norb
        if(iconfs(iorb,iconf).ne.iconfs(iorb,jconf))goto 1900
2000    continue
        write(6,'(/,a,2i4)')
     >    ' Fatal error - spatial VB configuration repeated :',
     >    jconf,iconf
        write(6,'(i8,a,20i3)')jconf,'   =>  ',(iocc(ii),ii=1,norb)
        write(6,'(i8,a,20i3)')iconf,'   =>  ',(iocc(ii),ii=1,norb)
        call abend_cvb()
1900    continue
      endif
1100  continue
      return
      end
