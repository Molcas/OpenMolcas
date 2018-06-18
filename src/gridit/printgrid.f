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

************************************************************************
      SubRoutine printline_nosupport(unit,line,len)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Include 'real.fh'
      Include 'WrkSpc.fh'
#include <SysDef.fh>
      Include 'grid.nosupport.fh'
      character line*80
      integer unit
      integer len

      if (isBinary .eq. 1) then
        write (unit) line(1:len)
      else
        write (unit,'(A)') line(1:len)
      endif
      Return
      End
***
************************************************************************
      SubRoutine PrintHeader_nosupport(nMOs,nShowMOs,
     &  nShowMOs_ab,nCoor,nInc,
     &  iiCoord,nTypes )
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Include 'real.fh'
      Include 'WrkSpc.fh'
#include <SysDef.fh>
      Include 'grid.nosupport.fh'
      character line*80
      Integer nTypes(7)
      LuVal_=LuVal
      nShowMOs_=nShowMOs
c      call bXML('Header')
c      call iXML('UHF',isUHF)
      do iiUHF=0,isUHF
      if(iiUHF.eq.1) then
      LuVal_=LuVal_ab
      nShowMOs_=nShowMOs_ab
      endif
      if(isLine.eq.0) then
       write (line,'(a,a)') 'VERSION=     ',VERSION
       call printline_nosupport(LuVal_,line,23)
c       write (line,'(a,a)') 'Extension=   ',0
c       call printline_nosupport(LuVal_,line,23)
       write (line,'(a,i10)') 'N_of_MO=     ',nMOs
c       call iXML('nMOs',nMOs)
       call printline_nosupport(LuVal_,line,23)
       write (line,'(a,i10)') 'N_of_Grids=  ',nShowMOs_
c       call iXML('nGrids',nShowMOs_)
       call printline_nosupport(LuVal_,line,23)
       write (line,'(a,i10)') 'N_of_Points= ',nCoor
c       call iXML('nPoints',nCoor)
       call printline_nosupport(LuVal_,line,23)
       write (line,'(a,i10)') 'Block_Size=  ',nInc
c       call iXML('Block Size',nInc)
       call printline_nosupport(LuVal_,line,23)
       ntmp=nCoor/nInc+1
       write (line,'(a,i10)')   'N_Blocks=    ',ntmp
c       call iXML('nBlocks',ntmp)
       call printline_nosupport(LuVal_,line,23)
c new cut off
       write (line,'(a,i10)')   'Is_cutoff=   ',isCutOff
       call printline_nosupport(LuVal_,line,23)
       write (line,'(a,f10.4)') 'CutOff=      ',CutOff
       call printline_nosupport(LuVal_,line,23)
       write (line,'(a,i10)')   'N_P=         ',iiCoord
       call printline_nosupport(LuVal_,line,23)
       write (line,'(a,7I5)')   'N_INDEX=     ',nTypes
       call printline_nosupport(LuVal_,line,48)
      else
       write(line,'(a,2i10)') '# ', nShowMOs_ , nCoor
       call printline_nosupport(LuVal_,line,23)
      endif
      if (isTheOne .eq. 1) goto 777

      write (line,'(a,3i5)') 'Net=         ',iGridNpt(1)-1,
     *        iGridNpt(2)-1,iGridNpt(3)-1
c      call iaXML('Net',iGridNpt,3)
      call printline_nosupport(LuVal_,line,28)

      write (line,'(a,3f12.3)') 'Origin= ',GridOrigin
c      call daXML('Origin',GridOrigin,3)
      call printline_nosupport(LuVal_,line,44)

      write (line,'(a,3f12.3)') 'Axis_1= ',GridAxis1
c      call daXML('Axis 1',GridAxis1,3)
      call printline_nosupport(LuVal_,line,44)
      write (line,'(a,3f12.3)') 'Axis_2= ',GridAxis2
c      call daXML('Axis 2',GridAxis2,3)
      call printline_nosupport(LuVal_,line,44)
      write (line,'(a,3f12.3)') 'Axis_3= ',GridAxis3
c      call daXML('Axis 3',GridAxis3,3)
      call printline_nosupport(LuVal_,line,44)

777   continue
      enddo
c      call eXML('Header')

      Return
      End
************************************************************************

************************************************************************
      SubRoutine PrintCubeHeader(nShowMOs)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Include 'real.fh'
      Include 'WrkSpc.fh'
#include <SysDef.fh>
      Include 'grid.nosupport.fh'
      Character namepack*3
      Character CubeName*80
      Real tmp(3)


      do i=1,nShowMOs
        ii=i+50
        write(namepack,'(i3.3)') i
        ij=index(Project,' ')-1
        CubeName=project(1:ij)//namepack//'.cube'
        call molcas_open(ii,CubeName)
c        open(ii,file=CubeName)
        write(ii,'(a)') 'file in pseudo-Gaussian cube format'
         nGA=nGatoms
        if(isDensity.eq.1.and.i.eq.nShowMOs) then
          write(ii,'(a)') 'Density'
        else
          write(ii,'(a)') 'MO coefficients'
          nGA=-nGA
        endif
c however, coordinates itself should be shifted!
        write(ii,'(I5,3F12.6)') nGA,GridOrigin(1),
     &   GridOrigin(2),GridOrigin(3)
        do i1=1,3
          tmp(i1)=GridAxis1(i1)/(iGridNpt(1)-1)
        enddo
        write(ii,'(I5,3F12.6)') iGridNpt(1),tmp
        do i1=1,3
          tmp(i1)=GridAxis2(i1)/(iGridNpt(2)-1)
        enddo
        write(ii,'(I5,3F12.6)') iGridNpt(2),tmp
        do i1=1,3
          tmp(i1)=GridAxis3(i1)/(iGridNpt(3)-1)
        enddo
        write(ii,'(I5,3F12.6)') iGridNpt(3),tmp
        call molcas_open(37,'coord.gtmp')
c        open(37,file='coord.gtmp')
        do ij=1,nGatoms
          read(37,'(I5,4F12.6)') iAlabel,xtmp,tmp
c          do i1=1,3
c            tmp(i1)=tmp(i1) - GridOrigin(i1)
c          enddo
          write(ii,'(I5,4F12.6)') iAlabel,xtmp,tmp
        enddo
        close(37)
        if(isDensity.ne.1.or.i.ne.nShowMOs) then
          write(ii,'(2i5)') 1,1
        endif
      enddo
      End
************************************************************************

************************************************************************
      SubRoutine PrintCubeMO(xMo,i,nInc)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Include 'real.fh'
      Include 'WrkSpc.fh'
#include <SysDef.fh>
      Include 'grid.nosupport.fh'
      Dimension xMo(*)

      ij=0
      do igau=1,iGridNpt(1)
        do jgau=1,iGridNpt(2)
          kk = 0
          npts3=iGridNpt(3)
          do while (kk.lt.npts3)
            n = 6
            if (npts3-kk.lt.n) n = npts3-kk
            write (i+50,'(6E13.5)') (xMo(k+ij),k=1,n)
            ij=ij+n
            if(ij.ge.nInc) goto 3939

            kk = kk + 6
          end do
        end do
      end do
3939  continue
      End
************************************************************************
      Subroutine OpenGrid(INPORB)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
      Include 'grid.nosupport.fh'
      Character FullName*256
      Character RealName*306
      Character ss*2
      Character Env*40
      Character INPORB*(*)
      logical exist
      logical is_error
      Character Slash
      Character*12 Alpha
      Slash='/'
      if(INPORB(1:1).eq.Slash) then
        write (6,*)
      endif
      LuOrb=isFreeUnit(46)
      iPRGM=0
c      Call Chk_Rdvec(INPORB,LuOrb,isUHF)
      isUHF=0
      LuVal_ab=-99999
      Alpha='.grid'
      if(isUHF.eq.1) then
       Alpha='_a.grid'
      endif
      do iiUHF=0,isUHF
      if(iiUHF.eq.1) then
       Alpha='_b.grid'
      endif
      Env='WorkDir '
      Call getenvf(Env,FullName)
      iPRGM=0
      if(TheName.ne.' ') then
         open(88,file='extra.prgm')
         iPRGM=1
      endif
      if(FullName.eq.' '.or.TheName(1:1).eq.' ') then
         RealName='M2MSI'
      if(isUHF.eq.1) then
       if(iiUHF.eq.0) RealName='AM2MSI'
       if(iiUHF.eq.1) RealName='BM2MSI'
      endif

      else
       l=1
       ii=1
888     i=index(FullName(l:),Slash)
        if(i.gt.0) then
          ii=i+l
          l=l+i+1
          goto 888
        endif
        Call getenvf('Project',Project)
c        Project=FullName(ii:)
       if(TheName.eq.'NEW'.or.TheName.eq.'new'.or.TheName.eq.'New') then
        do i=1,99
         write(ss,'(i2)') i
         if(i.lt.10) ss(1:1)='0'
         RealName=FullName(1:index(FullName,' ')-1)//Slash
     +    //Project(1:index(Project,' ')-1)//'.'
     +    //ss//Alpha
         call f_inquire(RealName,exist)
         if(.not.exist) goto 889
        enddo

       else

         RealName=FullName(1:index(FullName,' ')-1)//Slash
     +    //Project(1:index(Project,' ')-1)//'.'
     +    //TheName(1:index(TheName,' ')-1)//Alpha
       endif
       write(6,*) 'Grid file: ',RealName
      endif

889   continue
      if(TheName.ne.' ') then
c         open(88,file='extra.prgm')
         write(88,'(a,a,a)') ' (file) M2MSI ',
     *        RealName(1:index(RealName,' ')),'  rwsg'
c         close(88)
      endif
      if(iiUHF.eq.0) then
      LuVal=isFreeUnit(49)
      if(isBinary.eq.1) Then
      call molcas_open_ext2(LuVal,RealName,'sequential',
     & 'unformatted',iostat,.false.,irecl,'unknown',is_error)
c        open(unit=LuVal,access='sequential',
c     ,       form='unformatted', file=RealName)
c        write(6,*) '** Create Grid file:',
c     &        RealName(1:index(RealName,' '))
        write(LuVal) 'a'
        if (imoPack .ne. 0) then
          g=2003.9
          i=0
          write (LuVal) g
          write (LuVal) i
        else
          g=1999.0
          write(LuVal) g
        endif
        write(LuVal) Title1
      endif

      if(isBinary.eq.0) Then
        call molcas_open(LuVal,RealName)
c        open(unit=LuVal,file=RealName,Form='FORMATTED')
        if(isLine.eq.1) then
          Write(LuVal,'(a)') '# data in GNUplot format'
          goto 999
        endif
c        write(6,*) '** Create Grid file (in ASCII format):',
c     &        RealName(1:index(RealName,' '))
        if(iGauss.eq.1) then
c         Write(LuVal,'(a)') 'file in Gaussian cube format'

        else
          if (isTheOne.eq.1) then
            write(LuVal,'(a1)') '9'
          else
            if (imoPack .ne. 0) then
              write (LuVal,'(a1)') '1'
              write (LuVal,'(i5)') 0
            else
              write(LuVal,'(a1)') '0'
            endif
          endif
        endif
        if(isDebug.eq.0) then
        Write(Luval,'(a)') Title1
        else
        Write(Luval,'(a,a)') Title1,' DEBUG'
        endif
      endif
      else ! iiUHF
      LuVal_ab=isFreeUnit(51)
      if(isBinary.eq.1) Then
      call molcas_open_ext2(LuVal_ab,RealName,'sequential',
     &  'unformatted',iostat,.false.,irecl,'unknown',is_error)
c        open(unit=LuVal_ab,access='sequential',
c     ,       form='unformatted', file=RealName)
c        write(6,*) '** Create Grid file',
c     &        RealName(1:index(RealName,' '))
        write(LuVal_ab) 'a'
        if (imoPack .ne. 0) then
          g=2003.9
          i=0
          write (LuVal_ab) g
          write (LuVal_ab) i
        else
          g=1999.0
          write(LuVal_ab) g
        endif
        write(LuVal_ab) Title1
      endif

      if(isBinary.eq.0) Then
        call molcas_open(LuVal_ab,RealName)
c        open(unit=LuVal_ab,file=RealName,Form='FORMATTED')
c        write(6,*) '** Create Grid file (in ASCII format):',
c     &        RealName(1:index(RealName,' '))
        if(iGauss.eq.1) then
c         Write(LuVal,'(a)') 'file in Gaussian cube format'

        else
          if (isTheOne.eq.1) then
            write(LuVal_ab,'(a1)') '9'
          else
            if (imoPack .ne. 0) then
              write (LuVal_ab,'(a1)') '1'
              write (LuVal_ab,'(i5)') 0
            else
              write(LuVal_ab,'(a1)') '0'
            endif
          endif
        endif
        if(isDebug.eq.0) then
        Write(Luval_ab,'(a)') Title1
        else
        Write(Luval_ab,'(a,a)') Title1,' DEBUG'
        endif
      endif
      endif
      enddo
999   continue
      if(iPRGM.eq.1) close(88)
      return
      end

************************************************************************
      Subroutine PrintTitles_nosupport(LuVal,nShowMOs,isDensity,nMOs,
     &  iWipGRef, isEner, isVB, WipOcc, iWipType, Crypt,
     &  iWipNZ, WipE, VBocc, ifpartial,isLine,isSphere,
     &  isColor)
      Implicit Real*8 (A-H,O-Z)
c#include "itmax.fh"
c#include "info.fh"
c      Include 'grid.nosupport.fh'
      Character Line*80
      Character Crypt*7, bb
      Dimension iWipGRef(*), WipOcc(*), iWipType(*),
     &          iWipNZ(*), WipE(*)
       Character LineT*10
      iActOrb=0
c      if(nMOs.gt.10) then
c        call OnlyIMayUseIt('V.Veryazov')
c      endif
       LineT='GridName= '
       if(isLine.eq.1) LineT='#GridName='
      do i=1,nShowMOs-isDensity-isSphere-isColor
        j=iWipGRef(i)
        if(isEner.eq.1) then
          if(.not.(isVB.eq.1.and.
     >       WipOcc(j).gt.0d0.and.
     >       WipOcc(j).lt.2d0)) then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            write (line,'(a,i2,i5,f12.4,'' ('',f4.2,'')'',1x,a)')
     *                  LineT,
     *                  iWipNZ(j),iWipNZ(j+nMOs),
     *                  WipE(j),WipOcc(j),bb
            call printline_nosupport(LuVal,line,38)
          else
            iActOrb=iActOrb+1
            write (line,'(2a,i4,5x,a,f4.2,a)') LineT,
     *                  'VB orbital',iActOrb,' (',VBocc,')'
            call printline_nosupport(LuVal,line,38)
          endif
        else
          if(.not.(isVB.eq.1.and.
     >       WipOcc(j).gt.0d0 .and.
     >       WipOcc(j).lt.2d0)) then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            write (line,'(a,i2,i5,'' ('',f8.6,'')'',1x,a)')
     *                  LineT,iWipNZ(j),
     *                  iWipNZ(j+nMOs), WipOcc(j),bb
            call printline_nosupport(LuVal,line,30)
          else
            iActOrb=iActOrb+1
            write (line,'(2a,i4,5x,a,f4.2,a)') LineT,
     *                       'VB orbital',iActOrb,' (',VBocc,')'
            call printline_nosupport(LuVal,line,30)
          endif
        endif
      enddo
      if(isSphere.eq.1) Then
          write (line,'(a,a)') LineT,'  Sphere '
          call printline_nosupport(LuVal,line,19)
      endif
      if(isSphere.eq.1) Then
          write (line,'(a,a)') LineT,'0 Color  '
          call printline_nosupport(LuVal,line,19)
      endif
      if(isDensity.eq.1) Then
        if(ifpartial.eq.0) Then
          write (line,'(a,a)') LineT,'  Density'
          call printline_nosupport(LuVal,line,19)
        else
          write (line,'(a,a)') LineT,'  Density (partial)'
          call printline_nosupport(LuVal,line,29)
        endif
      endif
      return
      end
********************************
      Subroutine DumpM2Msi_nosupport(iRun,Luval,nShowMOs,isDensity,nMOs,
     &  iWipGRef, isVb, WipOcc, WipMO, WipOut, mCoor,
     &  WipVBmat, iGauss, nInc, imoPack, iWipPBlock,
     &  cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,
     &  isTheOne,isLine,isBinary, isEner, iWipType, iWipNZ,WipE,
     &  WLine,nLine,WCoor,iPrintCount,isDebug,
     &  isCutOff, iWipCutOff,isSphere,SphrDist,isColor,SphrColor)
      Implicit Real*8 (A-H,O-Z)
      Dimension xLimits(4)
      Integer iYDelta(3)
      Character*1  cMoBlock(*)
      Character Crypt*7, bb
      Character Line*80
      Dimension iWipGRef(*), WipOcc(*), WipMO(*), WipOut(*),
     &  WipVBmat(*), iWipPBlock(*),iWipType(*), iWipNZ(*),WipE(*),
     &  WLine(nLine,mCoor),WCoor(3,mCoor),iWipCutOff(*),
     &  SphrDist(mCoor),SphrColor(mCoor)
      Include 'WrkSpc.fh'
          iActOrb=0
          iPrintCount=iPrintCount+1
        do i=1, nShowMOs-isDensity-isSphere-isColor
          iMOs=iWipGRef(i)

          if(.not.(isVB.eq.1.and.
     >       WipOcc(iMOs).gt.0d0
     >        .and.
     >       WipOcc(iMOs).lt.2d0)) then
            call outmo_nosupport(iMOs,1,WipMO,dum,WipOut,mCoor,nMOs)
          else
            iActOrb=iActOrb+1
            call outmo_nosupport(0,1,WipMO,WipVBmat((iActOrb-1)*nMOs),
     >           WipOut,mCoor,nMOs)
          endif

          if (iGauss .ne. 0) then
            call PrintCubeMO(WipOut,i,nInc)
c note: BUG???
            goto 3939
          endif
          if(isLine.eq.0) then
            write (line,'(a,i4)') 'Title= ',iMOs
            call printline_nosupport(LuVal,line,12)
          endif
          if(isTheOne.eq.1)then
           if(isLine.eq.0) then
            write(LuVal,'(f18.12)')
     *            (WipOut(j),j=1,mCoor)
           else
             if(i+1.le.nLine) then
               do j=1,mCoor
                 WLine(i+1,j)=WipOut(j)
               enddo
             endif
           endif
            goto 3939
          endif

          if (imoPack .ne. 0) then
c            Call PackBlock(WipOut,iWipPBlock,mCoor,
c     >                     xLimits,iYDelta)
            write (line,9000)
     *                     0,
     *                     (xLimits(j),j=1,4),(iYDelta(j),j=1,3)
            call printline_nosupport(LuVal,line,73)

9000        format ('BHeader=',I2,1X,(4(E10.4,1X),3(I5,1X)))
            if (isBinary .ne. 0) then
c              call IArrToChar(iWipPBlock,cMoBlock,mCoor)
cvv ! NOT CODED YET
              write (LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
            else
              write (LuVal,'(I5)') (iWipPBlock(j), j=0,mCoor-1)
            endif
          else
            if (isBinary .eq. 1) then
c packing late
            if(isCutOff.eq.1) then
              Call GetMem('TMP','Allo','Real',ipCMP,mCoor)
              call dcopy_(mCoor,WipOut,1,Work(ipCMP),1)
              iii=0
              do ii=1,mCoor
                if(iWipCutOff(ii).eq.1) then
                 Work(ipCMP+iii)=WipOut(ii)
                 iii=iii+1
                endif
              enddo

            write (LuVal) (Work(ipCMP+j),j=0,iii-1)
            Call GetMem('TMP','Free','Real',ipCMP,mCoor)
           else
c no cut off
            write (LuVal) (WipOut(j),j=1,mCoor)

           endif
            else
              if(isDebug.eq.0) then
c normal output - just numbers
            if(isCutOff.eq.1) then
              do j=1,mCoor
              if(iWipCutOff(j).eq.1)
     &          write (LuVal, '(E10.4)') WipOut(j)
              enddo
             else
              write (LuVal,'(E10.4)')
     *                      (WipOut(j),j=1,mCoor)
             endif
              else
c extended output -
              write (LuVal,'(E10.4,3f8.4)')
     *       (WipOut(j),WCoor(1,j), WCoor(2,j), Wcoor(3,j) ,j=1,mCoor)
              endif
            endif
          endif
        j=iWipGRef(i)

        if(isEner.eq.1) then
          if(.not.(isVB.eq.1.and.
     >       WipOcc(j).gt.0d0.and.WipOcc(j).lt.2d0))then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            if(iRun.eq.1.and.iPrintCount.eq.1)
     *      write(6,'(a,i2,i5,f12.4,'' ('',f4.2,'') '',a)')
     *            'GridName= ',
     *            iWipNZ(j),iWipNZ(j+nMOs),
     *            WipE(j),WipOcc(j),bb
          else
c            iActOrb=iActOrb+1
            if(iRun.eq.1.and.iPrintCount.eq.1)
     *      write(6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ',
     *            'VB orbital',iActOrb,' (',VBocc,')'
          endif
        else
          if(.not.(isVB.eq.1.and.
     >       WipOcc(j).gt.0d0.and.WipOcc(j).lt.2d0))then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            if(iRun.eq.1.and.iPrintCount.eq.1)
     *      write(6,'(a,i2,i5,'' ('',f8.6,'') '',a)')
     *            'GridName= ',
     *            iWipNZ(j),iWipNZ(j+nMOs),
     *            WipOcc(j),bb
          else
c            iActOrb=iActOrb+1
            if(iRun.eq.1.and.iPrintCount.eq.1)
     *      write(6,'(2a,i4,5x,a,f4.2,a)') 'GridName= ',
     *            'VB orbital',iActOrb,' (',VBocc,')'
          endif
        endif
3939      continue
        enddo

        if(isSphere.eq.1) then
CVV FIXME: no packing.
            write (line,'(a,i4)') 'Title= ',-1
            call printline_nosupport(LuVal,line,12)
          if(isBinary.eq.0) then
           do j=1,mCoor
              write (LuVal,'(E18.12)')
     *        SphrDist(j)
           enddo
          else
              write(LuVal) (SphrDist(j),j=1,mCoor)
          endif
        endif

        if(isColor.eq.1) then
CVV FIXME: no packing.
            write (line,'(a,i4)') 'Title= ',-10
            call printline_nosupport(LuVal,line,12)
          if(isBinary.eq.0) then
           do j=1,mCoor
              write (LuVal,'(E18.12)')
     *        SphrColor(j)
           enddo
          else
              write(LuVal) (SphrColor(j),j=1,mCoor)
          endif
        endif

        if(isDensity.eq.1) then
          call outmo_nosupport(0,2,WipMO,WipOcc,WipOut,
     >               mCoor,nMOs)
          do j=1, mCoor
            dNorm=dNorm+WipOut(j)
          enddo
*****
          old=0.000d0
c          write(6,*) " mCoor=",mCoor
*****
          if(iGauss.ne.0) then
            call PrintCubeMO(WipOut,nShowMOs,nInc)
c BUG??
            goto 3940
          endif
          if(isLine.eq.0) then
            write (line,'(a,i4)') 'Title= ',0
            call printline_nosupport(LuVal,line,12)
          endif
          if (imoPack.ne.0) then
c            Call PackBlock(WipOut,iWipPBlock,mCoor,
c     >                     xLimits,iYDelta)
            write (line,9000)
     *                    0,
     *                    (xLimits(j),j=1,4),(iYDelta(j),j=1,3)
            call printline_nosupport(LuVal,line,73)

            if (isBinary .ne. 0) then
c              call IArrToChar(iWipPBlock,cMoBlock,mCoor)
cvv!!!!!!!
              write (LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
            else
              write (LuVal,'(I5)') (iWipPBlock(j), j=1,mCoor)
            endif
          else
            if (isBinary .eq. 1) then
c packing late
              if(isCutOff.eq.1) then
                Call GetMem('TMP','Allo','Real',ipCMP,mCoor)
                call dcopy_(mCoor,WipOut,1,Work(ipCMP),1)
                iii=0
                do ii=1,mCoor
                  if(iWipCutOff(ii).eq.1) then
                   Work(ipCMP+iii)=WipOut(ii)
                   iii=iii+1
                  endif
                enddo

            write (LuVal) (Work(ipCMP+j),j=0,iii-1)
            Call GetMem('TMP','Free','Real',ipCMP,mCoor)
           else
c no cut off

                 write (LuVal) (WipOut(j),j=1,mCoor)
            endif
            else

               if(isLine.eq.1) then
                do j=1,mCoor
                  WLine(1,j)=WipOut(j)
                enddo
               else
                 if (isDebug.eq.0) then
c normal output - just numbers
            if(isCutOff.eq.1) then
              do j=1,mCoor
                if(iWipCutOff(j).eq.1)
     &            write (LuVal, '(E18.12)') WipOut(j)
              enddo
             else

              write (LuVal,'(E18.12)') (WipOut(j),j=1,mCoor)
              endif
                 else
c extra output -
              write (LuVal,'(E18.12,3f8.4)')
     *       (WipOut(j),WCoor(1,j), WCoor(2,j), Wcoor(3,j) ,j=1,mCoor)
                 endif

               endif
CGG This is only for testing CASDFT functional. It will be restore.
CGG              write (LuVal,'(E10.4)') (Work(j),j=ipOut,ipOut+mCoor-1)
            endif
          endif
3940      continue
        endif
      if(isLine.eq.1) then
       do i=1,mCoor
        write(LuVal,'(3F10.6,22E20.12)')
     *  (WCoor(j,i),j=1,3),(WLine(j,i),j=1,nLine)
       enddo
      endif

*
      Return
      end
********************************
