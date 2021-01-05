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
      SubRoutine PrintLine(unit,line,len,isBinLuscus)
************************************************************************
* Adapted from SAGIT to work with OpenMolcas (October 2020)            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "grid.fh"
      character line*128
      integer unit
      integer len,ll,li

      if (ISLUSCUS .eq. 1) then
        ll=len
        li=isBinLuscus
c       print *,'before pl ',line,' ',ll,li
        call prt_lusc(unit, line, ll,li)
      else
        if (isBinary .eq. 1) then
          write (unit) line(1:len)
        else
          write (unit,'(A)') line(1:len)
        endif
      end if
      Return
      End
***

************************************************************************
      SubRoutine PrintHeader(nMOs,nShowMOs,nShowMOs_ab,nCoor,nInc,
     &  iiCoord,nTypes,iCRSIZE,NBYTES,NINLINE,nBlocks )
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "grid.fh"
      character line*128
      Integer nTypes(7)
      LuVal_=LuVal
      if(isLUSCUS.eq.1) LuVal_=LID
      nShowMOs_=nShowMOs
c      call bXML('Header')
c      call iXML('UHF',isUHF)
      do iiUHF=0,isUHF
      if(iiUHF.eq.1) then
      LuVal_=LuVal_ab
      if(isLUSCUS.eq.1) LuVal_=LID_ab
      nShowMOs_=nShowMOs_ab
      endif
      IF (ISLUSCUS .EQ. 1) THEN
C debug
C        write(*,*) 'N_of_MO=', nMOs
C        write(*,*) 'N_of_Grids=', nShowMOs_
C        write(*,*) 'N_of_Points=', nCoor
C        write(*,*) 'Block_Size=', nInc
C        write(*,*) 'N_Blocks=', nBlocks
C        write(*,*) 'Is_cutoff=', isCutOff
C        write(*,*) 'CutOff=', Cutoff
C        write(*,*) 'N_P=', iiCoord
C end debug


        IF (nMOs .GT. 99999) THEN
          Write(6,*) 'Number of MO''s can''t be larger that 99999'
          Call Quit_OnUserError()
        END IF
        IF (nShowMOs_ .GT. 9999) THEN
          Write(6,*) 'Number of grids can''t be larger that 9999'
          Call Quit_OnUserError()
        END IF
        IF (nCoor .GT. 99999999) THEN
          Write(6,*) 'Number of points can''t be larger that 99999999'
          Call Quit_OnUserError()
        END IF
        IF (nInc .GT. 999999) THEN
          Write(6,*) 'Block size can''t be larger that 999999'
          Call Quit_OnUserError()
        END IF
        nBlocks=nCoor/nInc+1
        IF (nBlocks .GT. 9999) THEN
          Write(6,*) 'Number of blocks can''t be larger that 9999'
c          write(*,*) 'nBlocks = ', nBlocks
          Call Quit_OnUserError()
        END IF
        IF (isCutOff .GT. 9) THEN
          Write(6,*) 'Wrong cutoff'
          Call Quit_OnUserError()
        END IF
        IF (iiCoord .GT. 99999999) THEN
          Write(6,*) 'N_P can''t be larger that 99999999'
          Call Quit_OnUserError()
        END IF
        WRITE(LINE,'(''<GRID>'')')
        CALL PRINTLINE(LUVAL_, LINE,6,0)

        WRITE(LINE,1000)
     +         nMOs, nShowMOs_, nCoor, nInc, nBlocks,
     +         isCutOff, Cutoff, iiCoord
 1000   FORMAT(1X,'N_of_MO=',I5,1X,'N_of_Grids=',I4,
     +         1X,'N_of_Points=',I8,1X,'Block_Size=',I6,
     +         1X,'N_Blocks=',I4,1X,'Is_cutoff=',I1,
     +         1X,'CutOff=',F8.4,1X,'N_P=',I8)
        CALL PRINTLINE(LUVAL_, LINE, 124,0)
        WRITE(LINE,1010) (nTypes(i),i=1,7)
 1010   FORMAT(1X,'N_INDEX=',7I5)
        CALL PRINTLINE(LUVAL_, LINE, 44,0)
        if(isLUSCUS.eq.0) then
        WRITE(LINE,1020) iGridNpt(1)-1, iGridNpt(2)-1, iGridNpt(3)-1
        else
        WRITE(LINE,1020) iGridNpt(1), iGridNpt(2), iGridNpt(3)
        endif
 1020   FORMAT(1X,'Net=',3I5)
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,1030) GridOrigin
 1030   FORMAT(1X,'Origin=',3(1X,F12.8))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
        WRITE(LINE,1040) GridAxis1
 1040   FORMAT(1X,'Axis_1=',3(1X,F12.3))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
        WRITE(LINE,1050) GridAxis2
 1050   FORMAT(1X,'Axis_2=',3(1X,F12.3))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
        WRITE(LINE,1060) GridAxis3
 1060   FORMAT(1X,'Axis_3=',3(1X,F12.3))
        CALL PRINTLINE(LUVAL_, LINE, 47,0)
c Here we dump all missing information
        if(isLUSCUS.eq.0) then
        WRITE(LINE,'(1x,A,I2)') 'CR_SIZE=',iCRSIZE
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,'(1x,A,I2)') 'PACK=',imoPACK
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,'(1x,A,I3)') 'BYTES=',NBYTES
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        WRITE(LINE,'(1x,A,I3)') 'N_in_Line=',NINLINE
        CALL PRINTLINE(LUVAL_, LINE, 20,0)
        NFIRST=nInc
        if(nBlocks.eq.1) NFIRST=nCoor
        WRITE(LINE,'(1x,A,I8)') 'N_FIRST=',NFIRST
        CALL PRINTLINE(LUVAL_, LINE, 20,0)

        NLAST=nCoor-(nBlocks-1)*nInc
        if(NLAST.eq.0) NLAST=nInc
        WRITE(LINE,'(1x,A,I8)') 'N_LAST=',NLAST
        CALL PRINTLINE(LUVAL_, LINE, 20,0)

        nn=(NFIRST/NINLINE)*(NBYTES*NINLINE+iCRSIZE)
        nn1=NFIRST-(NFIRST/NINLINE)*NINLINE
        if(nn1.gt.0) nn=nn+nn1*NBYTES+iCRSIZE

        WRITE(LINE,'(1x,A,I10)') 'N_OFFSET=',nn
        CALL PRINTLINE(LUVAL_, LINE, 30,0)

        nn=(NLAST/NINLINE)*(NBYTES*NINLINE+iCRSIZE)
        nn1=NLAST-(NLAST/NINLINE)*NINLINE
        if(nn1.gt.0) nn=nn+nn1*NBYTES+iCRSIZE
        WRITE(LINE,'(1x,A,I10)') 'N_LAST_OFFSET=',nn
        CALL PRINTLINE(LUVAL_, LINE, 30,0)
        endif

C
C skip file pointers
C skip end orbital section
C
      ELSE


      if(isLine.eq.0) then
       write (line,'(a,a)') 'VERSION=     ',VERSION
       call PrintLine(LuVal_,line,23,0)
c       write (line,'(a,a)') 'Extension=   ',0
c       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'N_of_MO=     ',nMOs
c       call iXML('nMOs',nMOs)
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'N_of_Grids=  ',nShowMOs_
c       call iXML('nGrids',nShowMOs_)
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'N_of_Points= ',nCoor
c       call iXML('nPoints',nCoor)
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)') 'Block_Size=  ',nInc
c       call iXML('Block Size',nInc)
       call PrintLine(LuVal_,line,23,0)
       nBlocks=nCoor/nInc+1
       write (line,'(a,i10)')   'N_Blocks=    ',nBlocks
c       call iXML('nBlocks',nBlocks)
       call PrintLine(LuVal_,line,23,0)
c new cut off
       write (line,'(a,i10)')   'Is_cutoff=   ',isCutOff
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,f10.4)') 'CutOff=      ',CutOff
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,i10)')   'N_P=         ',iiCoord
       call PrintLine(LuVal_,line,23,0)
       write (line,'(a,7I5)')   'N_INDEX=     ',nTypes
       call PrintLine(LuVal_,line,48,0)
      else
       write(line,'(a,2i10)') '# ', nShowMOs_ , nCoor
       call PrintLine(LuVal_,line,23,0)
      endif
      if (isTheOne .eq. 1) goto 777

      write (line,'(a,3i5)') 'Net=         ',iGridNpt(1)-1,
     *        iGridNpt(2)-1,iGridNpt(3)-1
c      call iaXML('Net',iGridNpt,3)
      call PrintLine(LuVal_,line,28,0)

      write (line,'(a,3f12.3)') 'Origin= ',GridOrigin
c      call daXML('Origin',GridOrigin,3)
      call PrintLine(LuVal_,line,44,0)

      write (line,'(a,3f12.3)') 'Axis_1= ',GridAxis1
c      call daXML('Axis 1',GridAxis1,3)
      call PrintLine(LuVal_,line,44,0)
      write (line,'(a,3f12.3)') 'Axis_2= ',GridAxis2
c      call daXML('Axis 2',GridAxis2,3)
      call PrintLine(LuVal_,line,44,0)
      write (line,'(a,3f12.3)') 'Axis_3= ',GridAxis3
c      call daXML('Axis 3',GridAxis3,3)
      call PrintLine(LuVal_,line,44,0)

777   continue
      END IF
      enddo
c      call eXML('Header')

      Return
      End

************************************************************************
      Subroutine OpenGrid(INPORB)
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "grid.fh"
      Character FullName*256
      Character RealName*306
      Character ss*2
      Character Env*40
      Character INPORB*(*)
      CHARACTER(LEN=512) TMPLUS
      INTEGER RC
c      logical exist
      Logical is_error
      Character Slash
      Character*12 Alpha
      INTEGER LUSOPEN
      EXTERNAL LUSOPEN
      Slash='/'

      LuOrb=isFreeUnit(46)
      iPRGM=0
      Call Chk_Vec_UHF(INPORB,LuOrb,isUHF)
      close(LuOrb)
      LuVal_ab=-99999
      if(isLuscus.eq.1) then
      Alpha='.lus'
      else
      Alpha='.grid'
      endif
      if(isUHF.eq.1) then
      if(isLuscus.eq.1) then
       Alpha='_a.lus'
      else
       Alpha='_a.grid'
      endif
      endif
      do iiUHF=0,isUHF
      if(iiUHF.eq.1) then
      if(isLuscus.eq.1) then
       Alpha='_b.lus'
      else
       Alpha='_b.grid'
      endif
      endif
      Env='WorkDir '
      Call getenvf(Env,FullName)
      iPRGM=0
      if(TheName.ne.' ') then
         call molcas_open(88,'extra.prgm')
         iPRGM=1
      endif
      if(FullName.eq.' '.or.TheName(1:1).eq.' ') then
         RealName='M2MSI'
      if(isUHF.eq.1) then
       if(isLuscus.eq.1) then
       if(iiUHF.eq.0) RealName='AM2L'
       if(iiUHF.eq.1) RealName='BM2L'
       else
       if(iiUHF.eq.0) RealName='AM2MSI'
       if(iiUHF.eq.1) RealName='BM2MSI'
       endif
      endif

      else
       l=1
c       ii=1
c       fullname=outf
888     i=index(FullName(l:),Slash)
        if(i.gt.0) then
c          ii=i+l
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
        enddo

       else

         RealName=FullName(1:index(FullName,' ')-1)//Slash
     +    //Project(1:index(Project,' ')-1)//'.'
     +    //TheName(1:index(TheName,' ')-1)//Alpha
       endif
       write(6,*) 'Grid file: ',RealName(:mylen(RealName))
      endif

      if(TheName.ne.' ') then
c         open(88,file='extra.prgm')
         write(88,'(a,a,a)') ' (file) M2MSI ',
     *        RealName(1:index(RealName,' ')),'  rwsg'
c         close(88)
      endif
      if(iiUHF.eq.0) then
      LuVal=isFreeUnit(49)
      if (ISLUSCUS .EQ. 1) THEN
C FIXME: User can't define luscus input file name
        RC=-1
        if(Thename.eq.' ') then
        if(isUHF.eq.0) TMPLUS(1:)="LUSCUS"
        if(isUHF.eq.1.and.iiUHF.eq.0) TMPLUS(1:8)="alph.lus"
        if(isUHF.eq.1.and.iiUHF.eq.1) TMPLUS(1:8)="beta.lus"
        mm=mylen(TMPLUS)
c        print *,' before 2 lusop', mm
        RC=lusopen(LID,TMPLUS,mm)
        else
        mm=mylen(RealName)
c        print *,' before 2 lusop', mm
        RC=lusopen(LID,RealName,mm)
        endif
C        rc=AixOpn(LID,"LUSCUS",.TRUE.)
        IF (RC .NE. 0) THEN
          write(6,*) 'ERROR: Can''t open luscus file!'
          CALL Abend()
        END IF
      ELSE ! not luscus
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
     +         nMOs, nShowMOs_, nCoor, nInc, nBlocks,
     +         isCutOff, Cutoff, iiCoord
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
          if (isTheOne.eq.1) then
            write(LuVal,'(a1)') '9'
          else
              write(LuVal,'(a1)') '0'
          endif
        if(isDebug.eq.0) then
        Write(Luval,'(a)') Title1
        else
        Write(Luval,'(a,a)') Title1,' DEBUG'
        endif
      endif
      END IF
      else ! iiUHF
      if (ISLUSCUS .EQ. 1) THEN
C FIXME: User can't define luscus input file name
        if(Thename.eq.' ') then
        if(isUHF.eq.1.and.iiUHF.eq.0) TMPLUS="AM2L"
        if(isUHF.eq.1.and.iiUHF.eq.1) TMPLUS="BM2L"
        mm=6
c        print *,' before 1 lusop', mm
        RC=lusopen(LID_ab,TMPLUS,mm)
        else
        mm=mylen(RealName)
c        print *,' before lusop', mm
        RC=lusopen(LID_ab,RealName,mm)
        endif
C        rc=AixOpn(LID,"LUSCUS",.TRUE.)
        IF (RC .NE. 0) THEN
          write(6,*) 'ERROR: Can''t open luscus file!'
          CALL Abend()
        END IF
      endif
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
        if(isDebug.eq.0) then
        Write(Luval_ab,'(a)') Title1
        else
        Write(Luval_ab,'(a,a)') Title1,' DEBUG'
        endif
      endif
      END IF
      enddo
999   continue
      if(iPRGM.eq.1) close(88)
      return
      end

************************************************************************
      Subroutine PrintTitles(LuVal,nShowMOs,isDensity,nMOs,
     &  iWipGRef, isEner,  WipOcc, iWipType, Crypt,
     &  iWipNZ, WipE, VBocc, ifpartial,isLine,isSphere,
     &  isColor, ISLUSCUS, nCoor,nBlocks,nInc)
      Implicit Real*8 (A-H,O-Z)

      Character Line*12000
      Character Crypt*7, bb
      Dimension iWipGRef(*), WipOcc(*), iWipType(*),
     &          iWipNZ(*), WipE(*)
       Character LineT*10
       integer Sizeof8
       Sizeof8=8
      iActOrb=0
       LineT='GridName= '
       if(isLine.eq.1) LineT='#GridName='
c       print *,'here',nInc, nBlocks, nCoor
      if(isLUSCUS.eq.1) then
c        NFIRST=nInc

c        if(nBlocks.eq.1) NFIRST=nCoor

c      ii=0
c      do i=1,nShowMOs
c      write(Line,'(A17,1000i12)') ' File_pointers = ',
c     * ((nFirst*(j-1)+ii)*Sizeof8,j=1,nBlocks)
c      CALL PRINTLINE(LUVAL,LINE,12*nBlocks+17,0)
c      ii=ii+nFirst*nBlocks
c      enddo
      write(Line,'(A,i22)') ' ORBOFF = ', nCoor*nShowMOs*Sizeof8
c     * + 28*nShowMOs - 20
c ?? -20?
      CALL PRINTLINE(LUVAL,LINE,32,0)
      endif
      do i=1,nShowMOs-isDensity-isSphere-isColor
        j=iWipGRef(i)
        if(isEner.eq.1) then
          if(.not.(0.eq.1.and.
     >       WipOcc(j).gt.0d0.and.
     >       WipOcc(j).lt.2d0)) then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1000) iWipNZ(j), iWipNZ(j+nMOs), WipE(j),
     +                         WipOcc(j), bb
              CALL PRINTLINE(LUVAL,LINE,72,0)
            ELSE
              write (line,'(a,i2,i5,f12.4,'' ('',f4.2,'')'',1x,a)')
     *                    LineT,
     *                    iWipNZ(j),iWipNZ(j+nMOs),
     *                    WipE(j),WipOcc(j),bb
              call PrintLine(LuVal,line,38,0)
            END IF
 1000       FORMAT(1X,'GridName= Orbital sym=',i2,' index=',i5,
     +             ' Energ=',F12.4,' occ=', F4.2,' type=',a1)
          else
            iActOrb=iActOrb+1
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1010) 'VB orbital',iActOrb,' (',VBocc,')'
 1010         FORMAT(1X,'GridName= VB_orbital iActOrb= ',I4,
     +               ' occ= ',F4.2)
              CALL PRINTLINE(LUVAL,LINE,45,0)
            ELSE
              write (line,'(2a,i4,5x,a,f4.2,a)') LineT,
     *                    'VB orbital',iActOrb,' (',VBocc,')'
              call PrintLine(LuVal,line,38,0)
            END IF
          endif
        else
          if(.not.(0.eq.1.and.
     >       WipOcc(j).gt.0d0 .and.
     >       WipOcc(j).lt.2d0)) then
            ib=iWipType(j)
            bb=' '
            if(ib.gt.0.and.ib.lt.8) bb=Crypt(ib:ib)
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1020) iWipNZ(j), iWipNZ(j+nMOs),
     +                         WipOcc(j), bb
 1020         FORMAT(1X,'GridName= Orbital sym=',I2,' index=',I5,
     +                  ' occ=',F4.2,' type=',A1)
              CALL PRINTLINE(LUVAL,LINE,53,0)
            ELSE
              write (line,'(a,i2,i5,'' ('',f8.6,'')'',1x,a)')
     *                    LineT,iWipNZ(j),
     *                    iWipNZ(j+nMOs), WipOcc(j),bb
              call PrintLine(LuVal,line,30,0)
            END IF
          else
            iActOrb=iActOrb+1
            IF (ISLUSCUS .EQ. 1) THEN
              WRITE(LINE,1010) iActOrb, VBocc
              CALL PRINTLINE(LUVAL,LINE,45,0)
            ELSE
              write (line,'(2a,i4,5x,a,f4.2,a)') LineT,
     *                         'VB orbital',iActOrb,' (',VBocc,')'
              call PrintLine(LuVal,line,30,0)
            END IF
          endif
        endif
      enddo
      if(isSphere.eq.1) Then
          write (line,'(a,a)') LineT,'  Sphere '
          call PrintLine(LuVal,line,19,0)
      endif
      if(isSphere.eq.1) Then
          write (line,'(a,a)') LineT,'0 Color  '
          call PrintLine(LuVal,line,19,0)
      endif
      if(isDensity.eq.1) Then
        if(ifpartial.eq.0) Then
          IF (ISLUSCUS .EQ. 1) THEN
            LINE=' GridName= Density'
            CALL PRINTLINE(LUVAL,LINE,18,0)
          ELSE
            write (line,'(a,a)') LineT,'  Density'
            call PrintLine(LuVal,line,19,0)
          END IF
        else
          IF (ISLUSCUS .EQ. 1) THEN
            LINE=' GridName= Density (partial)'
            CALL PRINTLINE(LUVAL,LINE,28,0)
          ELSE
            write (line,'(a,a)') LineT,'  Density (partial)'
            call PrintLine(LuVal,line,29,0)
          END IF
        endif
      endif
          IF (ISLUSCUS .EQ. 1) THEN
            LINE=' <DENSITY>'
            CALL PRINTLINE(LUVAL,LINE,10,0)
          endif

      return
* Avoid unused argumet warnings
      if (.false.) then
        call Unused_integer(nBlocks)
        call Unused_integer(nInc)
      end if
      end
********************************
      Subroutine DumpM2Msi(iRun,Luval,LID,nShowMOs,isDensity,nMOs,
     &  iWipGRef,  WipOcc, WipMO, WipOut, mCoor,
     &   iGauss, nInc, imoPack, iWipPBlock,
     &  cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,
     &  isTheOne,isLine,isBinary, isEner, iWipType, iWipNZ,WipE,
     &  WLine,nLine,WCoor,iPrintCount,isDebug,
     &  isCutOff, iWipCutOff,isSphere,SphrDist,isColor,SphrColor,
     +  ISLUSCUS, NBYTES,NINLINE)
      Implicit Real*8 (A-H,O-Z)
#include "SysDef.fh"
      Dimension xLimits(4)
      Integer iYDelta(3)
      Character*1  cMoBlock(*)
      Character Crypt*7, bb
      Character Line*128
c     Character fmt*20
C test
c     character*3 cint
c     character*1 cx(64)
c     data cx /
c    +         '0','1','2','3','4','5','6','7','8','9',
c    +         'a','b','c','d','e','f','g','h','i','j',
c    +         'k','l','m','n','o','p','q','r','s','t',
c    +         'u','v','w','x','y','z','A','B','C','D',
c    +         'E','F','G','H','I','J','K','L','M','N',
c    +         'O','P','Q','R','S','T','U','V','W','X',
c    +         'Y','Z','@','#' /
C test end
      Dimension iWipGRef(*), WipOcc(*), WipMO(*), WipOut(*),
     &   iWipPBlock(*),iWipType(*), iWipNZ(*),WipE(*),
     &  WLine(nLine,mCoor),WCoor(3,mCoor),iWipCutOff(*),
     &  SphrDist(mCoor),SphrColor(mCoor)
       Dimension DumArr(2)
#include "WrkSpc.fh"
c          write(*,*) 'entering DumpM2Msi'
           if(irun.gt.100) print *, iGauss, nbytes, ninc, ninline
          iActOrb=0
          iPrintCount=iPrintCount+1
        do i=1, nShowMOs-isDensity-isSphere-isColor
          iMOs=iWipGRef(i)

          if(.not.(0.eq.1.and.
     >       WipOcc(iMOs).gt.0d0
     >        .and.
     >       WipOcc(iMOs).lt.2d0)) then
            call outmo(iMOs,1,WipMO,dumArr,WipOut,mCoor,nMOs)
c          else
c            iActOrb=iActOrb+1
c            call outmo(0,1,WipMO,WipVBmat(1+(iActOrb-1)*nMOs),
c     >           WipOut,mCoor,nMOs)
          endif

        if(isLine.eq.0.and.isLuscus.eq.0) then
            write (line,'(a,i4)') 'Title= ',iMOs
            call PrintLine(LuVal,line,12,1)
          endif
          if(isTheOne.eq.1)then
            if(isLine.eq.0 .AND. ISLUSCUS .eq. 0) then
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
            call PrintLine(LuVal,line,73,0)

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
                Call GetMem('TMP','ALLO','REAL',ipCMP,mCoor)
                call dcopy_(mCoor,WipOut,1,Work(ipCMP),1)
                iii=0
                do ii=1,mCoor
                  if(iWipCutOff(ii).eq.1) then
                   Work(ipCMP+iii)=WipOut(ii)
                   iii=iii+1
                  endif
                enddo

                write (LuVal) (Work(ipCMP+j),j=0,iii-1)
                Call GetMem('TMP','FREE','REAL',ipCMP,mCoor)
              else
c no cut off
                write (LuVal) (WipOut(j),j=1,mCoor)

              endif
            else !isBinary
c             write(cint, '(i3.3)') i
              if(isDebug.eq.0) then
c normal output - just numbers
               if(isCutOff.eq.1) then
                  do j=1,mCoor
                    if(iWipCutOff(j).eq.1)
     &                write (LuVal, '(E10.4)') WipOut(j)
                  enddo
               else if (isLUSCUS .eq. 1) then
c NOPACKING
                if(imoPack.eq.0) then
                 call dump_lusc(LID, WipOut,mCoor)
c debug dump of data
c
c                  do j=1,mCoor,NINLINE
c                    num=mCoor-j
c                    if(num.gt.NINLINE) num=NINLINE
c                    write(fmt,'(A,I2,A,I2,A,A)')
c     &                       '(',NINLINE,'E',NBYTES,'.4',')'
c                    write(line, fmt) (WipOut(j-1+ij),ij=1,num)
c                    call printline(LID,line,NINLINE*NBYTES,0)
c                  enddo
                else
c PACKING NOT implemented
C                  open(unit=38, file='testgr_'//cint//'.txt')
                   do j = 1, min(mcoor,10)
                     iexpnt=int(log10(abs(wipout(j))))
                     write(6,*) 'num=', wipout(j), 'exp = ', iexpnt
                     dnum=(1.0D0 + wipout(j)/10.0D0**iexpnt)/2.0D0

                     in1 = int(dnum*64.0d0)
                     in2 = int((dnum-dble(in1)*64.0D0)*4096.0D0)
                     iexpnt = iexpnt + 50
                     write(6,'(1x,3(1x,e18.8),2x,3(1x,i3))')
     +                      wipout(j), dnum, dexpnt,
     +                      in1, in2, iexpnt
                     if (iexpnt .lt. 1) then
                       iexpnt=1
                     else if (iexpnt .gt. 64) then
                       iexpnt=64
                     end if
                     write(6,'(1x,3(1x,e18.8),2x,3(1x,i3))')
     +                      wipout(j), dnum, dexpnt,
     +                      in1, in2, iexpnt
c                     write(*,*) '-----------------------'

C     +                      cx(in1), cx(in2), cx(iexpnt)
                   end do
C                  close(38)
c                  xxxmin=9.99d+99
c                  xxxmax=-9.99d+99
c                  do j=1,mCoor
c                    if(wipout(j) .gt. xxxmax) xxxmax = wipout(j)
c                    if(wipout(j) .lt. xxxmin) xxxmin = wipout(j)
c                  end do
C                  write(*, *) 'test min/max', xxxmin, xxxmax


                endif !imoPack
               else!isCutOff
C writing of data
                  write (LuVal,'(E10.4)')
     *                      (WipOut(j),j=1,mCoor)
               endif !isCutOff, isLuscus
              else !isDebug
c extended output -
                write (LuVal,'(E10.4,3f8.4)')
     *       (WipOut(j),WCoor(1,j), WCoor(2,j), Wcoor(3,j) ,j=1,mCoor)
              endif !isDebug
            endif !isBinary
          endif !imoPack

        j=iWipGRef(i)

        if(isEner.eq.1) then
          if(.not.(0.eq.1.and.
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
          if(.not.(0.eq.1.and.
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
            call PrintLine(LuVal,line,12,1)
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
            call PrintLine(LuVal,line,12,1)
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
          call outmo(0,2,WipMO,WipOcc,WipOut,
     >               mCoor,nMOs)
          do j=1, mCoor
            dNorm=dNorm+WipOut(j)
          enddo
*****
c          write(6,*) " mCoor=",mCoor
*****
          if(isLine.eq.0.and.IsLuscus.eq.0) then
            write (line,'(a,i4)') 'Title= ',0
            call PrintLine(LuVal,line,12,1)
          endif
          if (imoPack.ne.0) then
c         print *,'pack code'
c            Call PackBlock(WipOut,iWipPBlock,mCoor,
c     >                     xLimits,iYDelta)
c            write (line,9000)
c     *                    0,
c     *                    (xLimits(j),j=1,4),(iYDelta(j),j=1,3)
c            call PrintLine(LuVal,line,73,0)
c
c            if (isBinary .ne. 0) then
c              call IArrToChar(iWipPBlock,cMoBlock,mCoor)
cvv!!!!!!!
c              IF (ISLUSCUS .EQ. 1) THEN
c                RC=C_WRITE(LID, CMOBLOCK, (mCoor*nBytesPackedVal)*RTOB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
c                IF (RC .EQ. 0) THEN
c                  WRITE(6,*) 'error in writing luscus file!'
c                  CALL Abend()
c                END IF
c              ELSE
c                write (LuVal) (cMoBlock(j),j=1,mCoor*nBytesPackedVal)
c              END IF
c            else
c              IF (ISLUSCUS .EQ. 1) THEN
c                RC=C_WRITE(LID, IWIPPBLOCK, (mCoor)*RTOB) !!!!!!!!!!!!!!!!!!!!check mCoor*nBytesPackedVal
c                IF (RC .EQ. 0) THEN
c                  WRITE(6,*) 'error in writing luscus file!'
c                  CALL Abend()
c                END IF
c              ELSE
c                write (LuVal,'(I5)') (iWipPBlock(j), j=1,mCoor)
c              END IF
c            endif
          else
              IF (ISLUSCUS .EQ. 1) THEN
                 call dump_lusc(LID, WipOut,mCoor)
               Endif

            if (isBinary .eq. 1) then
c packing late
              if(isCutOff.eq.1) then
                Call GetMem('TMP','ALLO','REAL',ipCMP,mCoor)
                call dcopy_(mCoor,WipOut,1,Work(ipCMP),1)
                iii=0
                do ii=1,mCoor
                  if(iWipCutOff(ii).eq.1) then
                   Work(ipCMP+iii)=WipOut(ii)
                   iii=iii+1
                  endif
                enddo

            IF (ISLUSCUS .EQ. 1) THEN
              !!!!!!!!!!!!!!!!!!!!check iii-1
              RC=C_WRITE(LID, WORK(IPCMP), (III-1)*RTOB)
              IF (RC .EQ. 0) THEN
                WRITE(6,*) 'error in writing luscus file!'
                CALL Abend()
              END IF
            ELSE
              write (LuVal) (Work(ipCMP+j),j=0,iii-1)
            END IF
            Call GetMem('TMP','FREE','REAL',ipCMP,mCoor)
           else
c no cut off
C              IF (ISLUSCUS .EQ. 1) THEN
C                 call dump_lusc(LID, WipOut,mCoor)
C           print *,'here'
c                RC=C_WRITE(LID, WIPOUT, MCOOR*RTOB) !!!!!!!!!!!!!!!!!!!!check MCOOR
c                IF (RC .EQ. 0) THEN
c                  WRITE(6,*) 'error in writing luscus file!'
c                  CALL Abend()
c                END IF
C              ELSE
                 write (LuVal) (WipOut(j),j=1,mCoor)
C              END IF
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

              if(ISLUSCUS.eq.0)
     &              write (LuVal,'(E18.12)') (WipOut(j),j=1,mCoor)
              endif
                 else
c extra output -
              if(ISLUSCUS.eq.0) write (LuVal,'(E18.12,3f8.4)')
     &       (WipOut(j),WCoor(1,j), WCoor(2,j), Wcoor(3,j) ,j=1,mCoor)
                 endif

               endif
CGG This is only for testing CASDFT functional. It will be restore.
CGG              write (LuVal,'(E10.4)') (Work(j),j=ipOut,ipOut+mCoor-1)
            endif
          endif
        endif

      if(isLine.eq.1.and.ISLUSCUS.eq.0) then
       do i=1,mCoor
        write(LuVal,'(3F10.6,22E20.12)')
     *  (WCoor(j,i),j=1,3),(WLine(j,i),j=1,nLine)
       enddo
      endif

*
      Return
      end
      SUBROUTINE PRTLUSENDGRID(LUVAL)
      CHARACTER LINE*128
      WRITE(LINE,'(A)') ' </INPORB>'
      CALL PRINTLINE(LUVAL, LINE,10,0)
      WRITE(LINE,'(A)') ' </GRID>'
      CALL PRINTLINE(LUVAL, LINE,8,0)
      END
