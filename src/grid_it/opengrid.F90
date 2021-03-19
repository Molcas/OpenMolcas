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
      Subroutine OpenGrid(INPORB)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************
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
!      logical exist
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
!       ii=1
!       fullname=outf
888     i=index(FullName(l:),Slash)
        if(i.gt.0) then
!          ii=i+l
          l=l+i+1
          goto 888
        endif
        Call getenvf('Project',Project)
!        Project=FullName(ii:)
       if(TheName.eq.'NEW'.or.TheName.eq.'new'.or.TheName.eq.'New') then
        do i=1,99
         write(ss,'(i2)') i
         if(i.lt.10) ss(1:1)='0'
         RealName=FullName(1:index(FullName,' ')-1)//Slash              &
     &    //Project(1:index(Project,' ')-1)//'.'                        &
     &    //ss//Alpha
        enddo

       else

         RealName=FullName(1:index(FullName,' ')-1)//Slash              &
     &    //Project(1:index(Project,' ')-1)//'.'                        &
     &    //TheName(1:index(TheName,' ')-1)//Alpha
       endif
       write(6,*) 'Grid file: ',RealName(:mylen(RealName))
      endif

      if(TheName.ne.' ') then
!         open(88,file='extra.prgm')
         write(88,'(a,a,a)') ' (file) M2MSI ',                          &
     &        RealName(1:index(RealName,' ')),'  rwsg'
!         close(88)
      endif
      if(iiUHF.eq.0) then
      LuVal=isFreeUnit(49)
      if (ISLUSCUS .EQ. 1) THEN
! FIXME: User can't define luscus input file name
        RC=-1
        if(Thename.eq.' ') then
        if(isUHF.eq.0) TMPLUS(1:)="LUSCUS"
        if(isUHF.eq.1.and.iiUHF.eq.0) TMPLUS(1:8)="alph.lus"
        if(isUHF.eq.1.and.iiUHF.eq.1) TMPLUS(1:8)="beta.lus"
        mm=mylen(TMPLUS)
!        print *,' before 2 lusop', mm
        RC=lusopen(LID,TMPLUS,mm)
        else
        mm=mylen(RealName)
!        print *,' before 2 lusop', mm
        RC=lusopen(LID,RealName,mm)
        endif
!        rc=AixOpn(LID,"LUSCUS",.TRUE.)
        IF (RC .NE. 0) THEN
          write(6,*) 'ERROR: Can''t open luscus file!'
          CALL Abend()
        END IF
      ELSE ! not luscus
      if(isBinary.eq.1) Then
      call molcas_open_ext2(LuVal,RealName,'sequential',                &
     & 'unformatted',iostat,.false.,irecl,'unknown',is_error)
!        open(unit=LuVal,access='sequential',
!     ,       form='unformatted', file=RealName)
!        write(6,*) '** Create Grid file:',
!     &        RealName(1:index(RealName,' '))
        write(LuVal) 'a'
!        if (imoPack .ne. 0) then
!          g=2003.9
!          i=0
!          write (LuVal) g
!     +         nMOs, nShowMOs_, nCoor, nInc, nBlocks,
!     +         isCutOff, Cutoff, iiCoord
!          write (LuVal) i
!        else
          g=1999.0
          write(LuVal) g
!        endif
        write(LuVal) Title1
      endif

      if(isBinary.eq.0) Then
        call molcas_open(LuVal,RealName)
!        open(unit=LuVal,file=RealName,Form='FORMATTED')
        if(isLine.eq.1) then
          Write(LuVal,'(a)') '# data in GNUplot format'
          goto 999
        endif
!        write(6,*) '** Create Grid file (in ASCII format):',
!     &        RealName(1:index(RealName,' '))
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
! FIXME: User can't define luscus input file name
        if(Thename.eq.' ') then
        if(isUHF.eq.1.and.iiUHF.eq.0) TMPLUS="AM2L"
        if(isUHF.eq.1.and.iiUHF.eq.1) TMPLUS="BM2L"
        mm=6
!        print *,' before 1 lusop', mm
        RC=lusopen(LID_ab,TMPLUS,mm)
        else
        mm=mylen(RealName)
!        print *,' before lusop', mm
        RC=lusopen(LID_ab,RealName,mm)
        endif
!        rc=AixOpn(LID,"LUSCUS",.TRUE.)
        IF (RC .NE. 0) THEN
          write(6,*) 'ERROR: Can''t open luscus file!'
          CALL Abend()
        END IF
      endif
      LuVal_ab=isFreeUnit(51)


      if(isBinary.eq.1) Then
      call molcas_open_ext2(LuVal_ab,RealName,'sequential',             &
     &  'unformatted',iostat,.false.,irecl,'unknown',is_error)
!        open(unit=LuVal_ab,access='sequential',
!     ,       form='unformatted', file=RealName)
!        write(6,*) '** Create Grid file',
!     &        RealName(1:index(RealName,' '))
        write(LuVal_ab) 'a'
!        if (imoPack .ne. 0) then
!          g=2003.9
!          i=0
!          write (LuVal_ab) g
!          write (LuVal_ab) i
!        else
          g=1999.0
          write(LuVal_ab) g
!        endif
        write(LuVal_ab) Title1
      endif

      if(isBinary.eq.0) Then
        call molcas_open(LuVal_ab,RealName)
!        open(unit=LuVal_ab,file=RealName,Form='FORMATTED')
!        write(6,*) '** Create Grid file (in ASCII format):',
!     &        RealName(1:index(RealName,' '))
          if (isTheOne.eq.1) then
            write(LuVal_ab,'(a1)') '9'
          else
!            if (imoPack .ne. 0) then
!              write (LuVal_ab,'(a1)') '1'
!              write (LuVal_ab,'(i5)') 0
!            else
              write(LuVal_ab,'(a1)') '0'
!            endif
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
