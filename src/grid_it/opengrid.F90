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

subroutine OpenGrid(INPORB)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "grid.fh"
character FullName*256
character RealName*306
character ss*2
character Env*40
character INPORB*(*)
character(LEN=512) TMPLUS
integer RC
!logical exist
logical is_error
character Slash
character*12 Alpha
integer LUSOPEN
external LUSOPEN

Slash = '/'

LuOrb = isFreeUnit(46)
iPRGM = 0
call Chk_Vec_UHF(INPORB,LuOrb,isUHF)
close(LuOrb)
LuVal_ab = -99999
if (isLuscus == 1) then
  Alpha = '.lus'
else
  Alpha = '.grid'
end if
if (isUHF == 1) then
  if (isLuscus == 1) then
    Alpha = '_a.lus'
  else
    Alpha = '_a.grid'
  end if
end if
do iiUHF=0,isUHF
  if (iiUHF == 1) then
    if (isLuscus == 1) then
      Alpha = '_b.lus'
    else
      Alpha = '_b.grid'
    end if
  end if
  Env = 'WorkDir '
  call getenvf(Env,FullName)
  iPRGM = 0
  if (TheName /= ' ') then
    call molcas_open(88,'extra.prgm')
    iPRGM = 1
  end if
  if (FullName == ' ' .or. TheName(1:1) == ' ') then
    RealName = 'M2MSI'
    if (isUHF == 1) then
      if (isLuscus == 1) then
        if (iiUHF == 0) RealName = 'AM2L'
        if (iiUHF == 1) RealName = 'BM2L'
      else
        if (iiUHF == 0) RealName = 'AM2MSI'
        if (iiUHF == 1) RealName = 'BM2MSI'
      end if
    end if

  else
    l = 1
    !ii = 1
    !fullname = outf
    888 i = index(FullName(l:),Slash)
    if (i > 0) then
      !ii = i+l
      l = l+i+1
      goto 888
    end if
    call getenvf('Project',Project)
    !Project = FullName(ii:)
    if (TheName == 'NEW' .or. TheName == 'new' .or. TheName == 'New') then
      do i=1,99
        write(ss,'(i2)') i
        if (i < 10) ss(1:1) = '0'
        RealName = FullName(1:index(FullName,' ')-1)//Slash//Project(1:index(Project,' ')-1)//'.'//ss//Alpha
      end do

    else

      RealName = FullName(1:index(FullName,' ')-1)//Slash//Project(1:index(Project,' ')-1)//'.'//TheName(1:index(TheName,' ')-1)// &
                 Alpha
    end if
    write(6,*) 'Grid file: ',RealName(:mylen(RealName))
  end if

  if (TheName /= ' ') then
    !open(88,file='extra.prgm')
    write(88,'(a,a,a)') ' (file) M2MSI ',RealName(1:index(RealName,' ')),'  rwsg'
    !close(88)
  end if
  if (iiUHF == 0) then
    LuVal = isFreeUnit(49)
    if (ISLUSCUS == 1) then
      ! FIXME: User can't define luscus input file name
      RC = -1
      if (Thename == ' ') then
        if (isUHF == 0) TMPLUS(1:) = 'LUSCUS'
        if (isUHF == 1 .and. iiUHF == 0) TMPLUS(1:8) = 'alph.lus'
        if (isUHF == 1 .and. iiUHF == 1) TMPLUS(1:8) = 'beta.lus'
        mm = mylen(TMPLUS)
        !write(6,*) ' before 2 lusop', mm
        RC = lusopen(LID,TMPLUS,mm)
      else
        mm = mylen(RealName)
        !write(6,*) ' before 2 lusop', mm
        RC = lusopen(LID,RealName,mm)
      end if
      !rc = AixOpn(LID,'LUSCUS',.TRUE.)
      if (RC /= 0) then
        write(6,*) 'ERROR: Can''t open luscus file!'
        call Abend()
      end if
    else ! not luscus
      if (isBinary == 1) then
        call molcas_open_ext2(LuVal,RealName,'sequential','unformatted',iostat,.false.,irecl,'unknown',is_error)
        !open(unit=LuVal,access='sequential',form='unformatted',file=RealName)
        !write(6,*) '** Create Grid file:',RealName(1:index(RealName,' '))
        write(LuVal) 'a'
        !if (imoPack /= 0) then
        !  g = 2003.9
        !  i = 0
        !  write(LuVal) g,nMOs,nShowMOs_,nCoor,nInc,nBlocks,isCutOff,Cutoff,iiCoord
        !  write(LuVal) i
        !else
        g = 1999.0
        write(LuVal) g
        !end if
        write(LuVal) Title1
      end if

      if (isBinary == 0) then
        call molcas_open(LuVal,RealName)
        !open(unit=LuVal,file=RealName,Form='FORMATTED')
        if (isLine == 1) then
          write(LuVal,'(a)') '# data in GNUplot format'
          goto 999
        end if
        !write(6,*) '** Create Grid file (in ASCII format):',RealName(1:index(RealName,' '))
        if (isTheOne == 1) then
          write(LuVal,'(a1)') '9'
        else
          write(LuVal,'(a1)') '0'
        end if
        if (isDebug == 0) then
          write(Luval,'(a)') Title1
        else
          write(Luval,'(a,a)') Title1,' DEBUG'
        end if
      end if
    end if
  else ! iiUHF
    if (ISLUSCUS == 1) then
      ! FIXME: User can't define luscus input file name
      if (Thename == ' ') then
        if (isUHF == 1 .and. iiUHF == 0) TMPLUS = 'AM2L'
        if (isUHF == 1 .and. iiUHF == 1) TMPLUS = 'BM2L'
        mm = 6
        write(6,*) ' before 1 lusop', mm
        RC = lusopen(LID_ab,TMPLUS,mm)
      else
        mm = mylen(RealName)
        write(6,*) ' before lusop', mm
        RC = lusopen(LID_ab,RealName,mm)
      end if
      !rc = AixOpn(LID,'LUSCUS',.TRUE.)
      if (RC /= 0) then
        write(6,*) 'ERROR: Can''t open luscus file!'
        call Abend()
      end if
    end if
    LuVal_ab = isFreeUnit(51)

    if (isBinary == 1) then
      call molcas_open_ext2(LuVal_ab,RealName,'sequential','unformatted',iostat,.false.,irecl,'unknown',is_error)
      !open(unit=LuVal_ab,access='sequential',form='unformatted',file=RealName)
      !write(6,*) '** Create Grid file',RealName(1:index(RealName,' '))
      write(LuVal_ab) 'a'
      !if (imoPack /= 0) then
      !  g = 2003.9
      !  i = 0
      !  write(LuVal_ab) g
      !  write(LuVal_ab) i
      !else
      g = 1999.0
      write(LuVal_ab) g
      !end if
      write(LuVal_ab) Title1
    end if

    if (isBinary == 0) then
      call molcas_open(LuVal_ab,RealName)
      !open(unit=LuVal_ab,file=RealName,Form='FORMATTED')
      !write(6,*) '** Create Grid file (in ASCII format):',RealName(1:index(RealName,' '))
      if (isTheOne == 1) then
        write(LuVal_ab,'(a1)') '9'
      else
        !if (imoPack /= 0) then
        !  write(LuVal_ab,'(a1)') '1'
        !  write(LuVal_ab,'(i5)') 0
        !else
        write(LuVal_ab,'(a1)') '0'
        !end if
      end if
      if (isDebug == 0) then
        write(Luval_ab,'(a)') Title1
      else
        write(Luval_ab,'(a,a)') Title1,' DEBUG'
      end if
    end if
  end if
end do
999 continue
if (iPRGM == 1) close(88)

return

end subroutine OpenGrid
