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

use, intrinsic :: iso_c_binding, only: c_null_char
use grid_it_globals, only: iBinary, isDebug, isLine, isLuscus, isTheOne, isUHF, LID, LID_ab, LuVal, LuVal_ab, TheName, Title1
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: INPORB
integer(kind=iwp) :: i, iiUHF, istatus, iPRGM, irecl, l, LuExtra, LuOrb, mm, RC
real(kind=wp) :: g
logical(kind=iwp) :: is_error
character(len=512) :: TMPLUS
character(len=256) :: FullName
character(len=306) :: RealName
character(len=64) :: Project
character(len=40) :: Env
character(len=12) :: Alpha
character(len=6) :: ShortName
character(len=2) :: ss
character, parameter :: Slash = '/'
integer(kind=iwp), external :: isFreeUnit
interface
  function lusopen(lid,fname,fname_len) bind(C,name='lusopen_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: lusopen
    integer(kind=MOLCAS_C_INT) :: lid, fname_len
    character(kind=c_char) :: fname(*)
  end function lusopen
end interface

LuOrb = isFreeUnit(46)
iPRGM = 0
call Chk_Vec_UHF(INPORB,LuOrb,iiUHF)
isUHF = (iiUHF == 1)
close(LuOrb)
LuVal_ab = -99999
if (isLuscus) then
  Alpha = '.lus'
else
  Alpha = '.grid'
end if
ShortName = 'M2MSI'
if (isUHF) then
  if (isLuscus) then
    Alpha = '_a.lus'
  else
    Alpha = '_a.grid'
  end if
end if
do iiUHF=0,merge(1,0,isUHF)
  if (iiUHF == 1) then
    if (isLuscus) then
      Alpha = '_b.lus'
    else
      Alpha = '_b.grid'
    end if
  end if
  Env = 'WorkDir '
  call getenvf(Env,FullName)
  if ((TheName /= ' ') .and. (iiUHF == 0)) then
    LuExtra = isFreeUnit(88)
    call molcas_open(LuExtra,'extra.prgm')
    iPRGM = 1
  end if
  if (isUHF) then
    if (isLuscus) then
      if (iiUHF == 0) ShortName = 'AM2L'
      if (iiUHF == 1) ShortName = 'BM2L'
    else
      if (iiUHF == 0) ShortName = 'AM2MSI'
      if (iiUHF == 1) ShortName = 'BM2MSI'
    end if
  end if
  if ((FullName == ' ') .or. (TheName(1:1) == ' ')) then
    RealName = ShortName

  else
    l = 1
    !ii = 1
    !fullname = outf
    do
      i = index(FullName(l:),Slash)
      if (i == 0) exit
      !ii = i+l
      l = l+i+1
    end do
    call getenvf('Project',Project)
    !Project = FullName(ii:)
    if ((TheName == 'NEW') .or. (TheName == 'new') .or. (TheName == 'New')) then
      do i=1,99
        write(ss,'(i2)') i
        if (i < 10) ss(1:1) = '0'
        RealName = trim(FullName)//Slash//trim(Project)//'.'//ss//Alpha
      end do

    else

      RealName = trim(FullName)//Slash//trim(Project)//'.'//trim(TheName)//Alpha
    end if
    write(u6,*) 'Grid file: ',trim(RealName)
  end if

  if (TheName /= ' ') then
    !LuExtra = isFreeUnit(88)
    !open(LuExtra,file='extra.prgm')
    write(LuExtra,'(a,a,a,a,a)') ' (file) ',trim(ShortName),' ',trim(RealName),'  rwsg'
    !close(LuExtra)
  end if
  if (iiUHF == 0) then
    LuVal = isFreeUnit(49)
    if (isLuscus) then
      ! FIXME: User can't define luscus input file name
      RC = -1
      if (Thename == ' ') then
        if (.not. isUHF) TMPLUS = 'LUSCUS'
        if (isUHF) TMPLUS = 'AM2L'
        mm = len_trim(TMPLUS)+1
        !write(u6,*) ' before 2 lusop',mm
        RC = lusopen(LID,trim(TMPLUS)//c_null_char,mm)
      else
        mm = len_trim(RealName)+1
        !write(u6,*) ' before 2 lusop',mm
        RC = lusopen(LID,trim(RealName)//c_null_char,mm)
      end if
      !rc = AixOpn(LID,'LUSCUS',.true.)
      if (RC /= 0) then
        write(u6,*) 'ERROR: Can''t open luscus file!'
        call Abend()
      end if
    else ! not luscus
      if (iBinary == 1) then
        call molcas_open_ext2(LuVal,RealName,'sequential','unformatted',istatus,.false.,irecl,'unknown',is_error)
        !open(unit=LuVal,access='sequential',form='unformatted',file=RealName)
        !write(u6,*) '** Create Grid file:',trim(RealName)
        write(LuVal) 'a'
        !if (isMOPack) then
        !  g = 2003.9_wp
        !  i = 0
        !  write(LuVal) g,nMOs,nShowMOs_,nCoor,nInc,nBlocks,merge(1,0,isCutOff),Cutoff,iiCoord
        !  write(LuVal) i
        !else
        g = 1999.0_wp
        write(LuVal) g
        !end if
        write(LuVal) Title1
      else if (iBinary == 0) then
        call molcas_open(LuVal,RealName)
        !open(unit=LuVal,file=RealName,Form='FORMATTED')
        if (isLine) then
          write(LuVal,'(a)') '# data in GNUplot format'
          exit
        end if
        !write(u6,*) '** Create Grid file (in ASCII format):',trim(RealName)
        if (isTheOne) then
          write(LuVal,'(a1)') '9'
        else
          write(LuVal,'(a1)') '0'
        end if
        if (isDebug) then
          write(Luval,'(a,a)') Title1,' DEBUG'
        else
          write(Luval,'(a)') Title1
        end if
      end if
    end if
  else ! iiUHF
    LuVal_ab = isFreeUnit(51)
    if (isLuscus) then
      ! FIXME: User can't define luscus input file name
      RC = -1
      if (Thename == ' ') then
        if (isUHF) TMPLUS = 'BM2L'
        mm = len_trim(TMPLUS)+1
        !write(u6,*) ' before 1 lusop',mm
        RC = lusopen(LID_ab,trim(TMPLUS)//c_null_char,mm)
      else
        mm = len_trim(RealName)+1
        !write(u6,*) ' before 1 lusop',mm
        RC = lusopen(LID_ab,trim(RealName)//c_null_char,mm)
      end if
      !rc = AixOpn(LID,'LUSCUS',.true.)
      if (RC /= 0) then
        write(u6,*) 'ERROR: Can''t open luscus file!'
        call Abend()
      end if
    else ! not luscus
      if (iBinary == 1) then
        call molcas_open_ext2(LuVal_ab,RealName,'sequential','unformatted',istatus,.false.,irecl,'unknown',is_error)
        !open(unit=LuVal_ab,access='sequential',form='unformatted',file=RealName)
        !write(u6,*) '** Create Grid file',trim(RealName)
        write(LuVal_ab) 'a'
        !if (isMOPack) then
        !  g = 2003.9_wp
        !  i = 0
        !  write(LuVal_ab) g
        !  write(LuVal_ab) i
        !else
        g = 1999.0_wp
        write(LuVal_ab) g
        !end if
        write(LuVal_ab) Title1
      else if (iBinary == 0) then
        call molcas_open(LuVal_ab,RealName)
        !open(unit=LuVal_ab,file=RealName,Form='FORMATTED')
        !write(u6,*) '** Create Grid file (in ASCII format):',trim(RealName)
        if (isTheOne) then
          write(LuVal_ab,'(a1)') '9'
        else
          !if (isMOPack) then
          !  write(LuVal_ab,'(a1)') '1'
          !  write(LuVal_ab,'(i5)') 0
          !else
          write(LuVal_ab,'(a1)') '0'
          !end if
        end if
        if (isDebug) then
          write(Luval_ab,'(a,a)') Title1,' DEBUG'
        else
          write(Luval_ab,'(a)') Title1
        end if
      end if
    end if
  end if
end do
if (iPRGM == 1) close(LuExtra)

return

end subroutine OpenGrid
