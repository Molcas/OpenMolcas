!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001-2005, Valera Veryazov                             *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine BasisType(Filename,inline,BasisTypes)

! translate basis file names

use BasisType_Mod, only: BasTypeAll, BasTypeCon, BasTypeNuc, BasTypeRel
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Filename
integer(kind=iwp), intent(in) :: inline
integer(kind=iwp), intent(out) :: BasisTypes(4)
integer(kind=iwp) :: i, i1, iflag, iL, ileft, irecl, istatus, iunit, j, jj, LenOrig
logical(kind=iwp) :: is_error, isKnown, found
character(len=256) :: DirName, OrigName, Line, LineComp
character(len=16) :: CONT, ALLE, RELA, NUCL
character(len=4) :: tmp
character(len=3) :: vCONT, vALLE, vRELA, vNUCL
integer(kind=iwp), external :: isFreeUnit, StrnLn

CONT = '#Contraction '
ALLE = '#AllElectron '
RELA = '#Hamiltonian '
NUCL = '#Nucleus '
vCONT = 'UNK'
vALLE = 'UNK'
vRELA = 'UNK'
vNUCL = 'UNK'
isKnown = .false.
if (inline == 1) then
  BasisTypes(:) = -1
  call SysWarnMsg('BasisType','inline basis is used','assuming all defaults for the basis types')
  return
end if
BasisTypes(3) = 0
iunit = isfreeunit(20)
ileft = StrnLn(FileName)
i = StrnLn(FileName)
found = .false.

do while ((.not. found) .and. (i > 1))
  if (FileName(i:i) == '/') then
    ileft = i
    found = .true.
  end if
  i = i-1
end do

if (.not. found) then
  ileft = StrnLn(FileName)
  i = StrnLn(FileName)
  do while ((.not. found) .and. (i > 1))
    if (FileName(i:i) == '_') then
      ileft = i
      found = .true.
    end if
    i = i-1
  end do
end if
DirName = Filename(1:ileft)
i = index(Filename,' ')
if (i == 0) i = len(Filename)
OrigName = Filename(ileft+1:i)
LenOrig = i-ileft

! first, check tbl.
call molcas_open_ext2(iunit,DirName(1:ileft)//'basistype.tbl','sequential','formatted',istatus,.false.,irecl,'old',is_error)
!open(iunit,file=DirName(1:ileft)//'basistype.tbl',form='formatted',iostat=istatus)
if (istatus /= 0) then
  close(iunit)
  call molcas_open_ext2(iunit,'BASLIB_basistype.tbl','sequential','formatted',istatus,.false.,irecl,'old',is_error)
  !open(iunit,file='BASLIB_basistype.tbl',form='formatted',iostat=istatus)
end if
if (istatus /= 0) then
  write(u6,*) 'basistype.tbl is not found'
  close(iunit)
else
  do
    read(iunit,'(a)',iostat=istatus) Line
    if (istatus /= 0) exit
    if (Line(1:1) == '#') cycle
    i = index(Line,OrigName(1:LenOrig-1))
    if ((i == 1) .and. (Line(LenOrig:LenOrig) == ' ')) exit
  end do
  if (istatus == 0) then
    i = LenOrig+1
    iL = 1
    iflag = 0
    LineComp = ' '
    do jj=i,len(Line)
      if (Line(jj:jj) == ' ') then
        if (iflag == 0) then
          iflag = 1
          LineComp(iL:iL) = ':'
          iL = iL+1
        end if
      else
        iflag = 0
        LineComp(iL:iL) = Line(jj:jj)
        iL = iL+1
      end if
    end do
    isKnown = .true.
    vCONT = LineComp(2:4)
    vALLE = LineComp(6:8)
    vRELA = LineComp(10:12)
    vNUCL = LineComp(14:16)
    if (vCONT == ' ') vCONT = 'UNK'
    if (vALLE == ' ') vALLE = 'UNK'
    if (vRELA == ' ') vRELA = 'UNK'
    if (vNUCL == ' ') vNUCL = 'UNK'
  end if
  close(iunit)
end if

! now let's check the header of the basis file

Found = .false.
if (OrigName(1:LenOrig-1) /= ' ') call f_Inquire(DirName(1:ileft)//OrigName(1:LenOrig-1),Found)
if (Found) then
  !write(u6,*) 'open >',DirName(1:ileft)//OrigName(1:LenOrig-1),'<'
  call molcas_open(iunit,DirName(1:ileft)//OrigName(1:LenOrig-1))
  do
    read(iunit,'(a)',iostat=istatus) Line
    if (istatus /= 0) exit
    if (Line(1:1) == '/') then
      ! we already reached body of the file.
      LineComp = ':'//vCONT//':'//vALLE//':'//vRELA//':'//vNUCL//':'
      isKnown = .true.
      exit
    end if
    ! finish reading, let's create the return string

    !write(u6,*) 'line=',Line
    if (index(Line,CONT(1:index(CONT,' '))) == 1) then
      i = index(Line,' ')
      do j=i,len(Line)
        if (Line(j:j) /= ' ') exit
      end do
      vCONT = Line(j:j+3)
    end if
    if (index(Line,ALLE(1:index(ALLE,' '))) == 1) then
      i = index(Line,' ')
      do j=i,len(Line)
        if (Line(j:j) /= ' ') exit
      end do
      vALLE = Line(j:j+3)
    end if
    if (index(Line,RELA(1:index(RELA,' '))) == 1) then
      i = index(Line,' ')
      do j=i,len(Line)
        if (Line(j:j) /= ' ') exit
      end do
      vRELA = Line(j:j+3)
    end if

    if (index(Line,NUCL(1:index(NUCL,' '))) == 1) then
      i = index(Line,' ')
      do j=i,len(Line)
        if (Line(j:j) /= ' ') exit
      end do
      vNUCL = Line(j:j+3)
    end if

  end do

  ! well, now use tbl

  close(iunit)
end if

if (.not. isKnown) LineComp = ':UNK:UNK:UNK:UNK:'
!write(u6,*) 'DEBUG=',LineComp
tmp = LineComp(2:5)
i1 = index(BasTypeCon,tmp)
if ((i1 == 0) .or. (tmp == 'UNK:')) then
  BasisTypes(1) = -1
else
  i1 = i1/4+1
  BasisTypes(1) = i1
end if

tmp = LineComp(6:9)
i1 = index(BasTypeAll,tmp)
if ((i1 == 0) .or. (tmp == 'UNK:')) then
  BasisTypes(2) = -1
else
  i1 = i1/4+1
  BasisTypes(2) = i1
  ! hack to map YES to AE_, and NO_ to NAE
  if (BasisTypes(2) == 3) BasisTypes(2) = 1
  if (BasisTypes(2) == 4) BasisTypes(2) = 2
end if

tmp = LineComp(10:13)
i1 = index(BasTypeRel,tmp)
if ((i1 == 0) .or. (tmp == 'UNK:')) then
  BasisTypes(3) = -1
else
  i1 = i1/4+1
  BasisTypes(3) = i1
end if

tmp = LineComp(14:17)
i1 = index(BasTypeNuc,tmp)
if ((i1 == 0) .or. (tmp == 'UNK:')) then
  BasisTypes(4) = -1
else
  i1 = i1/4+1
  BasisTypes(4) = i1
end if

return

end subroutine BasisType
