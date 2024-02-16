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

subroutine fetch_init_const(nneq,neqv,nmax,exch,nLoc,nCenter,nT,nH,nTempMagn,nDir,nDirZee,nMult,nPair,MxRank1,MxRank2, &
                            old_aniso_format,iReturn)
! this routine looks into the file "single_aniso.input" for the "RESTart" keyword

use Constants, only: Zero
use Definitions, only: u5, u6

implicit none
#include "warnings.h"
integer, intent(out) :: nneq, neqv, nmax, exch, nLoc, nCenter, nT, nH, nTempMagn, nDir, nDirZee, nMult, nPair, MxRank1, MxRank2, &
                        iReturn
! local variables:
integer :: NMAXC
parameter(NMAXC=99)
integer :: i, j, linenr, nTempMagn_HEXP, nTempMagn_TMAG, nH_HEXP, nH_HINT, nT_TEXP, nT_TINT
integer :: neqA(NMAXC), nexchA(NMAXC)
integer :: sfs_check(NMAXC)
integer :: sos_check(NMAXC)
integer :: imaxrank(NMAXC,2)
integer :: idummy
integer :: irank1, irank2, iline
real(kind=8) :: rdummy
real(kind=8) :: TempMagn(NMAXC)
logical :: ab_initio_all
logical :: KeyCoor, KeyPair, KeyHEXP, KeyTEXP
!logical :: KeyHINT, KeyTINT, KeyTMAG, KeyMLTP, KeyMVEC, KeyNNEQ, KeyZEEM, KeyITOJ
integer :: LUANISO, Isfreeunit
character :: itype(NMAXC)
character(len=280) :: line
character(len=180) :: namefile_aniso
logical :: ifHDF
logical :: DBG
external Isfreeunit
logical, intent(in) :: old_aniso_format

iReturn = 0
nH = 0
nT = 0
nTempMagn = 1
nneq = 0
neqv = 0
nmax = 0
exch = 0
nCenter = 0
nDirZee = 0
nDir = 0
nMult = 0
nLoc = 0
nPair = 0
MxRank1 = 0
MxRank2 = 0
luaniso = 0
rdummy = Zero
neqA(1:nmaxc) = 0
nexchA(1:nmaxc) = 0
sfs_check(1:nmaxc) = 0
sos_check(1:nmaxc) = 0
ab_initio_all = .false.
itype(1:nmaxc) = ' '
imaxrank(1:nmaxc,1:2) = 0

!namefile_aniso = ''
ifHDF = .false.

DBG = .false.

!KeyNNEQ = .false.
KeyPair = .false.
KeyCoor = .false.
KeyHEXP = .false.
KeyTEXP = .false.
!KeyTMAG = .false.
!KeyTINT = .false.
!KeyHINT = .false.
!KeyMLTP = .false.
!KeyMVEC = .false.
!KeyZEEM = .false.
!KeyITOJ = .false.
nH_HEXP = 0
nH_HINT = 0
nT_TEXP = 0
nT_TINT = 0
nTempMagn_HEXP = 0
nTempMagn_TMAG = 0
!=========== End of default settings====================================
rewind(u5)
50 read(u5,'(A280)',end=998) LINE
call NORMAL(LINE)
if (LINE(1:11) /= '&POLY_ANISO') Go To 50
LINENR = 0
100 read(u5,'(A280)',end=998) line
LINENR = LINENR+1
call NORMAL(LINE)
if (LINE(1:1) == '*') Go To 100
if (LINE == ' ') Go To 100
if ((LINE(1:4) /= 'NNEQ') .and. (LINE(1:4) /= 'TEXP') .and. (LINE(1:4) /= 'HEXP') .and. (LINE(1:4) /= 'END ') .and. &
    (LINE(1:4) /= '    ') .and. (LINE(1:4) /= 'HINT') .and. (LINE(1:4) /= 'TINT') .and. (LINE(1:4) /= 'TMAG') .and. &
    (LINE(1:4) /= 'MVEC') .and. (LINE(1:4) /= 'ZEEM') .and. (LINE(1:4) /= 'MLTP') .and. (LINE(1:4) /= 'LIN1') .and. &
    (LINE(1:4) /= 'LIN3') .and. (LINE(1:4) /= 'LIN9') .and. (LINE(1:4) /= 'PAIR') .and. (LINE(1:4) /= 'ALIN') .and. &
    (LINE(1:4) /= 'COOR') .and. (LINE(1:4) /= 'ITOJ')) Go To 100
if ((LINE(1:4) == 'END ') .or. (LINE(1:4) == '    ')) Go To 200

if (line(1:4) == 'NNEQ') then

  !KeyNNEQ = .true.
  read(u5,*) nneq,ab_initio_all,ifHDF

  if (DBG) write(u6,*) nneq,ab_initio_all,ifHDF

  if (nneq < 0) then
    write(u6,'(A)') 'nneq<0! Must be positive!'
    call Quit_OnUserError()
  else if (nneq == 0) then
    write(u6,'(A)') 'nneq=0! Must be larger than zero!'
    call Quit_OnUserError()
  else if (nneq > NMAXC) then
    write(u6,'(A)') 'nneq>99! Must be smaller than this!'
    call Quit_OnUserError()
  end if

  read(u5,*) (neqA(i),i=1,nneq)
  if (DBG) write(u6,*) (neqA(i),i=1,nneq)

  do i=1,nneq
    if (neqA(i) < 0) then
      write(u6,'(A,i2,A)') 'neq(',i,')<0! Must be positive!'
      call Quit_OnUserError()
    else if (neqA(i) == 0) then
      write(u6,'(A,i2,A)') 'neq(',i,')=0! Must be larger than zero!'
      call Quit_OnUserError()
    end if
  end do

  read(u5,*) (nexchA(i),i=1,nneq)
  if (DBG) write(u6,*) (nexchA(i),i=1,nneq)

  do i=1,nneq
    if (nexchA(i) < 0) then
      write(u6,'(A,i2,A)') 'nexch(',i,')<0! Must be positive!'
      call Quit_OnUserError()
    else if (nexchA(i) == 0) then
      write(u6,'(A,i2,A)') 'nexch(',i,')=0! Must be larger than zero!'
      call Quit_OnUserError()
    end if
  end do

  neqv = 0
  neqv = maxval(neqA(1:nneq))
  if (DBG) write(u6,*) 'neqv = ',neqv
  nmax = 0
  nmax = maxval(nexchA(1:nneq))
  if (DBG) write(u6,*) 'nmax = ',nmax

  ! compute "exch"
  exch = 1
  do i=1,nneq
    do j=1,neqA(i)
      exch = exch*nexchA(i)
    end do
  end do
  if (DBG) write(u6,*) 'exch=',exch

  if (.not. ab_initio_all) then
    read(u5,*,err=997) (itype(i),i=1,Nneq)
  else
    do i=1,nneq
      itype(i) = 'A'
    end do
  end if !ab_initio_all

  ! check the maximal number of local spin-orbit states
  sos_check = 0
  sfs_check = 0
  nCenter = 0
  do i=1,NNEQ
    nCenter = nCenter+neqA(I)
  end do

  do i=1,NNEQ
    if (itype(i) == 'A') then
      if (ifHDF) then
        ! generating the name of the "aniso_input file for
        ! each center. Maxmimum 10 centers. CHAR(48)=0 (zero)
        if (i < 10) then
          write(namefile_aniso,'(4A)') 'aniso_hdf_',char(48+mod(int(i),10)),'.input'
        else if ((i >= 10) .and. (i <= 99)) then
          write(namefile_aniso,'(4A)') 'aniso_hdf_',char(48+mod(int((i)/10),10)),char(48+mod(int(i),10)),'.input'
        end if

#       ifdef _HDF5_
        call read_hdf5_init(NAMEFILE_ANISO,sfs_check(I),sos_check(I))
        if (DBG) write(u6,*) ' sfs(I) ',sfs_check(I),' sos(I) ',sos_check(I)
#       else
        call WarningMessage(2,'File '//trim(NAMEFILE_ANISO)//' cannot be opened. Molcas was compiled without HDF5 option.')
        call Quit_OnUserError()
#       endif
      else
        ! generating the name of the "aniso_input file for
        ! each center. Maxmimum 10 centers. CHAR(48)=0 (zero)
        if (i < 10) then
          write(namefile_aniso,'(4A)') 'aniso_',char(48+mod(int(i),10)),'.input'
        else if ((i >= 10) .and. (i <= 99)) then
          write(namefile_aniso,'(4A)') 'aniso_',char(48+mod(int((i)/10),10)),char(48+mod(int(i),10)),'.input'
        end if
        LUANISO = Isfreeunit(20)
        call molcas_open(LUANISO,NAMEFILE_ANISO)

        if (old_aniso_format) then
          read(LUANISO,*) sfs_check(I),sos_check(I)
          if (DBG) write(u6,*) ' sfs(I) ',sfs_check(I),' sos(I) ',sos_check(I)
        else
          call read_nss(LUANISO,sos_check(i),dbg)
          call read_nstate(LUANISO,sfs_check(i),dbg)
          if (DBG) write(u6,*) ' sfs(I) ',sfs_check(I),' sos(I) ',sos_check(I)
        end if
        close(LUANISO)
      end if ! ifHDF

    else if ((itype(i) == 'B') .or. (itype(i) == 'C')) then
      sfs_check(I) = 1
      sos_check(I) = NexchA(i)
    end if

  end do ! NNEQ

  nLoc = maxval(sos_check(1:nneq))

  LINENR = LINENR+3
  Go To 100
end if

if (line(1:4) == 'TEXP') then

  KeyTexp = .true.
  read(u5,*) nT_TEXP

  if (nT_TEXP <= 0) then
    call WarningMessage(2,'TEXP: Number of temperature points <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'TEXP:: = nT_TEXP',nT_TEXP

  LINENR = LINENR+1
  Go To 100
end if

if (line(1:4) == 'HEXP') then

  KeyHexp = .true.
  read(u5,*) nTempMagn_HEXP,(TempMagn(i),i=1,nTempMagn)
  read(u5,*) nH_HEXP

  if (nH_HEXP <= 0) then
    call WarningMessage(2,'HEXP: Number of field points <= 0! ')
    call Quit_OnUserError()
  end if

  if (nTempMagn_HEXP <= 0) then
    call WarningMessage(2,'HEXP: Number of temperature points <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'HEXP:: = nH_HEXP ',nH_HEXP
  if (DBG) write(u6,*) 'HEXP:: = nTempMagn_HEXP ',nTempMagn_HEXP
  LINENR = LINENR+2
  Go To 100
end if

if (line(1:4) == 'HINT') then

  !KeyHINT = .true.
  read(u5,*) rdummy,rdummy,nH_HINT

  if (nH_HINT <= 0) then
    call WarningMessage(2,'HINT: Number of field points <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'HINT:: = nH_HINT ',nH_HINT
  LINENR = LINENR+1
  Go To 100
end if

if (line(1:4) == 'TINT') then

  !KeyTINT = .true.
  read(u5,*) rdummy,rdummy,nT_TINT

  if (nT_TINT <= 0) then
    call WarningMessage(2,'TINT: Number of temperature points <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'TINT:: = nT_TINT ',nT_TINT
  LINENR = LINENR+1
  Go To 100
end if

if (line(1:4) == 'TMAG') then

  !KeyTMAG = .true.
  read(u5,*) nTempMagn_TMAG

  if (nTempMagn_TMAG <= 0) then
    call WarningMessage(2,'TMAG: Number of temperatureMAGN points <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'TMAG:: = nTempMagn_TMAG ',nTempMagn_TMAG
  LINENR = LINENR+1
  Go To 100
end if

if (line(1:4) == 'MVEC') then

  !KeyMVEC = .true.
  read(u5,*) nDir

  if (nDir <= 0) then
    call WarningMessage(2,'MVEC: Number of nDir points <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'MVEC:: = nDir ',nDir
  LINENR = LINENR+1
  Go To 100
end if

if (line(1:4) == 'ZEEM') then
  !KeyZEEM = .false.
  read(u5,*) nDirZee

  if (nDirZee <= 0) then
    call WarningMessage(2,'ZEEM: Number of nDirZee points <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'ZEEM:: = nDirZee ',nDirZee
  LINENR = LINENR+1
  Go To 100
end if

if (line(1:4) == 'MLTP') then

  !KeyMLTP = .true.
  read(u5,*) nMult

  if (nMult <= 0) then
    call WarningMessage(2,'MLTP: Number of multiplets <= 0! ')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'MLTP:: =nMult ',nMult
  LINENR = LINENR+1
  Go To 100
end if

if ((LINE(1:4) == 'LIN9') .or. (LINE(1:4) == 'LIN3') .or. (LINE(1:4) == 'LIN1') .or. (LINE(1:4) == 'ALIN') .or. &
    (LINE(1:4) == 'PAIR') .or. (LINE(1:4) == 'ITOJ')) then

  KeyPair = .true.
  read(u5,*,err=997) nPair

  if (nPair <= 0) then
    call WarningMessage(2,'PAIR OR LINx:: Number of interacting pairs <= 0!')
    call Quit_OnUserError()
  end if

  if (DBG) write(u6,*) 'PAIR:: =nPair ',nPair

  if (LINE(1:4) == 'ITOJ') then
    iline = 0
    do i=1,npair
      imaxrank(i,1) = 0
      imaxrank(i,2) = 0

      read(u5,*,err=997) idummy,idummy,imaxrank(i,1),imaxrank(i,2)
      do irank1=1,2*imaxrank(i,1)+1
        do irank2=1,2*imaxrank(i,2)+1
          read(u5,*,err=997) idummy,idummy,idummy,idummy,rdummy,rdummy
          iline = iline+1
        end do
      end do
    end do ! i
    MxRank1 = maxval(imaxrank(1:npair,1))
    MxRank2 = maxval(imaxrank(1:npair,2))
    LINENR = LINENR+iline
  end if
  LINENR = LINENR+1
  Go To 100
end if

if (line(1:4) == 'COOR') then
  KeyCoor = .true.
  LINENR = LINENR+1
  Go To 100
end if

200 continue

if ((.not. KeyPair) .and. KeyCoor) nPair = nCenter*(nCenter-1)/2

if (KeyHexp) then
  nTempMagn = nTempMagn_HEXP
  nH = nH_HEXP
else
  nTempMagn = nTempMagn_TMAG
  nH = nH_HINT
end if

if (KeyTEXP) then
  nT = nT_TEXP
else
  nT = nT_TINT
end if

! in case the user did not set up some of the above keywords
! assume the following default ones:
if (nMult == 0) nMult = 1
if (nT == 0) nT = 31
if (nH == 0) nH = 11
if (nTempMagn == 0) nTempMagn = 1

! preliminary check the values:
if (DBG) then
  write(u6,*) 'nneq     =',nneq
  write(u6,*) 'neqv     =',neqv
  write(u6,*) 'exch     =',exch
  write(u6,*) 'nLoc     =',nLoc
  write(u6,*) 'nmax     =',nmax
  write(u6,*) 'nCenter  =',nCenter
  write(u6,*) 'nT       =',nT
  write(u6,*) 'nH       =',nH
  write(u6,*) 'nTempMagn=',ntempMagn
  write(u6,*) 'nDir     =',nDir
  write(u6,*) 'nDirZee  =',nDirZee
  write(u6,*) 'nMult    =',nMult
  write(u6,*) 'nPair    =',nPair
  write(u6,*) 'MxRank1  =',MxRank1
  write(u6,*) 'MxRank2  =',MxRank2
end if

Go To 190
!------ errors ------------------------------
997 continue
write(u6,*) ' READIN: Error reading "poly_aniso.input" '
write(u6,*) ' near line nr.',LINENR+1
Go To 999
998 continue
write(u6,*) ' READIN: Unexpected End of input file.'
999 continue
call quit(_RC_INPUT_ERROR_)

190 continue
return
#ifdef _WARNING_WORKAROUND_
if (.false.) then
  call Unused_integer(idummy)
  call Unused_real(rdummy)
  call Unused_real_array(TempMagn)
end if
#endif

end subroutine fetch_init_const
