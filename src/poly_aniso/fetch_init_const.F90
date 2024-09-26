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

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp, u5, u6

implicit none
integer(kind=iwp), intent(out) :: nneq, neqv, nmax, exch, nLoc, nCenter, nT, nH, nTempMagn, nDir, nDirZee, nMult, nPair, MxRank1, &
                                  MxRank2, iReturn
logical(kind=iwp), intent(in) :: old_aniso_format
integer(kind=iwp), parameter :: NMAXC = 99
integer(kind=iwp) :: i, idummy, iline, imaxrank(NMAXC,2), irank1, irank2, istatus, linenr, LUANISO, neqA(NMAXC), nexchA(NMAXC), &
                     nH_HEXP, nH_HINT, nT_TEXP, nT_TINT, nTempMagn_HEXP, nTempMagn_TMAG, sfs_check(NMAXC), sos_check(NMAXC)
real(kind=wp) :: rdummy, TempMagn(NMAXC)
logical(kind=iwp) :: ab_initio_all, ifHDF, KeyCoor, KeyHEXP, KeyPair, KeyTEXP !, KeyHINT, KeyITOJ, KeyMLTP, KeyMVEC, KeyNNEQ, &
                     !KeyTINT, KeyTMAG, KeyZEEM
character(len=280) :: line
character(len=180) :: namefile_aniso
character :: itype(NMAXC)
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: dbg = _DBG_
integer(kind=iwp), external :: Isfreeunit

#include "macros.fh"

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
neqA(:) = 0
nexchA(:) = 0
sfs_check(:) = 0
sos_check(:) = 0
ab_initio_all = .false.
itype(:) = ' '
imaxrank(:,:) = 0

!namefile_aniso = ''
ifHDF = .false.

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
do
  read(u5,'(A280)',iostat=istatus) LINE
  if (istatus < 0) call Error(2)
  call NORMAL(LINE)
  if (LINE(1:11) == '&POLY_ANISO') exit
end do
LINENR = 0
do
  read(u5,'(A280)',iostat=istatus) line
  if (istatus < 0) call Error(2)
  LINENR = LINENR+1
  call NORMAL(LINE)
  if ((LINE(1:1) == '*') .or. (LINE == ' ')) cycle

  select case (LINE(1:4))

    case ('END ','    ')
      exit

    case ('NNEQ')

      !KeyNNEQ = .true.
      read(u5,*) nneq,ab_initio_all,ifHDF

#     ifdef _DEBUGPRINT_
      write(u6,*) nneq,ab_initio_all,ifHDF
#     endif

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
#     ifdef _DEBUGPRINT_
      write(u6,*) (neqA(i),i=1,nneq)
#     endif

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
#     ifdef _DEBUGPRINT_
      write(u6,*) (nexchA(i),i=1,nneq)
#     endif

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
#     ifdef _DEBUGPRINT_
      write(u6,*) 'neqv = ',neqv
#     endif
      nmax = 0
      nmax = maxval(nexchA(1:nneq))
#     ifdef _DEBUGPRINT_
      write(u6,*) 'nmax = ',nmax
#     endif

      ! compute "exch"
      exch = 1
      do i=1,nneq
        exch = exch*nexchA(i)**neqA(i)
      end do
#     ifdef _DEBUGPRINT_
      write(u6,*) 'exch=',exch
#     endif

      if (.not. ab_initio_all) then
        read(u5,*,iostat=istatus) (itype(i),i=1,Nneq)
        if (istatus /= 0) call Error(1)
      else
        itype(1:nneq) = 'A'
      end if !ab_initio_all

      ! check the maximal number of local spin-orbit states
      sos_check = 0
      sfs_check = 0
      nCenter = sum(neqA(1:NNEQ))

      do i=1,NNEQ
        if (itype(i) == 'A') then
          if (ifHDF) then
            ! generating the name of the "aniso_input file for
            ! each center. Maxmimum 10 centers. CHAR(48)=0 (zero)
            if (i < 10) then
              write(namefile_aniso,'(4A)') 'aniso_hdf_',char(48+mod(i,10)),'.input'
            else if ((i >= 10) .and. (i <= 99)) then
              write(namefile_aniso,'(4A)') 'aniso_hdf_',char(48+mod(i/10,10)),char(48+mod(i,10)),'.input'
            end if

#           ifdef _HDF5_
            call read_hdf5_init(NAMEFILE_ANISO,sfs_check(I),sos_check(I))
#           ifdef _DEBUGPRINT_
            write(u6,*) ' sfs(I) ',sfs_check(I),' sos(I) ',sos_check(I)
#           endif
#           else
            call WarningMessage(2,'File '//trim(NAMEFILE_ANISO)//' cannot be opened. Molcas was compiled without HDF5 option.')
            call Quit_OnUserError()
#           endif
          else
            ! generating the name of the "aniso_input file for
            ! each center. Maxmimum 10 centers. CHAR(48)=0 (zero)
            if (i < 10) then
              write(namefile_aniso,'(4A)') 'aniso_',char(48+mod(i,10)),'.input'
            else if ((i >= 10) .and. (i <= 99)) then
              write(namefile_aniso,'(4A)') 'aniso_',char(48+mod(i/10,10)),char(48+mod(i,10)),'.input'
            end if
            LUANISO = Isfreeunit(20)
            call molcas_open(LUANISO,NAMEFILE_ANISO)

            if (old_aniso_format) then
              read(LUANISO,*) sfs_check(I),sos_check(I)
#             ifdef _DEBUGPRINT_
              write(u6,*) ' sfs(I) ',sfs_check(I),' sos(I) ',sos_check(I)
#             endif
            else
              call read_nss(LUANISO,sos_check(i),dbg)
              call read_nstate(LUANISO,sfs_check(i),dbg)
#             ifdef _DEBUGPRINT_
              write(u6,*) ' sfs(I) ',sfs_check(I),' sos(I) ',sos_check(I)
#             endif
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

    case ('TEXP')

      KeyTexp = .true.
      read(u5,*) nT_TEXP

      if (nT_TEXP <= 0) then
        call WarningMessage(2,'TEXP: Number of temperature points <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'TEXP:: = nT_TEXP',nT_TEXP
#     endif

      LINENR = LINENR+1

    case ('HEXP')

      KeyHexp = .true.
      read(u5,*) nTempMagn_HEXP,(TempMagn(i),i=1,nTempMagn)
      unused_var(TempMagn)
      read(u5,*) nH_HEXP

      if (nH_HEXP <= 0) then
        call WarningMessage(2,'HEXP: Number of field points <= 0! ')
        call Quit_OnUserError()
      end if

      if (nTempMagn_HEXP <= 0) then
        call WarningMessage(2,'HEXP: Number of temperature points <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'HEXP:: = nH_HEXP ',nH_HEXP
      write(u6,*) 'HEXP:: = nTempMagn_HEXP ',nTempMagn_HEXP
#     endif
      LINENR = LINENR+2

    case ('HINT')

      !KeyHINT = .true.
      read(u5,*) rdummy,rdummy,nH_HINT

      if (nH_HINT <= 0) then
        call WarningMessage(2,'HINT: Number of field points <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'HINT:: = nH_HINT ',nH_HINT
#     endif
      LINENR = LINENR+1

    case ('TINT')

      !KeyTINT = .true.
      read(u5,*) rdummy,rdummy,nT_TINT

      if (nT_TINT <= 0) then
        call WarningMessage(2,'TINT: Number of temperature points <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'TINT:: = nT_TINT ',nT_TINT
#     endif
      LINENR = LINENR+1

    case ('TMAG')

      !KeyTMAG = .true.
      read(u5,*) nTempMagn_TMAG

      if (nTempMagn_TMAG <= 0) then
        call WarningMessage(2,'TMAG: Number of temperatureMAGN points <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'TMAG:: = nTempMagn_TMAG ',nTempMagn_TMAG
#     endif
      LINENR = LINENR+1

    case ('MVEC')

      !KeyMVEC = .true.
      read(u5,*) nDir

      if (nDir <= 0) then
        call WarningMessage(2,'MVEC: Number of nDir points <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'MVEC:: = nDir ',nDir
#     endif
      LINENR = LINENR+1

    case ('ZEEM')

      !KeyZEEM = .false.
      read(u5,*) nDirZee

      if (nDirZee <= 0) then
        call WarningMessage(2,'ZEEM: Number of nDirZee points <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'ZEEM:: = nDirZee ',nDirZee
#     endif
      LINENR = LINENR+1

    case ('MLTP')

      !KeyMLTP = .true.
      read(u5,*) nMult

      if (nMult <= 0) then
        call WarningMessage(2,'MLTP: Number of multiplets <= 0! ')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'MLTP:: =nMult ',nMult
#     endif
      LINENR = LINENR+1

    case ('LIN9','LIN3','LIN1','ALIN','PAIR','ITOJ')

      KeyPair = .true.
      read(u5,*,iostat=istatus) nPair
      if (istatus /= 0) call Error(1)

      if (nPair <= 0) then
        call WarningMessage(2,'PAIR OR LINx:: Number of interacting pairs <= 0!')
        call Quit_OnUserError()
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'PAIR:: =nPair ',nPair
#     endif

      if (LINE(1:4) == 'ITOJ') then
        iline = 0
        do i=1,npair
          imaxrank(i,1) = 0
          imaxrank(i,2) = 0

          read(u5,*,iostat=istatus) idummy,idummy,imaxrank(i,1),imaxrank(i,2)
          if (istatus /= 0) call Error(1)
          do irank1=1,2*imaxrank(i,1)+1
            do irank2=1,2*imaxrank(i,2)+1
              read(u5,*,iostat=istatus) idummy,idummy,idummy,idummy,rdummy,rdummy
              if (istatus /= 0) call Error(1)
              iline = iline+1
            end do
          end do
        end do ! i
        MxRank1 = maxval(imaxrank(1:npair,1))
        MxRank2 = maxval(imaxrank(1:npair,2))
        LINENR = LINENR+iline
      end if
      LINENR = LINENR+1

    case ('COOR')

      KeyCoor = .true.
      LINENR = LINENR+1

  end select

end do

if ((.not. KeyPair) .and. KeyCoor) nPair = nTri_Elem(nCenter-1)

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
#ifdef _DEBUGPRINT_
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
#endif

return

unused_var(idummy)
unused_var(rdummy)

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

# include "warnings.h"

  select case (code)
    case (1)
      write(u6,*) ' FETCH_INIT_CONST: Error reading "poly_aniso.input" '
      write(u6,*) ' near line nr.',LINENR+1
    case (2)
      write(u6,*) ' FETCH_INIT_CONST: Unexpected End of input file.'
  end select
  call quit(_RC_INPUT_ERROR_)

end subroutine Error

end subroutine fetch_init_const
