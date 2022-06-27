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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************

subroutine RDVEC_(FName,LU_,LABEL,IUHF,NSYM,NBAS,NORB,CMO,CMO_ab,OCC,OCC_ab,EORB,EORB_ab,INDT,TITLE,iWarn,iErr,iWFtype)
!-----------------------------------------------------------------------
! Advanced RdVec (to remove all clones!)
!-----------------------------------------------------------------------
!iWFtype =  0  -- Unknown origin of orbitals
!           1  -- Orbitals for Guessorb
!           2  -- Orbitals for closed shell HF
!           3  -- Orbitals for closed shell DFT
!           4  -- Orbitals for unrestricted HF
!           5  -- Orbitals for unrestricted DFT
!           6  -- Natural orbitals for unrestricted HF
!           7  -- Natural orbitals for unrestricted DFT
!           8  --
!-----------------------------------------------------------------------

use InpOrbFmt, only: FmtEne, FmtInd, FmtOcc, FmtOrb, mxVer, Magic, nDivEne, nDivInd, nDivOcc, nDivOrb, nSkpInd
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
character(len=*), intent(in) :: FName, LABEL
integer(kind=iwp), intent(in) :: LU_, iUHF, NSYM, NBAS(NSYM), NORB(NSYM), iWarn
real(kind=wp), intent(_OUT_) :: CMO(*), CMO_ab(*), OCC(*), OCC_ab(*), EORB(*), EORB_ab(*)
integer(kind=iwp), intent(_OUT_) :: INDT(*)
character(len=*), intent(out) :: TITLE
integer(kind=iwp), intent(out) :: iErr, iWFType
integer(kind=iwp) :: i, iA, iB, IBAS, IBASEND, iBeta, iCMO, iEne, iInd, IND, iOcc, IORB, IORBEND, iShift, istatus, ISYM, iVer, &
                     jVer, KCMO, KOCC, Lu, myiUHF, myNSYM, nDiv
logical(kind=iwp) :: Exists
character(len=256) :: LINE
character(len=40) :: FRMT
character(len=10) :: Buff
integer(kind=iwp), allocatable :: myNBAS(:), myNORB(:)
character(len=*), parameter :: Crypt = 'fi123sd', &
                               CryptUP = 'FIXXXSD', &
                               Location = 'rdVec_'

Line = 'not defined yet'

! Analyze Label

iCMO = 0
iOCC = 0
iEne = 0
iInd = 0
iBeta = 0
if (index(Label,'C') /= 0) iCMO = 1
if (index(Label,'O') /= 0) iOcc = 1
if (index(Label,'E') /= 0) iEne = 1
if (index(Label,'I') /= 0) iInd = 1
if (index(Label,'A') /= 0) iBeta = -1
if (index(Label,'B') /= 0) iBeta = 1
!----------------------------------------------------------------------*
! Open file Name                                                       *
!----------------------------------------------------------------------*
iErr = 0
Lu = Lu_
call OpnFl(FName,Lu,Exists)
if (.not. Exists) then
  write(u6,*) 'RdVec: File ',trim(FName),' not found!'
  call Abend()
end if
rewind(LU)
!----------------------------------------------------------------------*
! Check version!                                                       *
!----------------------------------------------------------------------*
read(LU,'(A256)',iostat=istatus) Line
if (istatus /= 0) call Error()
iVer = 0
do jVer=1,mxVer
  if (Magic(jVer) == Line(1:len(Magic(jVer)))) iVer = jVer
end do

if (iVer == 0) then
  call SysWarnMsg(Location,'INPORB file in old format',' ')
  call Abend()
end if
!----------------------------------------------------------------------*
! INFO section, read it unconditionally                                *
!----------------------------------------------------------------------*
do
  read(LU,'(A256)',iostat=istatus) Line
  if (istatus /= 0) call Error()
  if (Line(1:5) == '#INFO') exit
end do
read(Lu,'(a)',iostat=istatus) Title
if (istatus /= 0) call Error()
read(Lu,'(a)',iostat=istatus) Line
if (istatus /= 0) call Error()
Line(76:80) = '0 0 0'
read(Line,*) myiUHF,myNSYM,iWFtype
!read(Lu,*,iostat=istatus) myiUHF,myNSYM
!if (istatus /= 0) call Error()
! In case of UHF mismatch:
if (myiUHF /= iUHF) then
  ! Stop if UHF requested, but the INPORB is not UHF
  if (myiUHF == 0) then
    call SysWarnFileMsg(Location,FName,'IUHF does not match',' ')
    call Abend()
  end if
  ! With a UHF INPORB, only go on if alpha or beta orbitals
  ! explicitly requested
  if ((iUHF == 0) .and. (iBeta == 0)) then
    call SysWarnFileMsg(Location,FName,'IUHF does not match',' ')
    call Abend()
  end if
end if
if (myNSYM /= NSYM) then
  call SysWarnFileMsg(Location,FName,'NSYM does not match',' ')
  call Abend()
end if
call mma_allocate(myNBAS,NSYM,label='myNBAS')
call mma_allocate(myNORB,NSYM,label='myNORB')
read(Lu,*,iostat=istatus) (myNBAS(i),i=1,NSYM)
if (istatus /= 0) call Error()
read(Lu,*,iostat=istatus) (myNORB(i),i=1,NSYM)
if (istatus /= 0) call Error()
!----------------------------------------------------------------------*
! Do checks                                                            *
!----------------------------------------------------------------------*
do i=1,NSYM
  if (myNBAS(i) /= NBAS(i)) then
    Line = 'NBAS does not match'
    if (iWarn == 1) then
      call SysWarnMsg(Location,Line,' ')
    else
      call SysWarnFileMsg(Location,FName,Line,' ')
      call Abend()
    end if
  end if
end do
if (iWarn > 0) then
  do i=1,NSYM
    if (myNORB(i) < NORB(i)) then
      Line = 'NORB does not match'
      if (iWarn == 1) then
        call SysWarnMsg(Location,Line,' ')
      else
        call SysWarnFileMsg(Location,FName,Line,' ')
        call Abend()
      end if
    end if
  end do
end if
call mma_deallocate(myNBAS)
!----------------------------------------------------------------------*
! ORB section                                                          *
!----------------------------------------------------------------------*
if (iCMO == 1) then
  nDiv = nDivOrb(iVer)
  FRMT = FmtOrb(iVer)
  rewind(LU)
  do
    read(LU,'(A256)',iostat=istatus) Line
    if (istatus /= 0) call Error()
    if (Line(1:4) == '#ORB') exit
  end do
  KCMO = 0
  do ISYM=1,NSYM
    do IORB=1,myNORB(ISYM)
      do IBAS=1,NBAS(ISYM),NDIV
        IBASEND = min(IBAS+NDIV-1,NBAS(ISYM))
        do
          read(LU,'(A256)',iostat=istatus) LINE
          if (istatus /= 0) call Error()
          if (LINE(1:1) /= '*') exit
        end do
        if (iOrb <= nOrb(iSym)) then
          read(LINE,FRMT,iostat=istatus) (CMO(I+KCMO),I=IBAS,IBASEND)
          if (istatus /= 0) call Error()
        end if
      end do
      if (iOrb <= nOrb(iSym)) KCMO = KCMO+NBAS(ISYM)
    end do
  end do
  if ((iUHF == 1) .or. (iBeta == 1)) then
    do
      read(LU,'(A256)',iostat=istatus) Line
      if (istatus /= 0) call Error()
      if (Line(1:5) == '#UORB') exit
    end do
    KCMO = 0
    do ISYM=1,NSYM
      do IORB=1,NORB(ISYM)
        do IBAS=1,NBAS(ISYM),NDIV
          IBASEND = min(IBAS+NDIV-1,NBAS(ISYM))
          do
            read(LU,'(A256)',iostat=istatus) LINE
            if (istatus /= 0) call Error()
            if (LINE(1:1) /= '*') exit
          end do
          if (iOrb <= nOrb(iSym)) then
            if (iBeta == 1) then
              read(LINE,FRMT,iostat=istatus) (CMO(I+KCMO),I=IBAS,IBASEND)
              if (istatus /= 0) call Error()
            else
              read(LINE,FRMT,iostat=istatus) (CMO_ab(I+KCMO),I=IBAS,IBASEND)
              if (istatus /= 0) call Error()
            end if
          end if
        end do
        if (iOrb <= nOrb(iSym)) KCMO = KCMO+NBAS(ISYM)
      end do
    end do
  end if ! iUHF
end if ! iCMO
!----------------------------------------------------------------------*
! OCC section                                                          *
!----------------------------------------------------------------------*
if (iOcc == 1) then
  nDiv = nDivOcc(iVer)
  FRMT = FmtOcc(iVer)
  rewind(LU)
  do
    read(LU,'(A256)',iostat=istatus) Line
    if (istatus /= 0) call Error()
    if (Line(1:4) == '#OCC') exit
  end do
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,myNORB(ISYM),NDIV
      IORBEND = min(IORB+NDIV-1,myNORB(ISYM))
      do
        read(LU,'(A256)',iostat=istatus) LINE
        if (istatus /= 0) call Error()
        if (LINE(1:1) /= '*') exit
      end do
      read(LINE,FRMT,iostat=istatus) (OCC(I+KOCC),I=IORB,IORBEND)
      if (istatus /= 0) call Error()
    end do
    !KOCC = KOCC+myNORB(ISYM)
    KOCC = KOCC+nOrb(iSym)
  end do
  if ((iUHF == 1) .or. (iBeta == 1)) then
    do
      read(LU,'(A256)',iostat=istatus) Line
      if (istatus /= 0) call Error()
      if (Line(1:5) == '#UOCC') exit
    end do
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,myNORB(ISYM),NDIV
        IORBEND = min(IORB+NDIV-1,myNORB(ISYM))
        do
          read(LU,'(A256)',iostat=istatus) LINE
          if (istatus /= 0) call Error()
          if (LINE(1:1) /= '*') exit
        end do
        if (iBeta == 1) then
          read(LINE,FRMT,iostat=istatus) (OCC(I+KOCC),I=IORB,IORBEND)
          if (istatus /= 0) call Error()
        else
          read(LINE,FRMT,iostat=istatus) (OCC_ab(I+KOCC),I=IORB,IORBEND)
          if (istatus /= 0) call Error()
        end if
      end do
      !KOCC = KOCC+myNORB(ISYM)
      KOCC = KOCC+nOrb(iSym)
    end do
  end if ! iUHF
end if ! iOCC
!----------------------------------------------------------------------*
! ONE section                                                          *
!----------------------------------------------------------------------*
if (iEne == 1) then
  nDiv = nDivEne(iVer)
  FRMT = FmtEne(iVer)
  rewind(LU)
  do
    read(LU,'(A256)',iostat=istatus) Line
    if (istatus /= 0) then
      call End2()
      return
    end if
    if (Line(1:4) == '#ONE') exit
  end do
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,myNORB(ISYM),NDIV
      IORBEND = min(IORB+NDIV-1,myNORB(ISYM))
      do
        read(LU,'(A256)',iostat=istatus) LINE
        if (istatus /= 0) call Error()
        if (LINE(1:1) /= '*') exit
      end do
      read(LINE,FRMT,iostat=istatus) (EORB(I+KOCC),I=IORB,IORBEND)
      if (istatus /= 0) call Error()
    end do
    !KOCC = KOCC+myNORB(ISYM)
    KOCC = KOCC+nOrb(iSym)
  end do
  if ((iUHF == 1) .or. (iBeta == 1)) then
    do
      read(LU,'(A256)',iostat=istatus) Line
      if (istatus /= 0) call Error()
      if (Line(1:5) == '#UONE') exit
    end do
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,myNORB(ISYM),NDIV
        IORBEND = min(IORB+NDIV-1,myNORB(ISYM))
        do
          read(LU,'(A256)',iostat=istatus) LINE
          if (istatus /= 0) call Error()
          if (LINE(1:1) /= '*') exit
        end do
        if (iBeta == 1) then
          read(LINE,FRMT,iostat=istatus) (EORB(I+KOCC),I=IORB,IORBEND)
          if (istatus /= 0) call Error()
        else
          read(LINE,FRMT,iostat=istatus) (EORB_ab(I+KOCC),I=IORB,IORBEND)
          if (istatus /= 0) call Error()
        end if
      end do
      !KOCC = KOCC+myNORB(ISYM)
      KOCC = KOCC+nOrb(iSym)
    end do
  end if ! iUHF
end if ! iOne
!----------------------------------------------------------------------*
! INDEX section                                                        *
!----------------------------------------------------------------------*
if (iInd == 1) then
  rewind(LU)
  do
    read(LU,'(A256)',iostat=istatus) Line
    if (istatus /= 0) then
      call End2()
      return
    end if
    if (Line(1:6) == '#INDEX') exit
  end do
  FRMT = FMTIND(iVer)
  nDiv = nDivInd(iVer)
  iShift = 1
  do ISYM=1,NSYM
    !iShift = (ISYM-1)*7
    do i=1,nSkpInd(iVer)
      read(LU,*)
    end do
    do IORB=1,myNORB(ISYM),nDiv
      read(LU,FRMT,iostat=istatus) Buff
      if (istatus /= 0) then
        call End2()
        return
      end if
      do i=1,nDiv
        IND = index(Crypt,Buff(i:i))+index(CryptUp,Buff(i:i))
        if (Buff(i:i) /= ' ') then
          if (IND == 0) then
            write(u6,*) '* ERROR IN RDVEC WHILE READING TypeIndex'
            write(u6,'(3A)') '* Type=',Buff(i:i),' is UNKNOWN'
            call End2()
            return
          end if
          IndT(iShift) = IND
          iShift = iShift+1
        end if
      end do
    end do
  end do

  iA = 1
  iB = 1
  do iSym=1,nSym
    IndT(iB:iB+nOrb(iSym)-1) = IndT(iA:iA+nOrb(iSym)-1)
    iA = iA+nBas(iSym)
    iB = iB+nOrb(iSym)
  end do
end if  ! Index
close(Lu)
call End1()

return

contains

subroutine End1()
  call mma_deallocate(myNORB)
end subroutine End1

subroutine End2()
  !--------------------------------------------------------------------*
  ! a special case - INDEX information is not found                    *
  !--------------------------------------------------------------------*
  iErr = 1
  write(u6,*) '* TypeIndex information is IGNORED *'
  close(Lu)
  call End1()
end subroutine End2

subroutine Error()
  call SysWarnFileMsg(Location,FName,'Error during reading INPORB\n',Line)
  call Abend()
end subroutine Error

end subroutine RDVEC_
