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

subroutine RDVEC_(Name,LU_,LABEL,IUHF,NSYM,NBAS,NORB,CMO,CMO_ab,OCC,OCC_ab,EORB,EORB_ab,INDT,TITLE,iWarn,iErr,iWFtype)
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

implicit real*8(A-H,O-Z)
#include "WrkSpc.fh"
dimension NBAS(NSYM), NORB(NSYM), CMO(*), OCC(*), INDT(*), EORB(*)
dimension CMO_ab(*), OCC_ab(*), EORB_ab(*)
character*(*) TITLE, Name, Label
character LINE*256, FMT*40
logical Exist
character*7 Crypt, CryptUp
character*10 Buff
! Note! the size of Magic must be exact (thanks to MS formatted inporb!)
character*8 Location
data Crypt/'fi123sd'/
data CryptUP/'FIXXXSD'/
#include "inporbfmt.fh"

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
call OpnFl(Name,Lu,Exist)
if (.not. Exist) then
  write(6,*) 'RdVec: File ',Name(1:index(Name,' ')),' not found!'
  call Abend()
end if
rewind(LU)
!----------------------------------------------------------------------*
! Check version!                                                       *
!----------------------------------------------------------------------*
read(LU,'(A256)',end=999,ERR=999) Line
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
50 continue
read(LU,'(A256)',end=999,ERR=999) Line
if (Line(1:5) /= '#INFO') goto 50
read(Lu,'(a)',end=999,err=999) Title
read(Lu,'(a)',end=999,Err=999) Line
Line(76:80) = '0 0 0'
read(Line,*) myiUHF,myNSYM,iWFtype
!read(Lu,*,end=999,err=999) myiUHF,myNSYM
! In case of UHF mismatch:
if (myiUHF /= iUHF) then
  ! Stop if UHF requested, but the INPORB is not UHF
  if (myiUHF == 0) then
    call SysWarnFileMsg(Location,Name,'IUHF does not match',' ')
    call Abend()
  end if
  ! With a UHF INPORB, only go on if alpha or beta orbitals
  ! explicitly requested
  if ((iUHF == 0) .and. (iBeta == 0)) then
    call SysWarnFileMsg(Location,Name,'IUHF does not match',' ')
    call Abend()
  end if
end if
if (myNSYM /= NSYM) then
  call SysWarnFileMsg(Location,Name,'NSYM does not match',' ')
  call Abend()
end if
call GetMem('MYNBAS','Allo','Inte',imyNBAS,NSYM)
call GetMem('MYNORB','Allo','Inte',imyNORB,NSYM)
read(Lu,*,end=999,err=999) (iWork(imyNBAS+i-1),i=1,NSYM)
read(Lu,*,end=999,err=999) (iWork(imyNORB+i-1),i=1,NSYM)
!----------------------------------------------------------------------*
! Do checks                                                            *
!----------------------------------------------------------------------*
do i=1,NSYM
  if (iWork(imyNBAS+i-1) /= NBAS(i)) then
    Line = 'NBAS does not match'
    if (iWarn == 1) then
      call SysWarnMsg(Location,Line,' ')
    else
      call SysWarnFileMsg(Location,Name,Line,' ')
      call Abend()
    end if
  end if
end do
if (iWarn > 0) then
  do i=1,NSYM
    if (iWork(imyNORB+i-1) < NORB(i)) then
      Line = 'NORB does not match'
      if (iWarn == 1) then
        call SysWarnMsg(Location,Line,' ')
      else
        call SysWarnFileMsg(Location,Name,Line,' ')
        call Abend()
      end if
    end if
  end do
end if
call GetMem('MYNBAS','Free','Inte',imyNBAS,NSYM)
!----------------------------------------------------------------------*
! ORB section                                                          *
!----------------------------------------------------------------------*
if (iCMO == 1) then
  nDiv = nDivOrb(iVer)
  FMT = FmtOrb(iVer)
  rewind(LU)
51 continue
  read(LU,'(A256)',end=999,ERR=999) Line
  if (Line(1:4) /= '#ORB') goto 51
  KCMO = 0
  do ISYM=1,NSYM
    do IORB=1,iWork(imyNORB+ISYM-1)
      do IBAS=1,NBAS(ISYM),NDIV
        IBASEND = min(IBAS+NDIV-1,NBAS(ISYM))
111     continue
        read(LU,'(A256)',end=999,ERR=999) LINE
        if (LINE(1:1) == '*') goto 111
        if (iOrb <= nOrb(iSym)) then
          read(LINE,FMT,err=888,end=888) (CMO(I+KCMO),I=IBAS,IBASEND)
        end if
      end do
      if (iOrb <= nOrb(iSym)) KCMO = KCMO+NBAS(ISYM)
    end do
  end do
  if ((iUHF == 1) .or. (iBeta == 1)) then
52  continue
    read(LU,'(A256)',end=999,ERR=999) Line
    if (Line(1:5) /= '#UORB') goto 52
    KCMO = 0
    do ISYM=1,NSYM
      do IORB=1,iWork(imyNORB+ISYM-1)
        do IBAS=1,NBAS(ISYM),NDIV
          IBASEND = min(IBAS+NDIV-1,NBAS(ISYM))
112       continue
          read(LU,'(A256)',end=999,ERR=999) LINE
          if (LINE(1:1) == '*') goto 112
          if (iOrb <= nOrb(iSym)) then
            if (iBeta == 1) then
              read(LINE,FMT,err=888,end=888) (CMO(I+KCMO),I=IBAS,IBASEND)
            else
              read(LINE,FMT,err=888,end=888) (CMO_ab(I+KCMO),I=IBAS,IBASEND)
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
  FMT = FmtOcc(iVer)
  rewind(LU)
53 continue
  read(LU,'(A256)',end=999,ERR=999) Line
  if (Line(1:4) /= '#OCC') goto 53
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,iWork(imyNORB+ISYM-1),NDIV
      IORBEND = min(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
113   continue
      read(LU,'(A256)',end=999,ERR=999) LINE
      if (LINE(1:1) == '*') goto 113
      read(LINE,FMT,err=888,end=888) (OCC(I+KOCC),I=IORB,IORBEND)
    end do
    !KOCC = KOCC+iWork(imyNORB+ISYM-1)
    KOCC = KOCC+nOrb(iSym)
  end do
  if ((iUHF == 1) .or. (iBeta == 1)) then
54  continue
    read(LU,'(A256)',end=999,ERR=999) Line
    if (Line(1:5) /= '#UOCC') goto 54
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,iWork(imyNORB+ISYM-1),NDIV
        IORBEND = min(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
114     continue
        read(LU,'(A256)',end=999,ERR=999) LINE
        if (LINE(1:1) == '*') goto 114
        if (iBeta == 1) then
          read(LINE,FMT,err=888,end=888) (OCC(I+KOCC),I=IORB,IORBEND)
        else
          read(LINE,FMT,err=888,end=888) (OCC_ab(I+KOCC),I=IORB,IORBEND)
        end if
      end do
      !KOCC = KOCC+iWork(imyNORB+ISYM-1)
      KOCC = KOCC+nOrb(iSym)
    end do
  end if ! iUHF
end if ! iOCC
!----------------------------------------------------------------------*
! ONE section                                                          *
!----------------------------------------------------------------------*
if (iEne == 1) then
  nDiv = nDivEne(iVer)
  FMT = FmtEne(iVer)
  rewind(LU)
55 continue
  read(LU,'(A256)',end=666,ERR=666) Line
  if (Line(1:4) /= '#ONE') goto 55
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,iWork(imyNORB+ISYM-1),NDIV
      IORBEND = min(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
115   continue
      read(LU,'(A256)',end=999,ERR=999) LINE
      if (LINE(1:1) == '*') goto 115
      read(LINE,FMT,err=888,end=888) (EORB(I+KOCC),I=IORB,IORBEND)
    end do
    !KOCC = KOCC+iWork(imyNORB+ISYM-1)
    KOCC = KOCC+nOrb(iSym)
  end do
  if ((iUHF == 1) .or. (iBeta == 1)) then
56  continue
    read(LU,'(A256)',end=999,ERR=999) Line
    if (Line(1:5) /= '#UONE') goto 56
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,iWork(imyNORB+ISYM-1),NDIV
        IORBEND = min(IORB+NDIV-1,iWork(imyNORB+ISYM-1))
116     continue
        read(LU,'(A256)',end=999,ERR=999) LINE
        if (LINE(1:1) == '*') goto 116
        if (iBeta == 1) then
          read(LINE,FMT,err=888,end=888) (EORB(I+KOCC),I=IORB,IORBEND)
        else
          read(LINE,FMT,err=888,end=888) (EORB_ab(I+KOCC),I=IORB,IORBEND)
        end if
      end do
      !KOCC = KOCC+iWork(imyNORB+ISYM-1)
      KOCC = KOCC+nOrb(iSym)
    end do
  end if ! iUHF
end if ! iOne
!----------------------------------------------------------------------*
! INDEX section                                                        *
!----------------------------------------------------------------------*
if (iInd == 1) then
  rewind(LU)
57 continue
  read(LU,'(A256)',end=666,ERR=666) Line
  if (Line(1:6) /= '#INDEX') goto 57
  FMT = FMTIND(iVer)
  nDiv = nDivInd(iVer)
  iShift = 1
  do ISYM=1,NSYM
    !iShift = (ISYM-1)*7
    do i=1,nSkpInd(iVer)
      read(LU,*)
    end do
    do IORB=1,iWork(imyNORB+ISYM-1),nDiv
      read(LU,FMT,err=666,end=666) Buff
      do i=1,nDiv
        IND = index(Crypt,Buff(i:i))+index(CryptUp,Buff(i:i))
        if (Buff(i:i) /= ' ') then
          if (IND == 0) then
            write(6,*) '* ERROR IN RDVEC WHILE READING TypeIndex'
            write(6,'(3A)') '* Type=',Buff(i:i),' is UNKNOWN'
            write(6,*) '* TypeIndex information is IGNORED'
            iErr = 1
            close(Lu)
            goto 777
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
    call iCopy(nOrb(iSym),IndT(iA),1,IndT(iB),1)
    iA = iA+nBas(iSym)
    iB = iB+nOrb(iSym)
  end do
end if  ! Index
close(Lu)
goto 777
!----------------------------------------------------------------------*
! a special case - INDEX information is not found                      *
!----------------------------------------------------------------------*
666 continue
iErr = 1
write(6,*) '* TypeIndex information is IGNORED *'
close(Lu)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
777 continue
call GetMem('MYNORB','Free','Inte',imyNORB,NSYM)

return
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
999 continue
call SysWarnFileMsg(Location,Name,'Error during reading INPORB\n',Line)
call Abend()
888 continue
call SysWarnFileMsg(Location,Name,'Error during reading INPORB\n',Line)
call Abend()

end subroutine RDVEC_
