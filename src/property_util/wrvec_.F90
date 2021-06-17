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

subroutine WRVEC_(Name,LU_,LABEL,IUHF,NSYM,NBAS,NORB,CMO,CMO_ab,OCC,OCC_ab,EORB,EORB_ab,INDT,TITLE,iWFtype)
! The routine to dump information to INPORB
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
dimension NBAS(NSYM), NORB(NSYM), CMO(*), OCC(*), EORB(*), INDT(*)
dimension CMO_ab(*), OCC_ab(*), EORB_ab(*)
character*(*) TITLE, Name, LABEL
character*8 Line
character FMT*40
logical Exist
character*10 Buff
character*12 lBuff
dimension iBuff(0:7)
character*7 Crypt
!-SVC: variable to hold birth certificate
character cDNA*256
character*120 inout
character*20 InpOrbVer
logical IsBorn
data Crypt/'fi123sd'/
integer, save :: iVer = 0
#include "inporbfmt.fh"

! Analyze Label
iCMO = 0
iOCC = 0
iEne = 0
iInd = 0
iTwoE = 0
iAppend = 0
iExtras = 0
iKoor = 0
iBasis = 0
if (index(Label,'C') /= 0) iCMO = 1
if (index(Label,'O') /= 0) iOcc = 1
if (index(Label,'E') /= 0) iEne = 1
if (index(Label,'I') /= 0) iInd = 1
if (index(Label,'T') /= 0) iTwoE = 1
if (index(Label,'A') /= 0) iAppend = 1
if (index(Label,'X') /= 0) iExtras = 1
if (index(Label,'K') /= 0) iKoor = 1
if (index(Label,'B') /= 0) iBasis = 1
isIndex = 0
Lu = Lu_
call OpnFl(Name,Lu,Exist)
rewind(Lu)
if (iAppend == 1) then
  iInd = 1
50 continue
  read(Lu,'(A)',end=102,err=102) Line
  if (Line == '#INDEX') then
    isIndex = 1
    goto 100
  end if
  goto 50
102 continue
  call Append_file(LU)
!# ifdef NAGFOR
!  backspace(LU)
!# endif
  goto 100
end if

! Get version
iDefault = iVer22
if (iVer == 0) then
  call getenvf('MOLCAS_INPORB_VERSION',InpOrbVer)
  if (InpOrbVer == '') then
    iVer = iDefault
  else
    InpOrbVer = '#INPORB '//trim(adjustl(InpOrbVer))
    do jVer=1,mxVer
      if (Magic(jVer) == InpOrbVer) iVer = jVer
    end do
    if (iVer == 0) then
      call WarningMessage(0,'Unknown INPORB version, using the default')
      iVer = iDefault
    end if
  end if
end if

! Write INFO header
write(LU,'(A)') Magic(iVer)
write(Lu,'(A)') '#INFO'

KCMO = 0
if (TITLE(1:1) /= '*') TITLE = '*'//TITLE(:len(TITLE)-1)
write(LU,'(A)') trim(TITLE)
write(LU,'(3i8)') IUHF,NSYM,iWFtype
write(LU,'(8i8)') (NBAS(I),I=1,NSYM)
write(LU,'(8i8)') (NORB(I),I=1,NSYM)
call qpg_cArray('BirthCertificate',IsBorn,nDNA)
if (.not. IsBorn) then
  write(6,*) 'RunFile has no Birth Certificate'
else
  cDNA = ' '
  call Get_cArray('BirthCertificate',cDNA(:nDNA),nDNA)
  write(LU,'(A)') '*BC:'//trim(cDNA)
end if

! Extras section
if (iTwoE == 1) iExtras = 1  ! so far only case

if (iExtras == 1) then
  write(Lu,'(A)') '#EXTRAS'
  if (iTwoE == 1) then
    write(LU,'(A)') '* ACTIVE TWO-EL ENERGY'
    write(LU,'(E19.12)') EORB_ab(1)
  end if
end if

! ORB section
if (iCMO == 1) then
  NDIV = nDivOrb(iVer)
  FMT = FmtOrb(iVer)
  write(Lu,'(A)') '#ORB'
  KCMO = 0
  do ISYM=1,NSYM
    do IORB=1,NORB(ISYM)
      write(LU,'(A,2I5)') '* ORBITAL',ISYM,IORB
      do IBAS=1,NBAS(ISYM),NDIV
        IBASEND = min(IBAS+NDIV-1,NBAS(ISYM))
        write(LU,FMT) (CMO(I+KCMO),I=IBAS,IBASEND)
      end do
      KCMO = KCMO+NBAS(ISYM)
    end do
  end do
  if (iUHF == 1) then
    write(Lu,'(A)') '#UORB'
    KCMO = 0
    do ISYM=1,NSYM
      do IORB=1,NORB(ISYM)
        write(LU,'(A,2I5)') '* ORBITAL',ISYM,IORB
        do IBAS=1,NBAS(ISYM),NDIV
          IBASEND = min(IBAS+NDIV-1,NBAS(ISYM))
          write(LU,FMT) (CMO_ab(I+KCMO),I=IBAS,IBASEND)
        end do
        KCMO = KCMO+NBAS(ISYM)
      end do
    end do
  end if  ! UHF
end if  ! iCMO

! OCC section
if (iOcc == 1) then
  NDIV = nDivOcc(iVer)
  FMT = FmtOcc(iVer)
  write(Lu,'(A)') '#OCC'
  write(LU,'(A)') '* OCCUPATION NUMBERS'
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,NORB(ISYM),NDIV
      IORBEND = min(IORB+NDIV-1,NORB(ISYM))
      write(LU,FMT) (OCC(I+KOCC),I=IORB,IORBEND)
    end do
    KOCC = KOCC+NORB(ISYM)
  end do

  if (iUHF == 1) then
    write(Lu,'(A)') '#UOCC'
    write(LU,'(A)') '* Beta OCCUPATION NUMBERS'
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,NORB(ISYM),NDIV
        IORBEND = min(IORB+NDIV-1,NORB(ISYM))
        write(LU,FMT) (OCC_ab(I+KOCC),I=IORB,IORBEND)
      end do
      KOCC = KOCC+NORB(ISYM)
    end do
  end if  ! UHF

  NDIV = nDivOccHR(iVer)
  if (NDIV > 0) then
    FMT = FmtOccHR(iVer)
    write(Lu,'(A)') '#OCHR'
    write(LU,'(A)') '* OCCUPATION NUMBERS (HUMAN-READABLE)'
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,NORB(ISYM),NDIV
        IORBEND = min(IORB+NDIV-1,NORB(ISYM))
        write(LU,FMT) (OCC(I+KOCC),I=IORB,IORBEND)
      end do
      KOCC = KOCC+NORB(ISYM)
    end do

    if (iUHF == 1) then
      write(Lu,'(A)') '#UOCHR'
      write(LU,'(A)') '* Beta OCCUPATION NUMBERS (HUMAN-READABLE)'
      KOCC = 0
      do ISYM=1,NSYM
        do IORB=1,NORB(ISYM),NDIV
          IORBEND = min(IORB+NDIV-1,NORB(ISYM))
          write(LU,FMT) (OCC_ab(I+KOCC),I=IORB,IORBEND)
        end do
        KOCC = KOCC+NORB(ISYM)
      end do
    end if  ! UHF
  end if
end if  ! iOcc

! ONE section
if (iEne == 1) then
  NDIV = nDivEne(iVer)
  FMT = FmtEne(iVer)
  write(Lu,'(A)') '#ONE'
  write(LU,'(A)') '* ONE ELECTRON ENERGIES'
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,NORB(ISYM),NDIV
      IORBEND = min(IORB+NDIV-1,NORB(ISYM))
      write(LU,FMT) (EORB(I+KOCC),I=IORB,IORBEND)
    end do
    KOCC = KOCC+NORB(ISYM)
  end do

  if (iUHF == 1) then
    write(Lu,'(A)') '#UONE'
    write(LU,'(A)') '* Beta ONE ELECTRON ENERGIES'
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,NORB(ISYM),NDIV
        IORBEND = min(IORB+NDIV-1,NORB(ISYM))
        write(LU,FMT) (EORB_ab(I+KOCC),I=IORB,IORBEND)
      end do
      KOCC = KOCC+NORB(ISYM)
    end do
  end if  ! UHF
end if  ! iEne
call getenvf('MOLCAS_SAGIT',inout)
if ((inout(1:1) == 'y') .or. (inout(1:1) == 'Y')) then
  iKoor = iCMO
  iBasis = iCMO
else
  iKoor = 0
  iBasis = 0
end if
if ((iKoor == 1) .and. (iBasis == 1)) then
  in = 16
  in = isfreeunit(in)
  call molcas_open(in,'ORB.std')
765 continue
  read(in,'(a)',end=766,err=766) InOut
  write(Lu,'(a)') trim(InOut)
  goto 765
766 continue
  close(in)
end if
! INDEX section. NOTE THIS SECTION SHOULD ALWAYS BE LAST (Gv constraint)
100 continue
if (iInd == 1) then
  if ((iAppend == 0) .or. ((iAppend == 1) .and. (isIndex == 0))) then
    write(Lu,'(A)') '#INDEX'
  end if
  iShift = 0
  nDiv = nDivInd(iVer)

  iBuff(0) = 1
  !do i=1,7
  i = 1
601 continue
  iBuff(i) = iBuff(i-1)+IndT(i+iShift)
  i = i+1
  if (i <= 7) goto 601

  !end do

  do ISYM=1,NSYM
    if (nSkpInd(iVer) > 0) write(Lu,'(A)') '* 1234567890'
    iLab = 0
    iBuff(0) = 1
    !do i=1,7
    i = 1
600 continue
    iBuff(i) = iBuff(i-1)+IndT(i+iShift)
    i = i+1
    if (i <= 7) goto 600

    !end do
    Ip = 1
    do IORB=1,NORB(ISYM),nDiv
      Buff = '          '
      do i=1,nDiv
        iBB = 1
        do iB=1,7
          if (Ip >= iBuff(iB)) iBB = iBB+1
        end do
        if (iBB == 8) then
          Buff(i:i) = ' '
        else
          Buff(i:i) = Crypt(iBB:iBB)
        end if
        Ip = Ip+1
      end do
      write(lBuff,FMTIND(iVer)) Buff
      if (index(FMTIND(iVer),'X') > 0) write(lBuff(1:1),'(i1)') iLab
      write(LU,'(A)') trim(lBuff)
      iLab = iLab+1
      if (iLab > 9) iLab = 0
    end do
    iShift = iShift+7
  end do

end if  ! iInd

close(Lu)

return

end subroutine WRVEC_
