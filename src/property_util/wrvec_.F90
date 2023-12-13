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

subroutine WRVEC_(FName,LU_,LABEL,IUHF,NSYM,NBAS,NORB,CMO,CMO_ab,OCC,OCC_ab,EORB,EORB_ab,INDT,TITLE,iWFtype)
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

use InpOrbFmt, only: FmtEne, FmtInd, FmtOcc, FmtOccHR, FmtOrb, iVer22, Magic, mxVer, nDivEne, nDivInd, nDivOcc, nDivOccHR, &
                     nDivOrb, nSkpInd
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: FName, LABEL
integer(kind=iwp), intent(in) :: LU_, IUHF, NSYM, NBAS(NSYM), NORB(NSYM), INDT(*), iWFtype
real(kind=wp), intent(in) :: CMO(*), CMO_ab(*), OCC(*), OCC_ab(*), EORB(*), EORB_ab(*)
character(len=*), intent(in) :: TITLE
integer(kind=iwp) :: I, iAppend, IBAS, IBASEND, iBasis, iBB, iBuff(0:7), iCMO, iDefault, iEne, iExtras, iInd, iKoor, iLab, iOCC, &
                     IORB, IORBEND, Ip, iShift, isIndex, istatus, ISYM, iTwoE, iVer = 0, jVer, KCMO, KOCC, lin, Lu, NDIV, nDNA
logical(kind=iwp) :: Exists, IsBorn
!-SVC: variable to hold birth certificate
character(len=256) :: cDNA
character(len=120) :: str
character(len=84) :: FRMT
character(len=20) :: InpOrbVer
character(len=12) :: lBuff
character(len=10) :: Buff
character(len=8) :: Line
character(len=*), parameter :: Crypt = 'fi123sd'
integer(kind=iwp) :: isFreeUnit

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
call OpnFl(FName,Lu,Exists)
rewind(Lu)
if (iAppend == 1) then
  iInd = 1
  do
    read(Lu,'(A)',iostat=istatus) Line
    if (istatus /= 0) exit
    if (Line == '#INDEX') then
      isIndex = 1
      call WrInd()
      return
    end if
  end do
  call Append_file(LU)
!# ifdef NAGFOR
!  backspace(LU)
!# endif
  call WrInd()
  return
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
if (TITLE(1:1) /= '*') then
  write(LU,'(A)') '*'//trim(TITLE)
else
  write(LU,'(A)') trim(TITLE)
end if
write(LU,'(3i8)') IUHF,NSYM,iWFtype
write(LU,'(8i8)') (NBAS(I),I=1,NSYM)
write(LU,'(8i8)') (NORB(I),I=1,NSYM)
call qpg_cArray('BirthCertificate',IsBorn,nDNA)
if (.not. IsBorn) then
  write(u6,*) 'RunFile has no Birth Certificate'
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
    write(LU,'(ES19.12)') EORB_ab(1)
  end if
end if

! ORB section
if (iCMO == 1) then
  NDIV = nDivOrb(iVer)
  FRMT = FmtOrb(iVer)
  write(Lu,'(A)') '#ORB'
  KCMO = 0
  do ISYM=1,NSYM
    do IORB=1,NORB(ISYM)
      write(LU,'(A,2I5)') '* ORBITAL',ISYM,IORB
      do IBAS=1,NBAS(ISYM),NDIV
        IBASEND = min(IBAS+NDIV-1,NBAS(ISYM))
        write(LU,FRMT) (CMO(I+KCMO),I=IBAS,IBASEND)
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
          write(LU,FRMT) (CMO_ab(I+KCMO),I=IBAS,IBASEND)
        end do
        KCMO = KCMO+NBAS(ISYM)
      end do
    end do
  end if  ! UHF
end if  ! iCMO

! OCC section
if (iOcc == 1) then
  NDIV = nDivOcc(iVer)
  FRMT = FmtOcc(iVer)
  write(Lu,'(A)') '#OCC'
  write(LU,'(A)') '* OCCUPATION NUMBERS'
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,NORB(ISYM),NDIV
      IORBEND = min(IORB+NDIV-1,NORB(ISYM))
      write(LU,FRMT) (OCC(I+KOCC),I=IORB,IORBEND)
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
        write(LU,FRMT) (OCC_ab(I+KOCC),I=IORB,IORBEND)
      end do
      KOCC = KOCC+NORB(ISYM)
    end do
  end if  ! UHF

  NDIV = nDivOccHR(iVer)
  if (NDIV > 0) then
    FRMT = FmtOccHR(iVer)
    write(Lu,'(A)') '#OCHR'
    write(LU,'(A)') '* OCCUPATION NUMBERS (HUMAN-READABLE)'
    KOCC = 0
    do ISYM=1,NSYM
      do IORB=1,NORB(ISYM),NDIV
        IORBEND = min(IORB+NDIV-1,NORB(ISYM))
        write(LU,FRMT) (OCC(I+KOCC),I=IORB,IORBEND)
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
          write(LU,FRMT) (OCC_ab(I+KOCC),I=IORB,IORBEND)
        end do
        KOCC = KOCC+NORB(ISYM)
      end do
    end if  ! UHF
  end if
end if  ! iOcc

! ONE section
if (iEne == 1) then
  NDIV = nDivEne(iVer)
  FRMT = FmtEne(iVer)
  write(Lu,'(A)') '#ONE'
  write(LU,'(A)') '* ONE ELECTRON ENERGIES'
  KOCC = 0
  do ISYM=1,NSYM
    do IORB=1,NORB(ISYM),NDIV
      IORBEND = min(IORB+NDIV-1,NORB(ISYM))
      write(LU,FRMT) (EORB(I+KOCC),I=IORB,IORBEND)
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
        write(LU,FRMT) (EORB_ab(I+KOCC),I=IORB,IORBEND)
      end do
      KOCC = KOCC+NORB(ISYM)
    end do
  end if  ! UHF
end if  ! iEne
call getenvf('MOLCAS_SAGIT',str)
if ((str(1:1) == 'y') .or. (str(1:1) == 'Y')) then
  iKoor = iCMO
  iBasis = iCMO
else
  iKoor = 0
  iBasis = 0
end if
if ((iKoor == 1) .and. (iBasis == 1)) then
  lin = 16
  lin = isfreeunit(lin)
  call molcas_open(lin,'ORB.std')
  do
    read(lin,'(a)',iostat=istatus) str
    if (istatus /= 0) exit
    write(Lu,'(a)') trim(str)
  end do
  close(lin)
end if

call WrInd()

return

contains

subroutine WrInd()
  integer(kind=iwp) :: I, IB, IORB, ISYM

  ! INDEX section. NOTE THIS SECTION SHOULD ALWAYS BE LAST (Gv constraint)
  if (iInd == 1) then
    if ((iAppend == 0) .or. ((iAppend == 1) .and. (isIndex == 0))) then
      write(Lu,'(A)') '#INDEX'
    end if
    iShift = 0
    nDiv = nDivInd(iVer)

    iBuff(0) = 1
    !do i=1,7
    i = 1
    do
      iBuff(i) = iBuff(i-1)+IndT(i+iShift)
      i = i+1
      if (i > 7) exit
    end do
    !end do

    do ISYM=1,NSYM
      if (nSkpInd(iVer) > 0) write(Lu,'(A)') '* 1234567890'
      iLab = 0
      iBuff(0) = 1
      !do i=1,7
      i = 1
      do
        iBuff(i) = iBuff(i-1)+IndT(i+iShift)
        i = i+1
        if (i > 7) exit
      end do
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

end subroutine WrInd

end subroutine WRVEC_
