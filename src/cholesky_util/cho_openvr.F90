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

subroutine CHO_OPENVR(IOPT,ID)
!
! Purpose: open (IOPT=1) or close (IOPT=2) files for vector
!          and reduced set storage as well as restart files.
!          ID=1: open local files (for parallel run)
!          ID=2: open global files (as in serial run)

use Cholesky, only: Cho_AdrVec, LuCho, LuMap, LuPri, LuRed, LuRst, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IOPT, ID
integer(kind=iwp) :: ISYM
character(len=6) :: FMAP, FNVEC(8), FRST
character(len=5) :: FNRED
character(len=*), parameter :: SECNAM = 'CHO_OPENVR'

if (IOPT == 1) then
  FMAP = 'CHOMAP'
  if (ID == 1) then
    FNRED = 'CHRDL'
    do ISYM=1,NSYM
      write(FNVEC(ISYM),'(A5,I1)') 'CHVCL',ISYM
    end do
    FRST = 'CHRSTL'
  else
    FNRED = 'CHRED'
    do ISYM=1,NSYM
      write(FNVEC(ISYM),'(A5,I1)') 'CHVEC',ISYM
    end do
    FRST = 'CHORST'
  end if
  LURED = 7
  call DANAME_MF_WA(LURED,FNRED)
  if (CHO_ADRVEC == 1) then
    do ISYM=1,NSYM
      LUCHO(ISYM) = 7
      call DANAME_MF_WA(LUCHO(ISYM),FNVEC(ISYM))
    end do
  else if (CHO_ADRVEC == 2) then
    do ISYM=1,NSYM
      LUCHO(ISYM) = 7
      call DANAME_MF(LUCHO(ISYM),FNVEC(ISYM))
    end do
  else
    call CHO_QUIT('CHO_ADRVEC out of bounds in '//SECNAM//'. Perhaps the NOCHO keyword is needed?',102)
  end if
  LURST = 7
  call DANAME_MF_WA(LURST,FRST)
  LUMAP = 7
  call DANAME(LUMAP,FMAP)
else if (IOPT == 2) then
  if (LURED > 0) then
    call DACLOS(LURED)
    LURED = 0
  end if
  do ISYM=1,NSYM
    if (LUCHO(ISYM) > 0) then
      call DACLOS(LUCHO(ISYM))
      LUCHO(ISYM) = 0
    end if
  end do
  if (LURST > 0) then
    call DACLOS(LURST)
    LURST = 0
  end if
  if (LUMAP > 0) then
    call DACLOS(LUMAP)
    LUMAP = 0
  end if
else
  write(LUPRI,*) SECNAM,': IOPT out of bounds: ',IOPT
  call CHO_QUIT('Error in '//SECNAM,104)
end if

end subroutine CHO_OPENVR
