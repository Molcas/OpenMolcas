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

subroutine CHO_SETADDR(INFRED,INFVEC,MRED,MVEC,M2,MSYM)
!
! Purpose: set first disk addresses for reduced set info and
!          vectors.

use Cholesky, only: Cho_AdrVec, LuCho, nnBstR, nnBstRT, nnShl, nSym, NumCho, XnPass
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: MRED, MVEC, M2, MSYM
integer(kind=iwp), intent(inout) :: INFRED(MRED), INFVEC(MVEC,M2,MSYM)
integer(kind=iwp) :: IADR, IOPT, IPASS, IRED, ISYM, JPASS, LSA
real(kind=wp), allocatable :: KSA(:)
character(len=*), parameter :: SECNAM = 'CHO_SETADDR'

! Set addresses.
! --------------

if (XNPASS == 0) then
  INFRED(1) = 0
  INFVEC(1,3:4,1:NSYM) = 0
else if (XNPASS > 0) then
  IRED = 3
  IPASS = XNPASS
  call CHO_GETRED(IPASS,IRED,.false.)
  call CHO_SETREDIND(IRED)
  if (IPASS == 1) then
    INFRED(IPASS+1) = INFRED(IPASS)+NSYM*NNSHL+2*NNBSTRT(IRED)+NNSHL
  else
    INFRED(IPASS+1) = INFRED(IPASS)+NSYM*NNSHL+NNBSTRT(IRED)
  end if
  do ISYM=1,NSYM
    if (NUMCHO(ISYM) == 0) then
      INFVEC(1,3:4,ISYM) = 0
    else if (NUMCHO(ISYM) > 0) then
      if (CHO_ADRVEC == 1) then
        JPASS = INFVEC(NUMCHO(ISYM),2,ISYM)
        if (JPASS == IPASS) then
          INFVEC(NUMCHO(ISYM)+1,3,ISYM) = INFVEC(NUMCHO(ISYM),3,ISYM)+NNBSTR(ISYM,IRED)
          INFVEC(NUMCHO(ISYM)+1,4,ISYM) = INFVEC(NUMCHO(ISYM),4,ISYM)+NNBSTR(ISYM,IRED)
        else if ((JPASS <= XNPASS) .and. (JPASS > 0)) then
          IPASS = JPASS
          call CHO_GETRED(IPASS,IRED,.false.)
          call CHO_SETREDIND(IRED)
          INFVEC(NUMCHO(ISYM)+1,3,ISYM) = INFVEC(NUMCHO(ISYM),3,ISYM)+NNBSTR(ISYM,IRED)
          INFVEC(NUMCHO(ISYM)+1,4,ISYM) = INFVEC(NUMCHO(ISYM),4,ISYM)+NNBSTR(ISYM,IRED)
        else
          call CHO_QUIT('[1] JPASS error in '//SECNAM,104)
        end if
      else if (CHO_ADRVEC == 2) then
        JPASS = INFVEC(NUMCHO(ISYM),2,ISYM)
        if (JPASS == IPASS) then
          LSA = NNBSTR(ISYM,IRED)
          call mma_allocate(KSA,LSA,Label='KSA')
          IOPT = 2
          IADR = INFVEC(NUMCHO(ISYM),3,ISYM)
          call DDAFILE(LUCHO(ISYM),IOPT,KSA,LSA,IADR)
          INFVEC(NUMCHO(ISYM)+1,3,ISYM) = IADR
          INFVEC(NUMCHO(ISYM)+1,4,ISYM) = INFVEC(NUMCHO(ISYM),4,ISYM)+NNBSTR(ISYM,IRED)
          call mma_deallocate(KSA)
        else if ((JPASS <= XNPASS) .and. (JPASS > 0)) then
          IPASS = JPASS
          call CHO_GETRED(IPASS,IRED,.false.)
          call CHO_SETREDIND(IRED)
          LSA = NNBSTR(ISYM,IRED)
          call mma_allocate(KSA,LSA,Label='KSA')
          IOPT = 2
          IADR = INFVEC(NUMCHO(ISYM),3,ISYM)
          call DDAFILE(LUCHO(ISYM),IOPT,KSA,LSA,IADR)
          INFVEC(NUMCHO(ISYM)+1,3,ISYM) = IADR
          INFVEC(NUMCHO(ISYM)+1,4,ISYM) = INFVEC(NUMCHO(ISYM),4,ISYM)+NNBSTR(ISYM,IRED)
          call mma_deallocate(KSA)
        else
          call CHO_QUIT('[2] JPASS error in '//SECNAM,104)
        end if
      else
        call CHO_QUIT('CHO_ADRVEC error in '//SECNAM,102)
      end if
    else
      call CHO_QUIT('NUMCHO error in '//SECNAM,104)
    end if
  end do
else
  call CHO_QUIT('XNPASS error in '//SECNAM,104)
end if

end subroutine CHO_SETADDR
