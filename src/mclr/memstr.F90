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
! Copyright (C) 1990,1994, Jeppe Olsen                                 *
!***********************************************************************

subroutine MEMSTR()
! Construct pointers for saving information about strings and their mappings
!
!========
! Input :
!========
! Number and types of strings defined by /STRINP/
! Symmetry information stored in         /CSM/
! String information stored in           /STINF/
!=========
! Output
!=========
! Pointers stored in common block /STRBAS/
!
! Jeppe Olsen, Winter of 1990
!
! Last Revision, Dec 24 1990, Almaden
!
! Updated with iuniqtp, dec 11, 1994

use Str_Info, only: ISTAC, ITYP_DUMMY, IUNIQMP, IUNIQTP, MNRS1, MNRS3, MXRS1, MXRS3, NELEC, NOCTYP, NSTFTP, NSTTYP, NSTTYP, Str, &
                    Str_Hidden
use MCLR_Data, only: NACOB, NOBPT, NORB1, NORB2, NORB3
use input_mclr, only: nIrrep
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IANEQ, ICREQ, IITYP, IMNEW, ITYP, JJTYP, LENGTH, LSTRIN, NSTRIN, NUMST3
integer(kind=iwp), external :: NCASTR_MCLR

LENGTH = 100
! Start of string information

! =====================================================================
!
! 1 : String information
!
! =====================================================================

! Some calls are done with points which are out of bounds. To make
! this to be strictly secure we add a dummy layer and point to that
! in the case that we are out of bounds. The code seems, however,
! not to touch these arrays.

ITYP_Dummy = NSTTYP+1

do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) then
    NSTRIN = NUMST3(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP))
    LSTRIN = NSTRIN*NELEC(ITYP)
    ! Offsets for occupation of strings and reordering array
    call mma_allocate(Str_Hidden(ITYP)%OCSTR,LSTRIN,Label='OCSTR')
    call mma_allocate(Str_Hidden(ITYP)%STREO,NSTRIN,Label='STREO')

    ! Symmetry and class of each string
    call mma_allocate(Str_Hidden(ITYP)%STSM,NSTRIN,Label='STSM')
    call mma_allocate(Str_Hidden(ITYP)%STCL,NSTRIN,Label='STCL')

    ! Number of strings per symmetry and occupation
    call mma_allocate(Str_Hidden(ITYP)%NSTSO,NOCTYP(ITYP)*nIrrep,Label='NSTSO')
    ! Offset of strings per symmetry and occupation
    call mma_allocate(Str_Hidden(ITYP)%ISTSO,NOCTYP(ITYP)*nIrrep,Label='ISTSO')
    ! Number of electrons in RAS1 and RAS3 per sub type, is sub-type active
    call mma_allocate(Str_Hidden(ITYP)%EL1,NOCTYP(ITYP),Label='EL1')
    call mma_allocate(Str_Hidden(ITYP)%EL3,NOCTYP(ITYP),Label='EL3')
    !MS: New array introduced according to Jeppe's new strinfo representation
    call mma_allocate(Str_Hidden(ITYP)%EL123,3*NOCTYP(ITYP),Label='EL123')
    ! Lexical addressing of arrays: NB! Not allocated here in Jeppe's new version!
    call mma_allocate(Str_Hidden(ITYP)%Z,NACOB*NELEC(ITYP),Label='Z')

    !MS: Introduced according to Jeppe's new concept.
    !MS: NB! Str(ITYP)%EL123 added to IEL13 parameter list!
    !MS: Be aware that IEL13 is also called in STRINF
    call IEL13(MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),NELEC(ITYP),NOCTYP(ITYP),Str_Hidden(ITYP)%EL1,Str_Hidden(ITYP)%EL3, &
               Str_Hidden(ITYP)%EL123)

    IITYP = ITYP
  else
    IITYP = -IUNIQTP(ITYP)
  end if
  Str(ITYP)%OCSTR => Str_Hidden(IITYP)%OCSTR
  Str(ITYP)%STREO => Str_Hidden(IITYP)%STREO
  Str(ITYP)%STSM => Str_Hidden(IITYP)%STSM
  Str(ITYP)%STCL => Str_Hidden(IITYP)%STCL
  Str(ITYP)%NSTSO => Str_Hidden(IITYP)%NSTSO
  Str(ITYP)%ISTSO => Str_Hidden(IITYP)%ISTSO
  Str(ITYP)%EL1 => Str_Hidden(IITYP)%EL1
  Str(ITYP)%EL3 => Str_Hidden(IITYP)%EL3
  Str(ITYP)%EL123 => Str_Hidden(IITYP)%EL123
  Str(ITYP)%Z => Str_Hidden(IITYP)%Z
end do

! Mappings between different string types
do ITYP=1,NSTTYP
  !write(u6,*) nelec(ityp),nstrin

  if ((ISTAC(ITYP,2) /= 0) .and. (ISTAC(ITYP,1) /= 0)) then
    ! creation on string allowed, use full orbital notation
    LENGTH = NACOB*NSTFTP(ITYP)
    NSTRIN = 1
  else if ((ISTAC(ITYP,1) /= 0) .and. (ISTAC(ITYP,2) == 0)) then
    ! only annihilation allowed, use compact scheme
    LENGTH = NELEC(ITYP)*NSTFTP(ITYP)
    NSTRIN = 1
  else if ((ISTAC(ITYP,1) == 0) .and. (ISTAC(ITYP,2) /= 0)) then
    !MS: New else block
    ! Only creation allowed, use compact scheme with offsets
    Str(ITYP)%NSTSO(:) = 0
    call NUMST4_MCLR(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP),Str(ITYP)%NSTSO)
    LENGTH = NCASTR_MCLR(2,Str(ITYP)%NSTSO,NOCTYP(ITYP),ITYP,NOBPT,3,Str(ITYP)%EL123)
    NSTRIN = NSTFTP(ITYP)
  else
    NSTRIN = 0
  end if
  ! Explicit offsets and lengths
  call mma_allocate(Str_Hidden(ITYP)%STSTMI,NSTRIN,Label='STSTMI')
  call mma_allocate(Str_Hidden(ITYP)%STSTMN,NSTRIN,Label='STSTMN')

  ! has this map been constructed before ?
  IMNEW = 0
  if (IUNIQTP(ITYP) == ITYP) then
    IMNEW = 1
    IUNIQMP(ITYP) = ITYP
  else
    ! check type of previous map
    do JJTYP=1,ITYP-1
      IITYP = -IUNIQTP(ITYP)
      if ((abs(IUNIQTP(JJTYP)) == IITYP) .and. (IUNIQMP(JJTYP) == JJTYP)) then
        if (((ISTAC(ITYP,1) == 0) .and. (ISTAC(JJTYP,1) == 0)) .or. &
            ((ISTAC(ITYP,1) /= 0) .and. (ISTAC(JJTYP,1) /= 0) .and. &
             (abs(IUNIQTP(ISTAC(ITYP,1))) == abs(IUNIQTP(ISTAC(JJTYP,1)))))) then
          IANEQ = 1
        else
          IANEQ = 0
        end if
        if (((ISTAC(ITYP,2) == 0) .and. (ISTAC(JJTYP,2) == 0)) .or. &
            ((ISTAC(ITYP,2) /= 0) .and. (ISTAC(JJTYP,2) /= 0) .and. &
             (abs(IUNIQTP(ISTAC(ITYP,2))) == abs(IUNIQTP(ISTAC(JJTYP,2)))))) then
          ICREQ = 1
        else
          ICREQ = 0
        end if
        if ((IANEQ == 1) .and. (ICREQ == 1)) then
          IUNIQMP(ITYP) = -JJTYP
          exit
        end if

      end if
    end do
    ! Normal exit from DO loop only if no identical map was found
    if (JJTYP == ITYP) then
      IMNEW = 1
      IUNIQMP(ITYP) = ITYP
    end if
  end if

  ! note: not IITYP here
  Str(ITYP)%STSTMI => Str_Hidden(ITYP)%STSTMI
  Str(ITYP)%STSTMN => Str_Hidden(ITYP)%STSTMN

  if (IMNEW == 1) then
    call mma_allocate(Str_Hidden(ITYP)%STSTM,LENGTH,2,Label='STSTM')
    IITYP = ITYP
  else
    IITYP = -IUNIQMP(ITYP)
  end if

  Str(ITYP)%STSTM => Str_Hidden(IITYP)%STSTM

end do

! Symmetry of conjugated orbitals and orbital excitations
!    COBSM,NIFSJ,IFSJ,IFSJO
!call mma_allocate(COBSM,NACOB,Label='COBSM')
!call mma_allocate(NIFSJ,NACOB*nIrrep,Label='NIFSJ')
!call mma_allocate(IFSJ,NACOB**2,Label='IFSJ')
!call mma_allocate(IFSJO,NACOB*nIrrep,Label='IFSJO')
! Symmetry of excitation connecting  strings of given symmetry
!call mma_allocate(STSTX,nIrrep*nIrrep,Label='STSTX')

! Some dummy allocations

ITYP = ITYP_Dummy
call mma_allocate(Str_Hidden(ITYP)%NSTSO,1,Label='NSTSO')
call mma_allocate(Str_Hidden(ITYP)%EL1,1,Label='EL1')
call mma_allocate(Str_Hidden(ITYP)%EL3,1,Label='EL3')
Str(ITYP)%NSTSO => Str_Hidden(ITYP)%NSTSO
Str(ITYP)%EL1 => Str_Hidden(ITYP)%EL1
Str(ITYP)%EL3 => Str_Hidden(ITYP)%EL3

end subroutine MEMSTR
