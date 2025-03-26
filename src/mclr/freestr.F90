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
! Copyright (C) 1990,1994,1995, Jeppe Olsen                            *
!***********************************************************************

subroutine FREESTR()
! Free pointers for saving information about strings and
! their mappings
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
! Modified for deallocation, Sept. 25, 2005.

use Str_Info, only: NSTTYP, STR, ITYP_DUMMY, IUNIQMP, INDMAP, INUMAP, ISTAC, IUNIQTP
use stdalloc, only: mma_deallocate

implicit none
integer ITYP, IITYP, IIIITEST, IMNEW, JJTYP, IANEQ, ICREQ

! Start of string information
! =====================================================================
!
! 1 : String information
!
! =====================================================================

do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) then
    ! Offsets for occupation of strings and reordering array
    call mma_deallocate(Str(ITYP)%OCSTR_Hidden)
    call mma_deallocate(Str(ITYP)%STREO_Hidden)
    ! Symmetry and class of each string
    call mma_deallocate(Str(ITYP)%STSM_Hidden)
    call mma_deallocate(Str(ITYP)%STCL_Hidden)
  end if
  nullify(Str(ITYP)%OCSTR,Str(ITYP)%STREO,Str(ITYP)%STSM,Str(ITYP)%STCL)
end do

! Number of strings per symmetry and occupation
do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) then
    call mma_deallocate(Str(ITYP)%NSTSO_Hidden)
    ! Offset of strings per symmetry and occupation
    call mma_deallocate(Str(ITYP)%ISTSO_Hidden)
    ! Number of electrons in RAS1 and RAS3 per sub type, is sub-type active
    call mma_deallocate(Str(ITYP)%EL1_Hidden)
    call mma_deallocate(Str(ITYP)%EL3_Hidden)
    call mma_deallocate(Str(ITYP)%ACTP_Hidden)
    !MS: New array introduced according to Jeppes new strinfo representation
    call mma_deallocate(Str(ITYP)%EL123_Hidden)
    ! Lexical addressing of arrays: NB! Not allocated here in Jeppes new version!
    call mma_deallocate(Str(ITYP)%Z_Hidden)
  else
    ! redirect
    IITYP = -IUNIQTP(ITYP)
  end if
  nullify(Str(ITYP)%NSTSO,Str(ITYP)%ISTSO,Str(ITYP)%EL1,Str(ITYP)%EL3,Str(ITYP)%ACTP,Str(ITYP)%EL123,Str(ITYP)%Z)
end do

! Mappings between different string types
do ITYP=1,NSTTYP
  if ((ISTAC(ITYP,2) /= 0) .and. (ISTAC(ITYP,1) /= 0)) then
    ! creation on string allowed, use full orbital notation
    call mma_deallocate(Str(ITYP)%STSTMI)
    call mma_deallocate(Str(ITYP)%STSTMN)
  else if ((ISTAC(ITYP,1) /= 0) .and. (ISTAC(ITYP,2) == 0)) then

    ! only annihilation allowed, use compact scheme
    call mma_deallocate(Str(ITYP)%STSTMI)
    call mma_deallocate(Str(ITYP)%STSTMN)
  else if ((ISTAC(ITYP,1) == 0) .and. (ISTAC(ITYP,2) /= 0)) then
    !MS: New else block
    ! Only creation allowed, use compact scheme with offsets

    ! Explicit offsets and lengths
    call mma_deallocate(Str(ITYP)%STSTMI)
    call mma_deallocate(Str(ITYP)%STSTMN)
  end if
  ! has this map been constructed before ?
  IIIITEST = 0
  if ((IUNIQTP(ITYP) == ITYP) .or. (IIIITEST == 1)) then
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
          IMNEW = 0
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
  if (IMNEW == 1) call mma_deallocate(Str(ITYP)%STSTM_Hidden)
  nullify(Str(ITYP)%STSTM)
end do
! Symmetry of conjugated orbitals and orbital excitations
!   COBSM,NIFSJ,IFSJ,IFSJO
!call mma_deallocate(COBSM)
!call mma_deallocate(NIFSJ)
!call mma_deallocate(IFSJ)
!call mma_deallocate(IFSJO)
! Symmetry of excitation connecting  strings of given symmetry
!call mma_deallocate(STSTX)

! Up and down mappings of strings containing the same number of electrons

do ITYP=1,NSTTYP
  if (INUMAP(ITYP) /= 0) call mma_deallocate(Str(ITYP)%NUMAP)
  if (INDMAP(ITYP) /= 0) call mma_deallocate(Str(ITYP)%NDMAP)
end do

! Some dummy dallocations

ITYP = ITYP_Dummy
call mma_deallocate(Str(ITYP)%NSTSO_Hidden)
call mma_deallocate(Str(ITYP)%EL1_Hidden)
call mma_deallocate(Str(ITYP)%EL3_Hidden)
nullify(Str(ITYP)%NSTSO,Str(ITYP)%EL1,Str(ITYP)%EL3)

end subroutine FREESTR
