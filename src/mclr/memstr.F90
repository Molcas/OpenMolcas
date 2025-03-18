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

use Str_Info, only: STR, ITYP_DUMMY, NSTTYP, IUNIQMP, INDMAP, INUMAP, ISTAC, MNRS1, MNRS3, MXRS1, MXRS3, NELEC, NOCTYP, NSTFTP, &
                    NSTTYP, IUNIQTP
use stdalloc, only: mma_allocate
use MCLR_Data, only: NACOB, NOBPT, NORB1, NORB2, NORB3
use input_mclr, only: nIrrep

implicit none
! Local variables
integer LENGTH, ITYP, NSTRIN, LSTRIN, IITYP, IIIITEST, IMNEW, JJTYP, IANEQ, ICREQ, NCASTR_MCLR, NUMST3

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
allocate(Str(1:NSTTYP+1))

do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) then
    NSTRIN = NUMST3(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP))
    LSTRIN = NSTRIN*NELEC(ITYP)
    ! Offsets for occupation of strings and reordering array
    call mma_allocate(Str(ITYP)%OCSTR_Hidden,LSTRIN,Label='OCSTR')
    Str(ITYP)%OCSTR => Str(ITYP)%OCSTR_Hidden
    call mma_allocate(Str(ITYP)%STREO_Hidden,NSTRIN,Label='STREO')
    Str(ITYP)%STREO => Str(ITYP)%STREO_Hidden

    ! Symmetry and class of each string
    call mma_allocate(Str(ITYP)%STSM_Hidden,NSTRIN,Label='STSM')
    Str(ITYP)%STSM => Str(ITYP)%STSM_Hidden
    call mma_allocate(Str(ITYP)%STCL_Hidden,NSTRIN,Label='STCL')
    Str(ITYP)%STCL => Str(ITYP)%STCL_Hidden
  else
    IITYP = -IUNIQTP(ITYP)
    Str(ITYP)%OCSTR => Str(IITYP)%OCSTR_Hidden
    Str(ITYP)%STREO => Str(IITYP)%STREO_Hidden
    Str(ITYP)%STSM => Str(IITYP)%STSM_Hidden
    Str(ITYP)%STCL => Str(IITYP)%STCL_Hidden
  end if
end do

! Number of strings per symmetry and occupation
do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) then
    call mma_allocate(Str(ITYP)%NSTSO_Hidden,NOCTYP(ITYP)*nIrrep,Label='NSTSO')
    Str(ITYP)%NSTSO => Str(ITYP)%NSTSO_Hidden
    ! Offset of strings per symmetry and occupation
    call mma_allocate(Str(ITYP)%ISTSO_Hidden,NOCTYP(ITYP)*nIrrep,Label='ISTSO')
    Str(ITYP)%ISTSO => Str(ITYP)%ISTSO_Hidden
    ! Number of electrons in RAS1 and RAS3 per sub type, is sub-type active
    call mma_allocate(Str(ITYP)%EL1_Hidden,NOCTYP(ITYP),Label='EL1')
    Str(ITYP)%EL1 => Str(ITYP)%EL1_Hidden
    call mma_allocate(Str(ITYP)%EL3_Hidden,NOCTYP(ITYP),Label='EL3')
    Str(ITYP)%EL3 => Str(ITYP)%EL3_Hidden
    call mma_allocate(Str(ITYP)%ACTP_Hidden,NOCTYP(ITYP),Label='ACTP')
    Str(ITYP)%ACTP => Str(ITYP)%ACTP_Hidden
    !MS: New array introduced according to Jeppes new strinfo representation
    call mma_allocate(Str(ITYP)%EL123_Hidden,3*NOCTYP(ITYP),Label='EL123')
    Str(ITYP)%EL123 => Str(ITYP)%EL123_Hidden
    ! Lexical addressing of arrays: NB! Not allocated here in Jeppes new version!
    call mma_allocate(Str(ITYP)%Z_Hidden,NACOB*NELEC(ITYP),Label='Z')
    Str(ITYP)%Z => Str(ITYP)%Z_Hidden
  else
    ! redirect
    IITYP = -IUNIQTP(ITYP)
    Str(ITYP)%NSTSO => Str(IITYP)%NSTSO_Hidden
    Str(ITYP)%ISTSO => Str(IITYP)%ISTSO_Hidden
    Str(ITYP)%EL1 => Str(IITYP)%EL1_Hidden
    Str(ITYP)%EL3 => Str(IITYP)%EL3_Hidden
    Str(ITYP)%ACTP => Str(IITYP)%ACTP_Hidden
    Str(ITYP)%EL123 => Str(IITYP)%EL123_Hidden
    Str(ITYP)%Z => Str(IITYP)%Z_Hidden
  end if
end do

!MS: Introduced according to Jeppes new concept.
!MS: NB! Str(ITYP)%EL123 added to IEL13 parameter list!
!MS: Be aware that IEL13 is also called in STRINF
do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) &
    call IEL13(MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),NELEC(ITYP),NOCTYP(ITYP),Str(ITYP)%EL1,Str(ITYP)%EL3, &
               Str(ITYP)%EL123,Str(ITYP)%ACTP)
end do

! Mappings between different string types
do ITYP=1,NSTTYP
  !write(6,*) nelec(ityp),nstrin
  NSTRIN = NSTFTP(ITYP)

  if ((ISTAC(ITYP,2) /= 0) .and. (ISTAC(ITYP,1) /= 0)) then
    ! creation on string allowed, use full orbital notation
    LENGTH = NACOB*NSTRIN
    call mma_allocate(Str(ITYP)%STSTMI,1,Label='STSTMI')
    call mma_allocate(Str(ITYP)%STSTMN,1,Label='STSTMN')
  else if ((ISTAC(ITYP,1) /= 0) .and. (ISTAC(ITYP,2) == 0)) then

    ! only annihilation allowed, use compact scheme
    LENGTH = NELEC(ITYP)*NSTRIN
    call mma_allocate(Str(ITYP)%STSTMI,1,Label='STSTMI')
    call mma_allocate(Str(ITYP)%STSTMN,1,Label='STSTMN')
    !MS: New else block
  else if ((ISTAC(ITYP,1) == 0) .and. (ISTAC(ITYP,2) /= 0)) then
    ! Only creation allowed, use compact scheme with offsets

    call ICopy(NOCTYP(ITYP),[0],0,Str(ITYP)%NSTSO,1)
    call NUMST4_MCLR(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP),Str(ITYP)%NSTSO)
    LENGTH = NCASTR_MCLR(2,Str(ITYP)%NSTSO,NOCTYP(ITYP),ITYP,NOBPT,3,Str(ITYP)%EL123)
    ! Explicit offsets and lengths
    call mma_allocate(Str(ITYP)%STSTMI,NSTRIN,Label='STSTMI')
    call mma_allocate(Str(ITYP)%STSTMN,NSTRIN,Label='STSTMN')
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
          goto 1211
        end if

      end if
    end do
    ! Normal exit from DO loop only if no identical map was found
    IMNEW = 1
    IUNIQMP(ITYP) = ITYP
1211 continue
  end if

  if (IMNEW == 1) then
    call mma_allocate(Str(ITYP)%STSTM_Hidden,LENGTH,2,Label='STSTM')
    Str(ITYP)%STSTM => Str(ITYP)%STSTM_Hidden
  else
    IITYP = -IUNIQMP(ITYP)
    Str(ITYP)%STSTM => Str(IITYP)%STSTM_Hidden
  end if

end do

! Symmetry of conjugated orbitals and orbital excitations
!    COBSM,NIFSJ,IFSJ,IFSJO
!call mma_allocate(COBSM,NACOB,Label='COBSM')
!call mma_allocate(NIFSJ,NACOB*nIrrep,Label='NIFSJ')
!call mma_allocate(IFSJ,NACOB**2,Label='IFSJ')
!call mma_allocate(IFSJO,NACOB*nIrrep,Label='IFSJO')
! Symmetry of excitation connecting  strings of given symmetry
!call mma_allocate(STSTX,nIrrep*nIrrep,Label='STSTX')

! Up and down mappings of strings containing the same number of electrons

do ITYP=1,NSTTYP
  if (INUMAP(ITYP) /= 0) call mma_allocate(Str(ITYP)%NUMAP,NSTFTP(ITYP),Label='NUMAP')
  if (INDMAP(ITYP) /= 0) call mma_allocate(Str(ITYP)%NDMAP,NSTFTP(ITYP),Label='NDMAP')
end do

! Some dummy allocations

ITYP = ITYP_Dummy
call mma_allocate(Str(ITYP)%NSTSO_Hidden,1,Label='NSTSO')
Str(ITYP)%NSTSO => Str(ITYP)%NSTSO_Hidden
call mma_allocate(Str(ITYP)%EL1_Hidden,1,Label='EL1')
Str(ITYP)%EL1 => Str(ITYP)%EL1_Hidden
call mma_allocate(Str(ITYP)%EL3_Hidden,1,Label='EL3')
Str(ITYP)%EL3 => Str(ITYP)%EL3_Hidden

! Last word of string information
return

end subroutine MEMSTR
