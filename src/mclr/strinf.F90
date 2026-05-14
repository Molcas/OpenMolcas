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

!#define _DEBUGPRINT_
subroutine STRINF()
! Strings for internal space.
! Information is stored in
! Largest allowed length is MSTINF
!
! Input
! /LUCINP/,/ORBINP/,/CSM/
! Output
! /STRINP/,/STINF/,/STRBAS/ and string information in STIN

use Str_Info, only: ISTAC, IUNIQMP, IUNIQTP, MNRS1, MNRS3, MXRS1, MXRS3, NELEC, NOCTYP, NSTFTP, NSTTYP, STR
use MCLR_Data, only: NACOB, NORB1, NORB2, NORB3
use input_mclr, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxINT
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IMAX, ISGSTI(1), ISGSTO(1), ITYP, JTYP, LROW
integer(kind=iwp), allocatable :: KFREEL(:)

! 2 : Number of classes per string type and mappings between string types (/STINF/)

call ZSTINF_MCLR()

! 3 : Static memory for string information

call MEMSTR()

! 4 :Reverse lexical adresing schemes for each type of string

! First free address
call mma_MaxINT(imax)
call mma_allocate(KFREEL,imax,Label='KFREEL')
do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) &
    call WEIGHT_mclr(Str(ITYP)%Z,NELEC(ITYP),NORB1,NORB2,NORB3,MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),KFREEL)
end do

! 5 : Number of electrons in RAS1 and RAS3 per string sub type

do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) &
    call IEL13(MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),NELEC(ITYP),NOCTYP(ITYP),Str(ITYP)%EL1,Str(ITYP)%EL3,Str(ITYP)%EL123)
end do

! 6 : Number of strings per type and symmetry for a given string type

do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) then
    call NSTRSO_MCLR(NELEC(ITYP),NORB1,NORB2,NORB3,MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),KFREEL,Str(ITYP)%NSTSO, &
                     NOCTYP(ITYP),nIrrep,ITYP)
    ! Corresponding offset array
    call ZBASE(Str(ITYP)%NSTSO,Str(ITYP)%ISTSO,nIrrep*NOCTYP(ITYP))
    ! Symmetry and class index for each string
    call ZSMCL(nIrrep,NOCTYP(ITYP),Str(ITYP)%NSTSO,Str(ITYP)%STSM,Str(ITYP)%STCL)
  end if
end do

! 7 Construct strings, ordered according to symmetry and class

do ITYP=1,NSTTYP
  if (IUNIQTP(ITYP) == ITYP) &
    call GENSTR_MCLR(NELEC(ITYP),MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),Str(ITYP)%ISTSO,NOCTYP(ITYP),nIrrep,Str(ITYP)%Z, &
                     KFREEL,Str(ITYP)%STREO,Str(ITYP)%OCSTR,KFREEL(1+NOCTYP(ITYP)*nIrrep),ITYP)
end do

! 8 Internal annihilation arrays between types of strings

do ITYP=1,NSTTYP
  if (IUNIQMP(ITYP) == ITYP) then
    if (ISTAC(ITYP,1) /= 0) then
      JTYP = ISTAC(ITYP,1)
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Annihilator arrays between types ',ITYP,JTYP
      write(u6,*) ' ==========================================='
#     endif
      if (ISTAC(ITYP,2) == 0) then
        LROW = NELEC(ITYP)
      else
        LROW = NACOB
      end if
      Str(ITYP)%STSTM(1:NSTFTP(ITYP)*LROW,1) = 0
      Str(ITYP)%STSTM(1:NSTFTP(ITYP)*LROW,2) = 0
      call ANNSTR(Str(ITYP)%OCSTR,NSTFTP(ITYP),NSTFTP(JTYP),NELEC(ITYP),NACOB,Str(JTYP)%Z,Str(JTYP)%STREO,LROW,0,ISGSTI,ISGSTO, &
                  Str(ITYP)%STSTM(:,1),Str(ITYP)%STSTM(:,2),JTYP)
    end if
  end if
end do

! 6 : Creation arrays

do ITYP=1,NSTTYP
  if (IUNIQMP(ITYP) == ITYP) then
    if (ISTAC(ITYP,2) /= 0) then
      ! Type of creation map
      if (ISTAC(ITYP,1) == 0) then
        ! Only creation map, compact scheme with offsets
        LROW = -1
      else if (ISTAC(ITYP,1) /= 0) then
        ! Both annihilation and creation, use full form
        LROW = NACOB
      end if
      JTYP = ISTAC(ITYP,2)
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Creator  arrays between types ',ITYP,JTYP
      write(u6,*) ' ==========================================='
#     endif
      call CRESTR(Str(ITYP)%OCSTR,NSTFTP(ITYP),NSTFTP(JTYP),NELEC(ITYP),NACOB,Str(JTYP)%Z,Str(JTYP)%STREO,0,ISGSTI,ISGSTO, &
                  Str(ITYP)%STSTM(:,1),Str(ITYP)%STSTM(:,2),Str(ITYP)%STSTMN,Str(ITYP)%STSTMI,LROW,JTYP)
    end if
  end if
end do
call mma_deallocate(KFREEL)

return

end subroutine STRINF
