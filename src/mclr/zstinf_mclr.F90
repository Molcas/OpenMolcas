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
subroutine ZSTINF_MCLR()
! Set up common block /STINF/ from information in /STINP/
!
!=========
! Input
!=========
! Information in /STINP/ and /ORBINP/
!
!====================
! Output (in /STINF/)
!====================
! ISTAC (MXPSTT,2) : string type obtained by creating (ISTAC(ITYP,2))
!                    or annihilating (ISTAC(ITYP,1)) an electron
!                    from a string of type  ITYP. A zero indicates
!                    that this mapping is not included
!                    Only strings having the same ISTTP index are
!                    mapped
! NOCTYP(ITYP) : Number of occupation classes for given type
!
!
! NSTFTP(ITYP) : Number of strings of this type
!
!   / \           Zero order space                         !
!    !            Double excitations from reference space  !  Down
! Up !            single excitation from reference space   !
!    !            reference space                         \ /

use Str_Info, only: ISTAC, MNRS1, MNRS3, MXRS1, MXRS3, NELEC, NOCTYP, NSTFTP, NSTTYP
use MCLR_Data, only: NORB1, NORB2, NORB3
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: ITYP, NUMST3

! *****************************************************************
! Mappings between strings with the same type ISTTP index, +/- 1 el
! *****************************************************************
ISTAC(:,:) = 0
do ITYP=1,NSTTYP-1
  if (NELEC(ITYP+1) == NELEC(ITYP)-1) then
    ISTAC(ITYP,1) = ITYP+1
    ISTAC(ITYP+1,2) = ITYP
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Type - type mapping array ISTAC'
write(u6,*) ' ==============================='
call IWRTMA(ISTAC,NSTTYP,2,size(ISTAC,1),2)
#endif
! *************************************************
! Number of occupation classes and strings per type
! *************************************************
NOCTYP(1:NSTTYP) = (MXRS1(1:NSTTYP)-MNRS1(1:NSTTYP)+1)*(MXRS3(1:NSTTYP)-MNRS3(1:NSTTYP)+1)
#ifdef _DEBUGPRINT_
write(u6,*) ' Number of occupation classes per type'
write(u6,*) ' ====================================='
call IWRTMA(NOCTYP,1,NSTTYP,1,NSTTYP)
#endif

do ITYP=1,NSTTYP
  NSTFTP(ITYP) = NUMST3(NELEC(ITYP),NORB1,MNRS1(ITYP),MXRS1(ITYP),NORB2,NORB3,MNRS3(ITYP),MXRS3(ITYP))
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' Number of strings per  type'
write(u6,*) ' ==========================='
call IWRTMA(NSTFTP,1,NSTTYP,1,NSTTYP)
#endif
! ****************************************************************
! Mappings between strings containing the same number of electrons
! ****************************************************************
!INUMAP(:) = 0
!INDMAP(:) = 0
! Mapping to and from zero order space
! Note: some lines are commented out here since IARTP and IBRTP
!       have never been defined. (R. Lindh 2006)
!INUMAP(IARTP(3,5)) = IAZTP
!if (IARTP(3,4) /= 0) INUMAP(IARTP(3,4)) = IAZTP+1
!if (IARTP(3,3) /= 0) INUMAP(IARTP(3,3)) = IAZTP+2

!INUMAP(IBRTP(3,5)) = IBZTP
!if (IBRTP(3,4) /= 0) INUMAP(IBRTP(3,4)) = IBZTP+1
!if (IBRTP(3,3) /= 0) INUMAP(IBRTP(3,3)) = IBZTP+2

!NAEL = NELEC(IAZTP)
!INDMAP(IAZTP) = IARTP(3,5)
!if (NAEL >= 1) INDMAP(IAZTP+1) = IARTP(3,4)
!if (NAEL >= 2) INDMAP(IAZTP+2) = IARTP(3,3)

!NBEL = NELEC(IBZTP)
!INDMAP(IBZTP) = IBRTP(3,5)
!if (NBEL >= 1) INDMAP(IBZTP+1) = IBRTP(3,4)
!if (NBEL >= 2) INDMAP(IBZTP+2) = IBRTP(3,3)

! Number of electrons compared to reference
!do IDEL=-4,2
!  do IEX=1,2
!    ! Up mappings
!    if (IARTP(IEX,IDEL+5) /= 0) INUMAP(IARTP(IEX,IDEL+5)) = IARTP(IEX+1,IDEL+5)
!    if (IBRTP(IEX,IDEL+5) /= 0) INUMAP(IBRTP(IEX,IDEL+5)) = IBRTP(IEX+1,IDEL+5)
!    ! Down mappings
!    if (IARTP(IEX+1,IDEL+5) /= 0) INDMAP(IARTP(IEX+1,IDEL+5)) = IARTP(IEX,IDEL+5)
!    if (IBRTP(IEX+1,IDEL+5) /= 0) INDMAP(IBRTP(IEX+1,IDEL+5)) = IBRTP(IEX,IDEL+5)
!  end do
!end do

!#ifdef _DEBUGPRINT_
!write(u6,*) ' Up mappings of string types'
!call IWRTMA(INUMAP,1,NSTTYP,1,NSTTYP)
!write(u6,*) ' Down mappings of string types'
!call IWRTMA(INDMAP,1,NSTTYP,1,NSTTYP)
!#endif

return

end subroutine ZSTINF_MCLR
