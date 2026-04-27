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
! Copyright (C) 1989,1998, Per Ake Malmqvist                           *
!***********************************************************************
!  PROGRAM RASSI        PER-AAKE MALMQVIST
!  SUBROUTINE CITRA     IBM-3090 RELEASE 89 01 31
!  USE THE COEFFICIENTS FOR A SEQUENCE OF SINGLE-ORBITAL TRANSFOR-
!  MATION, TRA, TO TRANSFORM THE CI EXPANSION COEFFICIENTS
!  IN-PLACE TO A NEW NON-ON ORBITAL BASIS.
!  NEW VERSION 981122, using user define types SGS,CIS,XS.
!***********************************************************************
!  CITRA
!
!> @brief
!>   Recompute a CI coefficient array to use another orbital basis
!> @author P. &Aring;. Malmqvist
!>
!> @details
!> For a given linear transformation of the orbitals, and a
!> CI array where the CSF basis is built from the old
!> orbitals, compute the CI array using the new orbitals
!> instead. The orbitals are transformed sequentially, and
!> for each active orbital, a call to ::SSOTRA performs the
!> single-orbital transformation.
!>
!> @param[in]     WFTP Wave function Type Name
!> @param[in]     SGS Split Graph Structure user defined type
!> @param[in]     CIS CI Structure user define type
!> @param[in,out] EXS  Excitation operator Structure user defined type
!> @param[in]     LSM  Wave function Symmetry Label
!> @param[in]     TRA  Transformation Matrix
!> @param[in]     NCO  Number of Configuration Functions
!> @param[in,out] CI   CI Array
!***********************************************************************

!ifdef _DEBUGPRINT_
subroutine CITRA(WFTP,SGS,CIS,EXS,LSM,TRA,NCO,CI)

use gugx, only: SGStruct, CIStruct, EXStruct
use Symmetry_Info, only: nIrrep
use rassi_data, only: NTRA, NOSH, NISH, NASH
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
character(len=8), intent(in) :: WFTP
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: LSM, NCO
real(kind=wp), intent(in) :: TRA(NTRA)
real(kind=wp), intent(inout) :: CI(NCO)
integer(kind=iwp) :: I, II, ISTA, ISYM, NA, NI, NO
real(kind=wp) :: CKK, FAC
real(kind=wp), allocatable :: TMP(:)
#ifdef _DEBUGPRINT_
real(kind=wp), external :: ddot_
#endif

#ifdef _DEBUGPRINT_
write(u6,*) ' Entering CITRA. norm=',ddot_(NCO,CI,1,CI,1)
#endif
!write(u6,*) ' Entering CITRA. TRA='
!write(u6,'(1x,5f16.8)') (TRA(I),I=1,NTRA)
!write(u6,*) ' Entering CITRA. CI='
!write(u6,'(1x,5f16.8)') (CI(I),I=1,NCO)
! TRA contains square matrices, one per symmetry
!  FIRST TRANSFORM THE INACTIVE ORBITALS:
FAC = One
ISTA = 1
do ISYM=1,nIrrep
  NO = NOSH(ISYM)
  do I=1,NISH(ISYM)
    II = ISTA+(NO+1)*(I-1)
    CKK = TRA(II)
    FAC = FAC*CKK
  end do
  ISTA = ISTA+NO**2
end do
!write(u6,*) 'FAC, FAC**2 ... ',FAC,FAC**2
FAC = FAC**2
CI(:) = FAC*CI(:)
!write(u6,*) ' CITRA. inactive done CI='
!write(u6,'(1x,5f16.8)') (CI(I),I=1,NCO)
! THEN THE ACTIVE ONES:
if (WFTP /= 'EMPTY') then
  if ((WFTP == 'HISPIN') .or. (WFTP == 'CLOSED')) then
    ! The HISPIN case may be buggy and is not presently used.
    ISTA = 1
    do ISYM=1,nIrrep
      NI = NISH(ISYM)
      NA = NASH(ISYM)
      NO = NOSH(ISYM)
      do I=NI+1,NI+NA
        II = ISTA+(NO+1)*(I-1)
        CKK = TRA(II)
        FAC = FAC*CKK
      end do
      ISTA = ISTA+NO**2
    end do
    if (WFTP == 'CLOSED') FAC = FAC**2
    CI(:) = FAC*CI(:)
  else
    ! The general case:
    call mma_allocate(TMP,NCO,Label='TMP')
    ISTA = 1
    do ISYM=1,nIrrep
      NA = NASH(ISYM)
      NO = NOSH(ISYM)
      if (NA /= 0) call SSOTRA(SGS,CIS,EXS,ISYM,LSM,NA,NO,TRA(ISTA),NCO,CI,TMP)
      ISTA = ISTA+NO**2
    end do
    call mma_deallocate(TMP)
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) ' DONE in  CITRA. norm=',ddot_(NCO,CI,1,CI,1)
# endif
  !write(u6,*) ' CITRA completely done. CI='
  !write(u6,'(1x,5f16.8)') (CI(I),I=1,NCO)

end if

end subroutine CITRA
