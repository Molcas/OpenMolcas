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

subroutine SSOTRA(SGS,CIS,EXS,ISYM,LSM,NA,NO,TRA,NCO,CI,TMP)

use gugx, only: CIStruct, EXStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Half
use Definitions, only: wp, iwp

implicit none
type(SGSTruct) :: SGS
type(CISTruct) :: CIS
type(EXSTruct) :: ExS
integer(kind=iwp) :: ISYM, LSM, NA, NO, NCO
real(kind=wp) :: TRA(NO,NO), CI(NCO), TMP(NCO)
integer(kind=iwp) :: IK, IKLEV, IL, IP, IPLEV, NI
real(kind=wp) :: CKK, CPK, X
integer(kind=iwp), allocatable :: ILEV(:)

! ILEV(IORB)=GUGA LEVEL CORRESPONDING TO A SPECIFIC ACTIVE ORBITAL
! OF SYMMETRY ISYM.
call mma_allocate(ILEV,NA,Label='ILEV')
NI = NO-NA
IL = 0
do IP=1,NA
  do
    IL = IL+1
    if (SGS%ISM(IL) == ISYM) exit
  end do
  ILEV(IP) = IL
end do
!TEST write(u6,*)' Check prints in SSOTRA.'
!TEST write(u6,*)' ISYM:',ISYM
do IK=1,NA
  IKLEV = ILEV(IK)
  TMP(:) = Zero
  do IP=1,NA
    IPLEV = ILEV(IP)
    CPK = TRA(NI+IP,NI+IK)
    if (IP == IK) CPK = CPK-One
    X = Half*CPK
    !TEST write(u6,*)' IP,IK,X:',IP,IK,X
    if (abs(X) < 1.0e-14_wp) cycle
    call SIGMA1(SGS,CIS,EXS,IPLEV,IKLEV,X,LSM,CI,TMP)
  end do
  CKK = TRA(NI+IK,NI+IK)
  X = Three-CKK
  CI(:) = CI(:)+X*TMP(:)
  do IP=1,NA
    IPLEV = ILEV(IP)
    CPK = TRA(NI+IP,NI+IK)
    if (IP == IK) CPK = CPK-One
    if (abs(CPK) < 1.0e-14_wp) cycle
    call SIGMA1(SGS,CIS,EXS,IPLEV,IKLEV,CPK,LSM,TMP,CI)
  end do
end do
call mma_deallocate(ILEV)

end subroutine SSOTRA
