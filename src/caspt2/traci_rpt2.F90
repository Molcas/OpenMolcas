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

subroutine TRACI_RPT2(ISTART,NDIM,XMAT,STSYM,NCI,CI)

use sguga, only: CIS, EXS, SGS
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISTART, NDIM, STSym, NCI
real(kind=wp), intent(inout) :: XMAT(NDIM,NDIM), CI(NCI)
integer(kind=iwp) :: I, IORB, J, JORB, LI, LJ, M
real(kind=wp) :: Fact, SCL, XJM
real(kind=wp), allocatable :: SGM(:), TVEC(:), XSAV(:,:)

if (NDIM <= 0) return

call mma_allocate(XSAV,NDIM,NDIM,Label='XSAV')
XSAV(:,:) = XMAT(:,:)
call mma_allocate(TVEC,NDIM,LABEL='TVEC')
call mma_allocate(SGM,NCI,LABEL='SGM')
SGM(:) = Zero

do J=1,NDIM
  FACT = One/XMAT(J,J)
  TVEC(:) = -FACT*XMAT(:,J)
  TVEC(J) = FACT
  XMAT(:,J) = Zero
  XMAT(J,J) = One
  ! Array T now contains a factor of XMAT of the form
  ! (e(1),..e(k-1),T,..,e(n)), where e(i) is the standard
  ! unit column vector, with elements Kronecker(l,i).
  ! Apply its inverse to XMAT.
  do M=J+1,NDIM
    XJM = XMAT(J,M)
    XMAT(:,M) = XMAT(:,M)+TVEC(:)*XJM
    XMAT(J,M) = TVEC(J)*XJM
  end do
  ! Transform CI array:
  ! CI:=( 1 + Sum(U(I)E(IJ)) + (1/2)Sum(U(I)U(M)E(IJ,MJ)) ) CI,
  ! where U(I) = T(I)-Kronecker(I,J).
  JORB = ISTART-1+J
  LJ = SGS%LEVEL(JORB)
  SGM(1:NCI) = (OneHalf-Half*TVEC(J))*CI(1:NCI)
  do I=1,NDIM
    IORB = ISTART-1+I
    LI = SGS%LEVEL(IORB)
    SCL = Half*TVEC(I)
    if (I == J) SCL = SCL-Half
    if (ABS(SCL)<1.0E-12_wp) cycle
    call SG_Epq_Psi(SGS,CIS,EXS,LI,LJ,SCL,STSYM,CI,SGM)
  end do
  do I=1,NDIM
    IORB = ISTART-1+I
    LI = SGS%LEVEL(IORB)
    SCL = TVEC(I)
    if (I == J) SCL = SCL-One
    if (ABS(SCL)<1.0E-12_wp) cycle
    call SG_Epq_Psi(SGS,CIS,EXS,LI,LJ,SCL,STSYM,SGM,CI)
  end do

end do

call mma_deallocate(SGM)
call mma_deallocate(TVEC)
XMAT(:,:) = XSAV(:,:)
call mma_deallocate(XSAV)

end subroutine TRACI_RPT2
