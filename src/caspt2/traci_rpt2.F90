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

use sguga, only: sg_epq_psi
use sguga_states, only: CIS, EXS, SGS
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, OneHalf
use Definitions, only: wp, iwp

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, MyRank, nProcs
#endif

implicit none
integer(kind=iwp), intent(in) :: ISTART, NDIM, STSym, NCI
real(kind=wp), intent(inout) :: XMAT(NDIM,NDIM), CI(NCI)
integer(kind=iwp) :: I, IORB, J, JORB, LI, LJ, M
real(kind=wp) :: Fact, SCL, XJM
real(kind=wp), allocatable :: SGM(:), TVEC(:), XSAV(:,:)
integer(kind=iwp), parameter :: istate = 1
real(kind=wp), parameter :: THRSCL = 1.0e-12_wp

#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: ITASK
real(kind=wp), allocatable :: CIACC(:)
#endif

if (NDIM <= 0) return

call mma_allocate(XSAV,NDIM,NDIM,Label='XSAV')
XSAV(:,:) = XMAT(:,:)
call mma_allocate(TVEC,NDIM,LABEL='TVEC')
call mma_allocate(SGM,NCI,LABEL='SGM')
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) call mma_allocate(CIACC,NCI,LABEL='CIACC')
#endif

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
  LJ = SGS(istate)%LEVEL(JORB)

# ifdef _MOLCAS_MPP_
  ITASK = 0
  SGM(1:NCI) = Zero
# else
  SGM(1:NCI) = (OneHalf-Half*TVEC(J))*CI(1:NCI)
# endif
  do I=1,NDIM
    IORB = ISTART-1+I
    LI = SGS(istate)%LEVEL(IORB)
    SCL = Half*TVEC(I)
    if (I == J) SCL = SCL-Half
    if (abs(SCL) < THRSCL) cycle
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      ITASK = ITASK+1
      if (mod(ITASK-1,nProcs) /= MyRank) cycle
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),LI,LJ,SCL,STSYM,CI,SGM)
    else
#   endif
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),LI,LJ,SCL,STSYM,CI,SGM)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end do
# ifdef _MOLCAS_MPP_
  call GADGOP(SGM,NCI,'+')
  SGM(1:NCI) = SGM(1:NCI)+(OneHalf-Half*TVEC(J))*CI(1:NCI)
# endif

  !--- Second half-transformation: SGM -> CI --------------------------
# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    CIACC(1:NCI) = Zero
    ITASK = 0
  end if
# endif
  do I=1,NDIM
    IORB = ISTART-1+I
    LI = SGS(istate)%LEVEL(IORB)
    SCL = TVEC(I)
    if (I == J) SCL = SCL-One
    if (abs(SCL) < THRSCL) cycle
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      ITASK = ITASK+1
      if (mod(ITASK-1,nProcs) /= MyRank) cycle
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),LI,LJ,SCL,STSYM,SGM,CIACC)
    else
#   endif
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),LI,LJ,SCL,STSYM,SGM,CI)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end do
# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    call GADGOP(CIACC,NCI,'+')
    CI(1:NCI) = CI(1:NCI)+CIACC(1:NCI)
  end if
# endif

end do

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) call mma_deallocate(CIACC)
#endif
call mma_deallocate(SGM)
call mma_deallocate(TVEC)
XMAT(:,:) = XSAV(:,:)
call mma_deallocate(XSAV)

end subroutine TRACI_RPT2
