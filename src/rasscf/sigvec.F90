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

subroutine SIGVEC(CIN,HC,HD,BM,SXN,G,H,DIA,F1,F2,X,C,NTRIAL)
! RASSCF Program version IBM-3090: SX section
!
! Purpose: Calculation of the SIGMA vector HC for a super-CI
! Hamiltonian.
! called from HMAT
!PAM01 Added miniscule constant times CIN to HC.
!
! ********** IBM-3090 Release 88 09 01 **********

use rasscf_global, only: ICICP, ITER, NDIMSX, NROOT, NSXS, SXSHFT, IROOT, ENER
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: NSYM, NASH, NISH, NSSH, SXDAMP
use Constants, only: Zero, One
use Definitions, only: wp, u6

implicit none
real*8 CIN(*), HC(*), BM(*), SXN(*), G(*), H(*), DIA(*), F1(*), F2(*), X(*), C(*)
integer NTRIAL
character(len=16), parameter :: ROUTINE = 'SIGVEC  '
real*8 HD(*)
integer :: I, iPrLev, ISTAE, ISTBM, ISTH, ISTIA, ISTZ, ISYM, ITRIAL, NAE, NAO, NEO, NIA, NIO, NNST, NST

! Local print level (if any)
IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering ',ROUTINE
NST = 0
do ITRIAL=1,NTRIAL
  NNST = NROOT+NST

  ! renormalize the C vector

  do I=1,NSXS
    C(I) = SXN(I)*CIN(I+NNST)
  end do

  ! Remove any unwanted rotations from C:
  do I=1,NSXS
    if (HD(I+NROOT) > 1.0e20_wp) C(I) = Zero
  end do

  ! Initialize sigma vector to zero.
  call DCOPY_(NROOT+NSXS,[Zero],0,HC(NST+1),1)

  ISTIA = 1
  ISTAE = 1
  ISTBM = 1
  ISTH = 1
  ISTZ = 0
  do ISYM=1,NSYM
    NIO = NISH(ISYM)
    NAO = NASH(ISYM)
    NEO = NSSH(ISYM)
    NIA = NIO+NAO
    NAE = NAO+NEO
    if ((NIA == 0) .or. (NAE == 0)) GO TO 98

    ! G-matrix contribution to HC (G*C)

    call DGEMM_('N','N',NIA,NAE,NIA,One,G(ISTIA),NIA,C(ISTBM+NST),NIA,One,HC(ISTBM+NNST),NIA)

    ! H-matrix contribution to HC (-C*H)

    if (NAO /= 0) call DGEMM_('N','N',NIA,NAO,NAE,-One,C(ISTBM+NST),NIA,H(ISTH),NAE,One,HC(ISTBM+NNST),NIA)
    if (NAO*NEO /= 0) call DGEMM_('N','T',NIA,NEO,NAO,-One,C(ISTBM+NST),NIA,H(ISTH+NAO),NAE,One,HC(ISTBM+NAO*NIA+NNST),NIA)

    ! First Fock matrix contribution D*C*FP

    call DGEMM_('N','N',NIA,NAE,NAE,One,C(ISTBM+NST),NIA,F2(ISTAE),NAE,Zero,X,NIA)
    call DGEMM_('N','N',NIA,NAE,NIA,One,DIA(ISTIA),NIA,X,NIA,One,HC(ISTBM+NNST),NIA)

    ! Second Fock matrix contribution FP*C*D

    if (NAO /= 0) then
      call DGEMM_('N','N',NIA,NAO,NIA,One,F1(ISTIA),NIA,C(ISTBM+NST),NIA,Zero,X,NIA)
      call DGEMM_('N','N',NIA,NAO,NAO,One,X,NIA,DIA(ISTIA+NIA*NIO+NIO),NIA,One,HC(ISTBM+NNST),NIA)
    end if

98  ISTIA = ISTIA+NIA**2
    ISTAE = ISTAE+NAE**2
    ISTBM = ISTBM+NIA*NAE
    ISTH = ISTH+NAO*NAE
    ISTZ = ISTZ+(NAO**2-NAO)/2
  end do

  ! Add diagonal contributions to the CI part

  if (ICICP /= 0) then
    do I=1,NROOT
      HC(NST+I) = HC(NST+I)+CIN(NST+I)*(ENER(I,ITER)-ENER(IROOT(1),ITER))
    end do
  end if

  if (NSXS /= 0) then

    ! BM contributions to CI part of sigma

    !call DGEMTX(NSXS,NROOT,One,BM,NSXS,C,1,HC(NST+1),1)
    call DGEMV_('T',NSXS,NROOT,One,BM,NSXS,C,1,One,HC(NST+1),1)

    ! BM contributions to SX part of sigma

    call DAXPY_(NSXS,CIN(NST+1),BM,1,HC(NNST+1),1)

  end if

  ! Remove any unwanted rotations:
  do I=1,NSXS
    if (HD(I+NROOT) > 1.0e20_wp) HC(I+NNST) = Zero
  end do

  ! Adding a constant times C to the HC vectors at this point
  ! is equivalent to modifying the underlying minimization problem
  ! by adding a penalty for orbital rotations.
  ! This should be large enough that rotations get confined to values
  ! for which the Taylor expansion of rotations (which underlies the
  ! orbital optimization theory) is not too unrealistic.
  ! Example: if SXDAMP=0.0001, then a total rotation 'angle' of 0.1
  ! radians is penalized as an energy increase of 0.0001*(0.1)**2
  ! i.e., 10**(-6) a.u.
  ! Larger SXDAMP may make iterations sluggier.
  ! Too low SXDAMP may give erratic convergence in some exceptional
  ! cases.
  ! SXDAMP is set in READIN.
  call DAXPY_(NSXS,SXDAMP,C,1,HC(NNST+1),1)

  ! Renormalize the sigma vector

  do I=1,NSXS
    HC(I+NNST) = HC(I+NNST)*SXN(I)
  end do

  ! Add Level shift part of the diagonal

  call DAXPY_(NSXS,SXSHFT,CIN(1+NNST),1,HC(1+NNST),1)
  NST = NST+NDIMSX
end do

! Test print out of the sigma vector

if (IPRLEV >= DEBUG) write(u6,1000) (HC(I),I=1,NDIMSX)
1000 format(/1X,'Sigma vector in SIGVEC'/(1X,10F11.6))

return

end subroutine SIGVEC
