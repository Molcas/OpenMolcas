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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine PCG_RES(ICONV)

use PrintLevel, only: TERSE, USUAL
use EQSOLV, only: iRHS, iVecc, iVecc2, iVecR, iVecX
use caspt2_global, only: iPrGlb
use caspt2_module, only: HZERO, MaxIt, MxCase, rNorm, ThrConv
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ICONV
integer(kind=iwp) :: ITER, IVECP, IVECT, IVECU
real(kind=wp) :: ALPHA, BETA, DSCALE, E2NONV, EAIVX, EATVX, EBJAI, EBJAT, EBVAT, ECORR(0:8,0:MXCASE), EVJAI, EVJTI, EVJTU, &
                 OVLAPS(0:8,0:MXCASE), PR, PT, UR
logical(kind=iwp) :: converged

! Flag to tell wether convergence was obtained
ICONV = 0
converged = .false.

! Mnemonic names for vectors stored on LUSOLV, see EQCTL.
! Here, we use the local names IVECP, IVECT, IVECU which are thus
! to be seen as overlayed areas. The true vectors IVECC and IVECC2
! are computed on return from this routine, so for a while we use them
! for temporaries.
IVECP = IVECC
IVECT = IVECC2
IVECU = IVECC2

ITER = 0
RNORM = Zero

! Solve equations for the diagonal case, in eigenbasis:
! Current solution vector X, Current residual vector R
call PSCAVEC(-One,IRHS,IVECR)
call PRESDIA(IVECR,IVECX,OVLAPS)
if (MAXIT == 0) then
  if (HZERO /= 'DYALL') then
    if (IPRGLB >= TERSE) then
      write(u6,*)
      write(u6,'(A)') repeat('-',115)
      write(u6,*) ' DIAGONAL CASPT2 APPROXIMATION:'
    end if
  end if
  converged = .true.
else

  ! Pre-conditioned conjugate gradient:
  ! R <- R - (H0-E0)*X
  call SIGMA_CASPT2(-One,One,IVECX,IVECR)
  call POVLVEC(IVECR,IVECR,OVLAPS)
  RNORM = sqrt(OVLAPS(0,0))
  if (RNORM < THRCONV) then
    converged = .true.
  else

    if (IPRGLB >= USUAL) then
      write(u6,*)
      write(u6,*) 'Solving the Lambda equation for analytic gradients'
      write(u6,'(A)') repeat('-',125)
      write(u6,'(2X,A)') 'IT.      VJTU        VJTI        ATVX        AIVX        VJAI        BVAT        BJAT        BJAI'// &
                         '        TOTAL       RNORM  '
      write(u6,'(A)') repeat('-',125)
    end if
    call PRESDIA(IVECR,IVECP,OVLAPS)
    ! PCG iteration loops:
    !---------------------
    do
      call POVLVEC(IVECP,IVECP,OVLAPS)
      DSCALE = One/sqrt(OVLAPS(0,0))
      call PSCAVEC(DSCALE,IVECP,IVECP)
      call POVLVEC(IVECP,IVECR,OVLAPS)
      PR = OVLAPS(0,0)
      call SIGMA_CASPT2(One,Zero,IVECP,IVECT)
      call POVLVEC(IVECP,IVECT,OVLAPS)
      PT = OVLAPS(0,0)
      ALPHA = PR/PT
      call PLCVEC(ALPHA,One,IVECP,IVECX)
      call PLCVEC(-ALPHA,One,IVECT,IVECR)
      call POVLVEC(IVECR,IVECR,OVLAPS)
      RNORM = sqrt(OVLAPS(0,0))
      if (RNORM < THRCONV) then
        converged = .true.
        exit
      end if
      ITER = ITER+1
      call POVLVEC(IRHS,IVECX,ECORR)
      EVJTU = ECORR(0,1)
      EVJTI = ECORR(0,2)+ECORR(0,3)
      EATVX = ECORR(0,4)
      EAIVX = ECORR(0,5)
      EVJAI = ECORR(0,6)+ECORR(0,7)
      EBVAT = ECORR(0,8)+ECORR(0,9)
      EBJAT = ECORR(0,10)+ECORR(0,11)
      EBJAI = ECORR(0,12)+ECORR(0,13)
      E2NONV = ECORR(0,0)
      if (IPRGLB >= USUAL) then
        write(u6,'(1X,I3,1X,10F12.6)') ITER,EVJTU,EVJTI,EATVX,EAIVX,EVJAI,EBVAT,EBJAT,EBJAI,E2NONV,RNORM
        call XFLUSH(u6)
      end if
      if (ITER >= MAXIT) exit
      call PRESDIA(IVECR,IVECU,OVLAPS)
      UR = OVLAPS(0,0)
      BETA = PR/UR
      call PLCVEC(BETA,One,IVECU,IVECP)
    end do
  end if
  !---------------------

  if (.not. converged) then
    if (IPRGLB >= TERSE) then
      write(u6,*)
      write(u6,*) ' NOT CONVERGED AFTER MAX ITERATIONS.'
    end if
    ICONV = 16
  end if
end if

if (IPRGLB >= TERSE) then
  write(u6,'(A)') repeat('-',125)
  write(u6,*)
end if

end subroutine PCG_RES
