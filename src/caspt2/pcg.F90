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
! Copyright (C) 1996,1999, Per Ake Malmqvist                           *
!***********************************************************************
!--------------------------------------------*
! 1996  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
! 1999: GEMINAL-R12 ENABLED                  *
!--------------------------------------------*

subroutine PCG(ICONV)

use INPUTDATA, only: INPUT
use caspt2_global, only: EMP2
use caspt2_global, only: iPrGlb
use caspt2_global, only: sigma_p_epsilon, imag_shift, real_shift
use caspt2_global, only: do_grad, nStpGrd
use caspt2_global, only: LISTS
use PrintLevel, only: TERSE, USUAL
use stdalloc, only: mma_allocate, mma_deallocate
use EQSOLV, only: iRHS, iVecc, iVecc2, iVecR, iVecX, NLSTOT
use caspt2_module, only: DeNorm, E2Corr, E2Tot, MxCase, ERef, IfChol, MaxIt, nSym, rNorm, ThrConv, Cases
use constants, only: Zero, One, Two
use definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp), intent(out) :: ICONV
integer(kind=iwp) IC, IS, ITER
integer(kind=iwp) IVECP, IVECT, IVECU
integer(kind=iwp) LAXITY
integer(kind=iwp), external :: Cho_X_GetTol
real(kind=wp) ALPHA, BETA, PR, PT, UR
real(kind=wp) ECORR(0:8,0:MXCASE)
real(kind=wp) EAIVX, EATVX, EBJAI, EBJAT, EBVAT, EVJAI, EVJTI, EVJTU
real(kind=wp) E2NONV, ESHIFT
real(kind=wp) OVLAPS(0:8,0:MXCASE)
real(kind=wp) SAV, SAVI, savreg, DSCALE, REFWGT

! Flag to tell whether convergence was obtained
ICONV = 0

! Lists of coupling coefficients, used for sigma vector
! generation from non-diagonal blocks of H0.
call mma_allocate(LISTS,NLSTOT,LABEL='LISTS')
call MKLIST(LISTS,NLSTOT)

! Mnemonic names for vectors stored on LUSOLV, see EQCTL.
! Here, we use the local names IVECP, IVECT, IVECU which are thus
! to be seen as overlayed areas. The true vectors IVECC and IVECC2
! are computed on return from this routine, so for a while we use them
! for temporaries.
IVECP = IVECC
IVECT = IVECC2
IVECU = IVECC2

! Solve equations for the diagonal case, in eigenbasis:
! Current solution vector X, Current residual vector R
call PSCAVEC(-One,IRHS,IVECR)
call PRESDIA(IVECR,IVECX,OVLAPS)

if (MAXIT == 0) then
  if (IPRGLB >= TERSE) then
    write(u6,*)
    write(u6,'(A)') repeat('-',115)
    write(u6,*) ' DIAGONAL CASPT2 APPROXIMATION:'
  end if
else

  RNORM = Zero
  ! Pre-conditioned conjugate gradient:
  ! R <- R - (H0-E0)*X
  call SIGMA_CASPT2(-One,One,IVECX,IVECR)
  call POVLVEC(IVECR,IVECR,OVLAPS)
  RNORM = sqrt(OVLAPS(0,0))

  if (RNORM >= THRCONV) then

    if (IPRGLB >= USUAL) then
      write(u6,*)
      write(u6,*) 'The contributions to the second order correlation energy in atomic units.'
      write(u6,'(A)') repeat('-',125)
      write(u6,'(2X,A)') 'IT.      VJTU        VJTI        ATVX        AIVX        VJAI        BVAT        BJAT        BJAI'// &
                         '        TOTAL       RNORM  '
      write(u6,'(A)') repeat('-',125)
    end if

    call PRESDIA(IVECR,IVECP,OVLAPS)
    ! PCG iteration loops:
    !---------------------
    ITER = 0
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

      if (RNORM < THRCONV) exit

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
      if (IPRGLB >= USUAL) write(u6,'(1X,I3,1X,10F12.6)') ITER,EVJTU,EVJTI,EATVX,EAIVX,EVJAI,EBVAT,EBJAT,EBJAI,E2NONV,RNORM

      if (ITER >= MAXIT) then
        if (IPRGLB >= TERSE) then
          write(u6,*)
          write(u6,*) ' NOT CONVERGED AFTER MAX ITERATIONS.'
        end if
        ICONV = 16
        exit
      end if

      call PRESDIA(IVECR,IVECU,OVLAPS)
      UR = OVLAPS(0,0)
      BETA = PR/UR
      call PLCVEC(BETA,One,IVECU,IVECP)
    end do
    !---------------------
  end if

end if

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
call POVLVEC(IVECX,IVECX,OVLAPS)
DENORM = One+OVLAPS(0,0)
REFWGT = One/DENORM
!PAM Insert: Compute the variational second-order energy.
!PAM Use unshifted H0. Save any shifts, then restore them.
SAV = real_shift
SAVI = imag_shift
savreg = sigma_p_epsilon
real_shift = Zero
imag_shift = Zero
sigma_p_epsilon = Zero
call SIGMA_CASPT2(One,Zero,IVECX,IVECT)
real_shift = SAV
imag_shift = SAVI
sigma_p_epsilon = savreg
call POVLVEC(IVECX,IVECT,OVLAPS)
E2CORR = Two*E2NONV+OVLAPS(0,0)
!PAM End of insert.
ESHIFT = E2CORR-E2NONV
E2TOT = EREF+E2CORR

if (IPRGLB > USUAL) then
  write(u6,*)
  write(u6,*) ' Correlation energy /Case, /Symm, and sums:'
  do IC=1,13
    write(u6,'(1X,A8,9F12.8)') CASES(IC),(ECORR(IS,IC),IS=1,NSYM),ECORR(0,IC)
  end do
  write(u6,'(1X,A8,9F12.8)') 'Summed: ',(ECORR(IS,0),IS=1,NSYM),ECORR(0,0)
end if

if ((nStpGrd == 1) .or. ((nStpGrd == 2) .and. (.not. do_grad))) then
  if (IPRGLB >= TERSE) then
    write(u6,'(A)') repeat('-',125)
    write(u6,*)
    write(u6,*) ' FINAL CASPT2 RESULT:'
    write(u6,*)

    if (.not. Input%LovCASPT2) then
      write(u6,'(6x,a,f18.10)') 'Reference energy:     ',EREF
      write(u6,'(6x,a,f18.10)') 'E2 (Non-variational): ',E2NONV
      if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero)) &
        write(u6,'(6x,a,f18.10)') 'Shift correction:     ',ESHIFT
      write(u6,'(6x,a,f18.10)') 'E2 (Variational):     ',E2CORR
      if (.not. Input%FnoCASPT2) then
        write(u6,'(6x,a,f18.10)') 'Total energy:         ',E2TOT
      else
        write(u6,'(6x,a,f18.10,a)') 'FNO correction:       ',EMP2,'   (estimate)   '
        write(u6,'(6x,a,f13.5)')
        E2TOT = E2TOT+EMP2
        write(u6,'(6x,a,f18.10,a)') 'Total energy:         ',E2TOT,'   (FNO-CASPT2) '
      end if
      write(u6,'(6x,a,f18.10)') 'Residual norm:        ',RNORM
      write(u6,'(6x,a,f13.5)') 'Reference weight:     ',REFWGT
    else
      write(u6,'(6x,a,f18.10)') 'Reference energy:                 ',EREF
      write(u6,'(6x,a,f18.10)') 'Active-Site E2 (Non-variational): ',E2NONV
      if ((real_shift /= Zero) .or. (imag_shift /= Zero)) write(u6,'(6x,a,f18.10)') 'Shift correction:                 ',ESHIFT
      write(u6,'(6x,a,f18.10)') 'Active-Site E2 (Variational):     ',E2CORR
      write(u6,'(6x,a,f18.10)') 'Frozen region E2 :                ',EMP2
      write(u6,'(6x,a,f18.10)') 'Residual norm:                    ',RNORM
      write(u6,'(6x,a,f13.5)') 'Reference weight:                 ',REFWGT
      write(u6,'(6x,a,f13.5)')
      E2TOT = E2TOT+EMP2
      write(u6,'(6x,a,f18.10)') 'Total energy (LovCASPT2):         ',E2TOT
    end if
  end if
end if

! In automatic verification calculations, the precision is lower
! in case of Cholesky calculation.
LAXITY = 8
if (IfChol) LAXITY = Cho_X_GetTol(LAXITY)
call Add_Info('E_CASPT2',[E2TOT],1,LAXITY)

if ((nStpGrd == 1) .or. ((nStpGrd == 2) .and. (.not. do_grad))) then
  if (IPRGLB >= USUAL) then
    write(u6,*)
    write(u6,'(6x,a)') 'Contributions to the CASPT2 correlation energy'
    write(u6,'(6x,a,F18.10)') 'Active & Virtual Only:    ',EATVX+EBVAT
    write(u6,'(6x,a,F18.10)') 'One Inactive Excited:     ',EVJTU+EAIVX+EBJAT
    write(u6,'(6x,a,F18.10)') 'Two Inactive Excited:     ',EVJTI+EVJAI+EBJAI
    write(u6,*)
  end if
end if
call mma_deallocate(LISTS)

end subroutine PCG
