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
! Copyright (C) 2005, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2005  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine EQCTL2(ICONV)
! On return, the following data sets will be defined and stored
! on LUSOLV.
! At position IVEC=IRHS, the RHS array, in SR representation.
! At position IVEC=IVECX, the solution array, in SR representation.
! At position IVEC=IVECR, the residual array, in SR representation.
! At position IVEC=IVECC, the solution array, in contravariant rep.
! At position IVEC=IVECC2, the solution array, in covariant repr.
! At position IVEC=IVECW, the RHS array, in contravariant repr.

use PrintLevel, only: INSANE, USUAL, VERBOSE
use EQSOLV, only: IRHS, IVECC, IVECC2, IVECR, IVECW, IVECX
use ChoCASPT2, only: iALGO
use caspt2_global, only: do_grad, iPrGlb, iStpGrd, nStpGrd
use caspt2_module, only: CPUEIG, CPULCS, CPUNAD, CPUOVL, CPUPCG, CPURHS, CPUSBM, CPUSCA, CPUSER, CPUSGM, CPUVEC, E2TOT, HZERO, &
                         IfChol, NASUP, NINDEP, NISUP, NSYM, RHSDIRECT, SDECOM, SMATRIX, TIOEIG, TIOLCS, TIONAD, TIOOVL, TIOPCG, &
                         TIORHS, TIOSBM, TIOSCA, TIOSER, TIOSGM, TIOVEC
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: ICONV
integer(kind=iwp) :: ICASE, ISYM, LAXITY
real(kind=wp) :: CPU, CPU0, CPU1, TIO, TIO0, TIO1
integer(kind=iwp), external :: Cho_X_GetTol

if (iStpGrd == 1) then

  if (IPRGLB >= VERBOSE) then
    write(u6,'(1X,A)')
    write(u6,'(1X,A)') 'Computing the S/B matrices'
    write(u6,'(1X,A)') '--------------------------'
  end if

  ! Compute S and (possibly) B matrices and write them on LUSBT
  ! this uses previously stored data on LUSOLV!!
  call GASync()
  call TIMING(CPU0,CPU,TIO0,TIO)

  ! Necessary to reset NINDEP to conservative estimate. It gets its
  ! final value after SBDIAG call, but must have its original value when
  ! calling MKSMAT and MKBMAT.
  do ICASE=1,13
    do ISYM=1,NSYM
      if (NISUP(ISYM,ICASE) == 0) then
        NINDEP(ISYM,ICASE) = 0
      else
        NINDEP(ISYM,ICASE) = NASUP(ISYM,ICASE)
      end if
    end do
  end do

  if (SMATRIX /= 'NO') call MKSBMAT()

  ! Modify B matrices, if necessary:
  if (HZERO == 'CUSTOM') call NEWB()

  call GASync()
  call TIMING(CPU1,CPU,TIO1,TIO)
  CPUSBM = CPU1-CPU0
  TIOSBM = TIO1-TIO0

  ! Linear dependence removal, ON transformation of S matrix,
  ! and spectral resolution of H0:
  call GASync()
  call TIMING(CPU0,CPU,TIO0,TIO)

  if (SDECOM /= 'NO') call SBDIAG()

  call GASync()
  call TIMING(CPU1,CPU,TIO1,TIO)
  CPUEIG = CPU1-CPU0
  TIOEIG = TIO1-TIO0
end if

! The transformation matrices have now been computed and
! written at disk addresses given by IDTMAT().
! The B matrices were destroyed here. In their place,
! at the IDBMAT() addresses, we find two sets of diagonal
! energy values. Each consists of: first the energies of
! the active combination, which are eigenvalues of the
! B matrix using overlap S; then the energies of the
! non-active combinations. Usually, only the first set of
! energies is used. The second is used, in a manner specific
! to future modifications, for additional energy parameters
! of more sophisticated H0 or resolvent definitions.
! Modif PAM 961022: If BMATRIX='YES     ', diagonal energies
! are put at IDBMAT(). If BTRANS='YES     ', these are the
! diagonal values of B after transformation to ON basis.
! If furthermore BSPECT='YES     ', the ON basis consists
! of eigenvectors. This is the usual CASPT2 situation.
! However, if only BMATRIX is 'YES     ', then the values
! are the diagonal values of B divided by diagonal values
! of S.

call GASync()
call TIMING(CPU0,CPU,TIO0,TIO)

! Non-active part of diagonal elements of H0 are computed
! and written to LUSBT:
call NADIAG()
! Modify diagonal elements, if requested:
if (HZERO == 'CUSTOM') call NEWDIA()

! A second set of energy parameters may now have been
! computed and written to LUSBT.
call GASync()
call TIMING(CPU1,CPU,TIO1,TIO)
CPUNAD = CPU1-CPU0
TIONAD = TIO1-TIO0

if (IPRGLB >= VERBOSE) then
  write(u6,'(1X,A)')
  write(u6,'(1X,A)') 'Computing the right-hand side (RHS) elements'
  write(u6,'(1X,A)') '--------------------------------------------'
end if

IRHS = 1
IVECX = 2
IVECR = 3
IVECC = 4
IVECC2 = 5
IVECW = 6

call GASync()
call TIMING(CPU0,CPU,TIO0,TIO)

!-SVC: initialize the RHS array offsets
call RHS_INIT()
! at this point LUSOLV is not used as a temporary disk anymore, so it's
! safe to initialize it (safe to write zeros to LUSOLV or delete it)

! Set up right-hand side matrix elements.
if (IfChol .and. (iALGO == 1)) then
  if (RHSDIRECT) then
    if (NSYM == 1) then
      call RHSOD_NOSYM(IVECW)
    else
      call RHSOD(IVECW)
    end if
  else
    call RHS_ZERO(IVECW)
    call RHSALL2(IVECW)
  end if
else
  call MKRHS(IVECW)
end if

call GASync()
call TIMING(CPU1,CPU,TIO1,TIO)
CPURHS = CPU1-CPU0
TIORHS = TIO1-TIO0

if (IPRGLB >= INSANE) then
  write(u6,'("DEBUG> ")')
  write(u6,'("DEBUG> ",A)') 'Norms of the RHS blocks:'
  call RHS_FPRINT('C',IVECW)
end if

!-SVC: start PCG routine, set timers.
call GASync()
call TIMING(CPU0,CPU,TIO0,TIO)
CPUSCA = 0
CPULCS = 0
CPUOVL = 0
CPUSGM = 0
CPUVEC = 0
TIOSCA = 0
TIOLCS = 0
TIOOVL = 0
TIOSGM = 0
TIOVEC = 0

! Transform RHS of CASPT2 equations to eigenbasis for H0:
call PTRTOSR(1,IVECW,IRHS)

if (IPRGLB >= INSANE) then
  write(u6,'("DEBUG> ")')
  write(u6,'("DEBUG> ",A)') 'Norms of the RHS blocks (H0 eigenbasis):'
  call RHS_FPRINT('SR',IRHS)
end if

if (iStpGrd == 1) call PCG(ICONV)
if (iStpGrd == 2) then
  call RHS_ZERO(IVECR)
  ICONV = 0
  !! Just for verification
  LAXITY = 8
  if (IfChol) LAXITY = Cho_X_GetTol(LAXITY)
  call Add_Info('E_CASPT2',[E2TOT],1,LAXITY)
end if

call PTRTOC(0,IVECX,IVECC)
call PTRTOC(1,IVECX,IVECC2)

!-SVC: end of PCG routine, compute total time.
call GASync()
call TIMING(CPU1,CPU,TIO1,TIO)
CPUPCG = CPU1-CPU0
TIOPCG = TIO1-TIO0

!-SVC: collect and print information on coefficients/denominators
if ((nStpGrd == 1) .or. ((nStpGrd == 2) .and. (.not. do_grad))) then
  if (IPRGLB >= USUAL) call H0SPCT()
end if

call GASync()
call TIMING(CPU0,CPU,TIO0,TIO)
call GASync()
call TIMING(CPU1,CPU,TIO1,TIO)
CPUSER = CPU1-CPU0
TIOSER = TIO1-TIO0

end subroutine EQCTL2
