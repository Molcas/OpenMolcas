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

subroutine CHO_CHKCONF(NCONFL,VERBOSE)
!
! Purpose: check configuration, return the number of errors NCONFL.

use Index_Functions, only: nTri_Elem
use Cholesky, only: BLOCKSIZE, CHO_1CENTER, CHO_DECALG, CHO_IOVEC, CHO_NDECALG, CHO_NO2CENTER, CHO_PRESCREEN, Cho_Real_Par, &
                    CHO_SIMRI, Cho_SScreen, Damp, FRAC_CHVBUF, IALQUA, IFCSEW, lBuf, LuPri, MaxQual, MaxRed, MaxVec, MinQual, &
                    MXSHPR, N1_Qual, N1_VecRd, N2_Qual, N2_VecRd, N_Subtr, NALQUA, NBAST, nShell, nSym, RstCho, SCDIAG, Span, &
                    SSTau, Thr_PreScreen, ThrCom, ThrDef, ThrDiag, ThrNeg, TOONEG, WARNEG
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: NCONFL
logical(kind=iwp), intent(in) :: VERBOSE
integer(kind=iwp) :: INEGRR, MMM, NNN
logical(kind=iwp) :: REPORT
real(kind=wp) :: XLBUF, XMBUF

! Initialize.
! -----------

NCONFL = 0
REPORT = VERBOSE

! Check that output unit is appropriately set.
! (Upper bound not checked, as it may vary.)
! --------------------------------------------

if (LUPRI < 1) then
  NCONFL = NCONFL+1
  REPORT = .false.
end if

! Check decomposition algorithm.
! ------------------------------

if ((CHO_DECALG < 1) .or. (CHO_DECALG > CHO_NDECALG)) then
  if (REPORT) write(LUPRI,'(A,I4)') 'Illegal decomposition algorithm, CHO_DECALG = ',CHO_DECALG
  NCONFL = NCONFL+1
end if
if (Cho_Real_Par) then
  if ((CHO_DECALG /= 4) .and. (CHO_DECALG /= 5) .and. (CHO_DECALG /= 6)) then
    if (REPORT) then
      write(LUPRI,'(A,I4)') 'Illegal decomposition algorithm, CHO_DECALG = ',CHO_DECALG
      write(LUPRI,'(A)') 'Only parallel algorithm is allowed for parallel execution'
    end if
    NCONFL = NCONFL+1
  end if
end if

! Cancel exclusion of 2-center diagonals in case of symmetry.
! ----------------------------------------------------------

if (CHO_NO2CENTER) then
  if (NSYM /= 1) then
    if (REPORT) write(LUPRI,'(A)') 'Exclusion of 2-center diagonals only implemented for C1 point group.'
    NCONFL = NCONFL+1
  end if
end if

! Cancel 1-center decomposition in case of symmetry.
! --------------------------------------------------

if (CHO_1CENTER) then
  if (NSYM /= 1) then
    if (REPORT) write(LUPRI,'(A)') '1-center decomposition only implemented for C1 point group.'
    NCONFL = NCONFL+1
  end if
end if

! Checks specific to RI simulation.
! ---------------------------------

if (CHO_SIMRI) then
  if (.not. CHO_1CENTER) then
    if (REPORT) write(LUPRI,'(A)') '1-center decomposition required for RI simulation.'
    NCONFL = NCONFL+1
  end if
  if (CHO_DECALG /= 2) then
    if (REPORT) write(LUPRI,'(A)') 'RI simulation can only be executed with the two-step algorithm.'
    NCONFL = NCONFL+1
  end if
  if (RSTCHO) then
    if (REPORT) write(LUPRI,'(A)') 'RI simulation cannot be executed with decomposition restart.'
    NCONFL = NCONFL+1
  end if
end if

! Check for conflicts for specific algorithms.
! --------------------------------------------

if ((CHO_DECALG == 2) .or. (CHO_DECALG == 5)) then
  if (.not. SCDIAG) then
    if (REPORT) write(LUPRI,'(A)') 'Screening must be turned on for two-step decomposition algorithm.'
    NCONFL = NCONFL+1
  end if
  if (IFCSEW /= 2) then
    if (REPORT) then
      write(LUPRI,'(A)') 'The interface to Seward must be "2" (reduced set communication) for two-step algorithm.'
      write(LUPRI,'(A,I4)') 'Current value is IFCSEW=',IFCSEW
    end if
    NCONFL = NCONFL+1
  end if
  if (RSTCHO) then
    if (REPORT) write(LUPRI,'(A)') 'Decomposition restart is not possible for two-step algorithm.'
    NCONFL = NCONFL+1
  end if
end if

if (CHO_DECALG == 5) then
  if (BLOCKSIZE < 1) then
    if (REPORT) write(LUPRI,'(A,I8)') 'BlockSize must be strictly positive.Current value is BlockSize=',BlockSize
    NCONFL = NCONFL+1
  end if
end if

! Check max. number of Cholesky vectors and reduced sets.
! -------------------------------------------------------

if (MAXVEC < 1) then
  if (REPORT) write(LUPRI,'(A,I8)') 'Max. number of vectors < 1: ',MAXVEC
  NCONFL = NCONFL+1
end if
if (MAXRED < 1) then
  if (REPORT) write(LUPRI,'(A,I8)') 'Max. number of reduced sets < 1: ',MAXRED
  NCONFL = NCONFL+1
end if

! Check qualification algorithm.
! ------------------------------

if (IALQUA < 0) then
  if (REPORT) write(LUPRI,'(A,I4,A)') 'Qualification algorithm reset from ',IALQUA,' to 1'
  IALQUA = 1
else if (IALQUA > NALQUA) then
  if (REPORT) write(LUPRI,'(A,I4,A,I4)') 'Qualification algorithm reset from ',IALQUA,' to ',NALQUA
  IALQUA = NALQUA
end if

! Decomposition threshold.
! ------------------------

if (THRCOM < Zero) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6,A,D15.6,A)') 'Decomposition threshold not positive: ',THRCOM,' (default value: ',THRDEF,')'
  NCONFL = NCONFL+1
end if

! Diagonal prescreening threshold.
! --------------------------------

if (CHO_PRESCREEN) then
  if (THR_PRESCREEN > THRCOM) then
    if (REPORT) write(LUPRI,'(A,1P,D15.6,A,D15.6)') 'Diagonal prescreening threshold is greater than decomposition threshold: ', &
                                                    THR_PRESCREEN,' > ',THRCOM
    NCONFL = NCONFL+1
  end if
end if

! Screening threshold for vector subtraction.
! -------------------------------------------

if (CHO_SSCREEN) then
  if (SSTAU < Zero) then
    if (REPORT) write(LUPRI,'(A,1P,D15.6)') 'Screening threshold for vector subtraction not positive: ',SSTAU
    NCONFL = NCONFL+1
  end if
end if

! Buffer length.
! --------------

if (LBUF < 1) then
  if (REPORT) write(LUPRI,'(A,I8)') 'Buffer length < 1: ',LBUF
  NCONFL = NCONFL+1
else
  XLBUF = real(LBUF,kind=wp)
  XMBUF = real(NBAST)*(real(NBAST)+One)*Half
  if (XLBUF > XMBUF) LBUF = nint(XMBUF) ! make sure LBUF is not too large
end if

! Memory split for qualified columns.
! -----------------------------------

if (N1_QUAL >= N2_QUAL) then
  if (REPORT) write(LUPRI,'(A,2I5)') 'N1_QUAL >= N2_QUAL: ',N1_QUAL,N2_QUAL
  NCONFL = NCONFL+1
end if
if ((N1_QUAL < 1) .or. (N2_QUAL < 1)) then
  if (REPORT) write(LUPRI,'(A,2I5)') 'N1_QUAL and/or N2_QUAL < 1: ',N1_QUAL,N2_QUAL
  NCONFL = NCONFL+1
end if

! Memory split for buffered reading of previous vectors.
! Max. #vectors in each call to DGEMM in subtraction.
! (CHO_IOVEC=3,4 only.)
! ------------------------------------------------------

if ((CHO_IOVEC == 3) .or. (CHO_IOVEC == 4)) then
  if (N2_VECRD < N1_VECRD) then
    if (REPORT) write(LUPRI,'(A,2I5)') 'N1_VECRD >= N2_VECRD: ',N1_VECRD,N2_VECRD
    NCONFL = NCONFL+1
  end if
  if ((N1_VECRD < 1) .or. (N2_VECRD < 1)) then
    if (REPORT) write(LUPRI,'(A,2I5)') 'N1_VECRD and/or N2_VECRD < 1: ',N1_VECRD,N2_VECRD
    NCONFL = NCONFL+1
  end if
  if (N_SUBTR < 1) then
    if (REPORT) write(LUPRI,'(A,I8)') 'N_SUBTR: ',N_SUBTR
    NCONFL = NCONFL+1
  end if
end if

! Memory fraction used for vector buffer.
! ---------------------------------------

if (FRAC_CHVBUF < Zero) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6,A)') 'FRAC_CHVBUF=',FRAC_CHVBUF,' resetting value to 0.0'
  FRAC_CHVBUF = Zero
else if (FRAC_CHVBUF > 0.9_wp) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6,A)') 'FRAC_CHVBUF=',FRAC_CHVBUF,' resetting value to 0.9'
  FRAC_CHVBUF = 0.9_wp
end if

! Threshold for discarding elements of initial diagonal.
! ------------------------------------------------------

if (THRDIAG > THRCOM) then
  if (REPORT) then
    write(LUPRI,'(A,A,1P,D15.6)') 'Threshold for discarding initial diagonals',': ',THRDIAG
    write(LUPRI,'(A,1P,D15.6)') 'is larger than decomposition threshold: ',THRCOM
  end if
  NCONFL = NCONFL+1
end if

! Damping factors.
! ----------------

if (DAMP(1) < One) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6)') 'First damping factor < 1: ',DAMP(1)
  NCONFL = NCONFL+1
end if
if (DAMP(2) < One) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6)') 'Second damping factor < 1: ',DAMP(2)
  NCONFL = NCONFL+1
end if

! Span factor.
! ------------

if (SPAN > One) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6)') 'Span factor > 1: ',SPAN
  NCONFL = NCONFL+1
else if (abs(One-SPAN) < 1.0e-4_wp) then
  if (REPORT) write(LUPRI,'(A)') 'Span factor is too close to 1. Will use 0.9999 instead.'
  SPAN = 0.9999_wp
end if

! Max. #shell pairs allowed before proceeding to deco.
! Max. #qualifieds allowed to proceed to decomposition.
! Min. #qualifieds needed to proceed to decomposition.
! -----------------------------------------------------

NNN = nTri_Elem(NSHELL)
if ((MXSHPR < 0) .or. (MXSHPR > NNN)) then
  if (REPORT) then
    write(LUPRI,'(A,I8)') 'Max. #shell pairs allowed before proceeding to deco. is: ',MXSHPR
    write(LUPRI,'(A)') 'Resetting generic: 0'
  end if
  MXSHPR = 0
end if
if (MAXQUAL < 1) then
  if (REPORT) then
    write(LUPRI,'(A,I8)') 'Max. number of qualifieds is non-positive: ',MAXQUAL
    write(LUPRI,'(A,I8)') 'Using instead: ',abs(MAXQUAL)
  end if
  MAXQUAL = abs(MAXQUAL)
end if
MMM = NSYM*MAXQUAL
if ((MINQUAL < 1) .or. (MINQUAL > MMM)) then
  if (REPORT) then
    write(LUPRI,'(A,I8)') 'Min. #qualified needed before proceeding to deco. is: ',MINQUAL
    write(LUPRI,'(A,I8)') 'resetting to: ',MMM
  end if
  MINQUAL = MMM
end if

! Handling of negative diagonals.
! -------------------------------

INEGRR = 0
if (THRNEG > Zero) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6)') 'Threshold for zeroing neg. diag. > 0: ',THRNEG
  INEGRR = INEGRR+1
  NCONFL = NCONFL+1
end if
if (WARNEG > Zero) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6)') 'Threshold for warning about neg. diag.  > 0: ',WARNEG
  INEGRR = INEGRR+1
  NCONFL = NCONFL+1
end if
if (TOONEG > Zero) then
  if (REPORT) write(LUPRI,'(A,1P,D15.6)') 'Threshold for shutdown due to neg. diag.  > 0: ',TOONEG
  INEGRR = INEGRR+1
  NCONFL = NCONFL+1
end if
if (INEGRR == 0) then
  if (THRNEG < WARNEG) then
    if (REPORT) write(LUPRI,'(A,1P,D15.6,D15.6)') 'Threshold for zeroing neg. diag. > threshold for warning about neg. diag.: ', &
                                                  THRNEG,WARNEG
    NCONFL = NCONFL+1
  end if
  if (WARNEG < TOONEG) then
    if (REPORT) write(LUPRI,'(A,1P,D15.6,D15.6)') 'Threshold for warning about neg. diag. > threshold for shutdown due to neg. '// &
                                                  'diag.: ',WARNEG,TOONEG
    NCONFL = NCONFL+1
  end if
end if

end subroutine CHO_CHKCONF
