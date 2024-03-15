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

subroutine CHO_PRTHEAD(SKIPH)
!
! Purpose: print Cholesky header.

use Cholesky, only: BlockSize, Cho_1Center, Cho_AdrVec, Cho_DecAlg, Cho_IOVec, Cho_No2Center, Cho_PreScreen, Cho_SimRI, &
                    Cho_SScreen, Cho_UseAbs, Damp, Frac_ChvBuf, iAlQua, INF_INIT, IPRINT, LuPri, MaxQual, MaxRed, MaxVec, MinQual, &
                    MxShPR, N1_Qual, N2_Qual, RstCho, RstDia, ScDiag, Span, SSTau, Thr_PreScreen, Thr_SimRI, ThrCom, ThrDiag, &
                    ThrNeg, TooNeg, WarNeg
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: SKIPH
integer(kind=iwp) :: I, IADRMODE, IALG, IUSE
real(kind=wp) :: X1, X2, XF, XF2
integer(kind=iwp), parameter :: NADRMODE = 2, NALG = 6
character(len=*), parameter :: ADRMODE(0:NADRMODE) = ['      unknown','   word addr.','  direct acc.'], &
                               ALGORITHM(0:NALG) = ['     unknown','    one-step','    two-step','       naive','par one-step', &
                                                    'par two-step','   par naive'], &
                               USED(2) = ['(screening off)','(screening on) ']

if (LUPRI < 1) call CHO_QUIT('LUPRI undefined in Cholesky decomposition',101)

if (.not. SKIPH) then
  write(LUPRI,'(//,80A)') ('*',I=1,80)
  write(LUPRI,'(A,78X,A)') ('*',I=1,2)
  write(LUPRI,'(A,10X,A,10X,A)') '*','Cholesky Decomposition of Two-Electron Repulsion Integrals','*'
  write(LUPRI,'(A,78X,A)') ('*',I=1,2)
  write(LUPRI,'(80A)') ('*',I=1,80)
  write(LUPRI,*)
  write(LUPRI,*)

  if (RSTDIA) write(LUPRI,'(/,A)') '***** Using Restart Integral Diagonal *****'
  if (RSTCHO) then
    if (RSTDIA) then
      write(LUPRI,'(A)') '***** Using Restart Cholesky Vectors  *****'
    else
      write(LUPRI,'(/,A)') '***** Using Restart Cholesky Vectors  *****'
    end if
  end if
end if

if (IPRINT >= INF_INIT) then
  if (SCDIAG) then
    IUSE = 2
  else
    IUSE = 1
  end if
  if ((CHO_DECALG < 1) .or. (CHO_DECALG > NALG)) then
    IALG = 0
  else
    IALG = CHO_DECALG
  end if
  if (.not. SKIPH) call CHO_HEAD('Configuration','=',80,LUPRI)
  write(LUPRI,'(A,A)') 'Decomposition algorithm                   : ',ALGORITHM(IALG)
  if (CHO_1CENTER) then
    write(LUPRI,'(A)') '1-center decomposition                    :          Yes'
    if (CHO_NO2CENTER) then
      write(LUPRI,'(A)') 'Exclusion of 2-center diagonals           :          Yes'
    else
      write(LUPRI,'(A)') 'Exclusion of 2-center diagonals           :           No'
    end if
    if (CHO_SIMRI) write(LUPRI,'(A,ES12.4)') 'Simulation of RI, threshold               : ',THR_SIMRI
  else
    write(LUPRI,'(A)') '1-center decomposition                    :           No'
  end if
  write(LUPRI,'(A,ES12.4)') 'Decomposition threshold                   : ',THRCOM
  if (CHO_PRESCREEN) write(LUPRI,'(A,ES12.4)') 'Initial diagonal prescreening             : ',THR_PRESCREEN
  write(LUPRI,'(A,ES12.4)') 'Initial diagonal screening                : ',THRDIAG
  write(LUPRI,'(A,ES12.4,1X,A)') 'First  screening damping                  : ',DAMP(1),USED(IUSE)
  write(LUPRI,'(A,ES12.4,1X,A)') 'Second screening damping                  : ',DAMP(2),USED(IUSE)
  if (CHO_USEABS) then
    write(LUPRI,'(A)') 'Absolute values used in diagonal screening:          Yes'
  else
    write(LUPRI,'(A)') 'Absolute values used in diagonal screening:           No'
  end if
  write(LUPRI,'(A,ES12.4)') 'Threshold for negative  diagonal zeroing  : ',THRNEG
  write(LUPRI,'(A,ES12.4)') 'Threshold for warning about neg. diagonal : ',WARNEG
  write(LUPRI,'(A,ES12.4)') 'Threshold for too negative diagonal       : ',TOONEG
  write(LUPRI,'(A,ES12.4)') 'Span factor                               : ',SPAN
  write(LUPRI,'(A,I12)') 'Max. #Cholesky vectors per symmetry       : ',MAXVEC
  write(LUPRI,'(A,I12)') 'Max. #reduced sets (i.e., integral passes): ',MAXRED
  write(LUPRI,'(A,I12)') 'Min. #qualified required for decomposition: ',MINQUAL
  write(LUPRI,'(A,I12)') 'Max. #qualified per symmetry              : ',MAXQUAL
  if (N2_QUAL == 0) then
    XF = -9.99999999e15_wp
  else
    X1 = real(N1_QUAL,kind=wp)
    X2 = real(N2_QUAL,kind=wp)
    XF = 1.0e2_wp*X1/X2
  end if
  write(LUPRI,'(A,5X,F7.4,A)') 'Max. memory fraction for qualified columns: ',XF,'%'
  if (MXSHPR == 0) then
    write(LUPRI,'(A)') 'Max. #shell pair allowed per integral pass:      generic'
  else
    write(LUPRI,'(A,I12)') 'Max. #shell pair allowed per integral pass: ',MXSHPR
  end if
  if (IALQUA == 0) then
    write(LUPRI,'(A)') 'Qualification algorithm                   : dalton-style'
  else if (IALQUA == 1) then
    write(LUPRI,'(A)') 'Qualification algorithm                   :   sequential'
  else
    write(LUPRI,'(A)') 'Qualification algorithm                   :      sorting'
  end if
  if (CHO_IOVEC == 1) then
    write(LUPRI,'(A)') 'Algorithm for Cholesky vector I/O         :  rs2rs/batch'
  else if (CHO_IOVEC == 2) then
    write(LUPRI,'(A)') 'Algorithm for Cholesky vector I/O         : buffer/rs2rs'
  else if (CHO_IOVEC == 3) then
    write(LUPRI,'(A)') 'Algorithm for Cholesky vector I/O         : lrgbuf/rs2rs'
  else if (CHO_IOVEC == 4) then
    write(LUPRI,'(A)') 'Algorithm for Cholesky vector I/O         : fxdbuf/rs2rs'
  else
    write(LUPRI,'(A)') 'Algorithm for Cholesky vector I/O         : copy via rs1'
  end if
  IADRMODE = max(min(CHO_ADRVEC,NADRMODE),0)
  write(LUPRI,'(A,A13)') 'Address mode for Cholesky vector I/O      : ',ADRMODE(IADRMODE)
  XF2 = 1.0e2_wp*FRAC_CHVBUF
  write(LUPRI,'(A,5X,F7.4,A)') 'Memory fraction used as vector buffer     : ',XF2,'%'
  if (CHO_SSCREEN) write(LUPRI,'(A,ES12.4)') 'Screening threshold for vector subtraction: ',SSTAU
  if (CHO_DECALG == 5) write(LUPRI,'(A,I12)') 'Block size (blocked Z vector array)       : ',BLOCKSIZE
end if

end subroutine CHO_PRTHEAD
