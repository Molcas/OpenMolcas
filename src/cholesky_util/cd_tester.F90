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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  CD_Tester
!
!> @brief
!>   Test the decomposition modules
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Runs and tests the output from the general Cholesky decomposers
!> ::ChoDec (out-of-core) and ::CD_InCore (in-core). The positive
!> definite matrix to which \p ip_PDM points should be stored as a
!> lower triangle.
!>
!> Error codes:
!> - \p irc = ``-1``: \p n non-positive (&rarr; nothing done)
!> - \p irc =  ``0``: all OK
!> - \p irc =  ``1``: error in ::ChoDec
!> - \p irc =  ``2``: error in ::CD_InCore
!> - \p irc =  ``3``: error in both
!>
!> @param[out] irc     Return code
!> @param[in]  PDM     matrix in \p Work
!> @param[in]  n       Dimension of matrix (\p n &times; \p n)
!> @param[in]  Verbose Print flag
!***********************************************************************

subroutine CD_Tester(irc,PDM,n,Verbose)

use CDTHLP
use stdalloc

implicit real*8(a-h,o-z)
integer irc, n
real*8 PDM(n*(n+1)/2)
external CD_Tester_Col
external CD_Tester_Vec
logical Verbose
character(len=9), parameter :: SecNam = 'CD_Tester'
logical Restart
real*8, allocatable :: Diag(:), Buf(:), Qual(:), ES(:)
integer, allocatable :: Pivot(:), iQual(:)

irc = 0
if (n < 1) then
  if (Verbose) then
    write(6,*) SecNam,': nothing tested! Dimension is: n = ',n
    call xFlush(6)
  end if
  Go To 1
end if

! Test ChoDec.
! ============

Restart = .false.
Thr = 1.0d-12
Span = 1.0d-2
MxQual = max(n/10,1)
NumCho = 0

irc_sav = irc
if (Verbose) then
  write(6,*)
  write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<'
  write(6,*) '          >>>> Testing ChoDec <<<<'
  write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<'
  write(6,*)
  call xFlush(6)
end if

call mma_allocate(Mat,n*n,Label='Mat')
call mma_allocate(Vec,n*n,Label='Vec')

l_Buf = n*(n+MxQual)
l_Qual = n*(MxQual+1)
l_ES = 6
l_Pivot = n
l_iQual = MxQual

call mma_allocate(Buf,l_Buf,Label='Buf')
call mma_allocate(Diag,n,Label='Diag')
call mma_allocate(Qual,l_Qual,Label='Qual')
call mma_allocate(ES,l_ES,Label='ES')
call mma_allocate(Pivot,l_Pivot,Label='Pivot')
call mma_allocate(iQual,l_iQual,Label='iQual')

call CD_Tester_CPPF(PDM,Mat,n)
call CD_Tester_Diag(PDM,Diag,n)
jrc = 0
call ChoDec(CD_Tester_Col,CD_Tester_Vec,Restart,Thr,Span,MxQual,Diag,Qual,Buf,Pivot,iQual,n,l_Buf,ES,NumCho,jrc)
if (jrc == 0) then
  call DGEMM_('N','T',n,n,NumCho,-1.0d0,Vec,n,Vec,n,1.0d0,Mat,n)
  call CD_Tester_ES(Mat,n,ES)
  call CD_Tester_Diff(Mat,n,ES(4))
  call CD_Tester_Final(jrc,NumCho,n,Thr,ES,Verbose)
  if (jrc /= 0) irc = irc+1
else
  if (Verbose) then
    write(6,*) SecNam,': ChoDec returned error code ',jrc
    write(6,*) '<=> decomposition failed, no test performed!'
    call xFlush(6)
  end if
  irc = irc+1
end if
if (Verbose) then
  if (irc == irc_sav) then
    write(6,*) 'ChoDec succeeded....'
  else
    write(6,*) 'ChoDec failed....'
  end if
  call xFlush(6)
end if

! Test in-core decomposition.
! ===========================

if (Verbose) then
  irc_sav = irc
  write(6,*)
  write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<<<<'
  write(6,*) '          >>>> Testing CD_InCore <<<<'
  write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<<<<'
  write(6,*)
  call xFlush(6)
end if

call CD_Tester_CPPF(PDM,Mat,n)
jrc = 0
NumCho = 0
call CD_InCore(Mat,n,Vec,n,NumCho,Thr,jrc)
if (jrc == 0) then
  call CD_Tester_CPPF(PDM,Mat,n)
  call DGEMM_('N','T',n,n,NumCho,-1.0d0,Vec,n,Vec,n,1.0d0,Mat,n)
  call CD_Tester_ES(Mat,n,ES)
  call CD_Tester_Diff(Mat,n,ES(4))
  call CD_Tester_Final(jrc,NumCho,n,Thr,ES,Verbose)
  if (jrc /= 0) irc = irc+2
else
  if (Verbose) then
    write(6,*) SecNam,': CD_InCore returned error code ',jrc
    write(6,*) '<=> decomposition failed, no test performed!'
    call xFlush(6)
  end if
  irc = irc+2
end if
if (Verbose) then
  if (irc == irc_sav) then
    write(6,*) 'CD_InCore succeeded....'
  else
    write(6,*) 'CD_InCore failed....'
  end if
  write(6,*)
  call xFlush(6)
end if

call mma_deallocate(iQual)
call mma_deallocate(Pivot)
call mma_deallocate(ES)
call mma_deallocate(Qual)
call mma_deallocate(Diag)
call mma_deallocate(Buf)
call mma_deallocate(Vec)
call mma_deallocate(Mat)

1 continue

end subroutine CD_Tester
