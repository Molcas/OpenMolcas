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
! Copyright (C) 2013, Victor P. Vysotskiy                              *
!***********************************************************************

subroutine NIdiag_New(H,U,n,nv,iOpt)
!***********************************************************************
!                                                                      *
! This routine is a wrapper that calls appropriate LAPACK routines to  *
! perform diagonalization of symmetric matrices stored in lower        *
! triangular form.                                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Victor P. Vysotskiy                                         *
!          Lund university, Sweden                                     *
! Written  2013                                                        *
!                                                                      *
!***********************************************************************

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
! n    - Dimension of matrix                                           *
! nv   - Length of eigenvectors nv>=n                                  *
! H    - Matrix to be diagonalized                                     *
! U    - Eigenvectors                                                  *
! iOpt - Option flag, for future improvements.                         *
!----------------------------------------------------------------------*
#include "WrkSpc.fh"
external OrbPhase
integer n, nv, iOpt
real*8 H(*), U(nv,n)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer ipDIA, ipOFF, ipTAU, ipIWRK, ipWork, ipEVL, ipIPSZ, ipHDUP
integer lrwrk, liwrk, lh, info, I, M
real*8 abstol, dlamch_, Tmp, OrbPhase
external dlamch_
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if (n == 0) return
#ifdef _DEBUGPRINT_
write(6,*) 'New nidiag'
#endif
call FZero(U,nv*n)

lh = n*(n+1)/2
liwrk = 10*n
lrwrk = 20*n

call GetMem('DIA','ALLO','REAL',ipDIA,n)
call GetMem('EVL','ALLO','REAL',ipEVL,n)
call GetMem('OFF','ALLO','REAL',ipOFF,(n-1))
call GetMem('TAU','ALLO','REAL',ipTAU,(n-1))
call GetMem('IPSZ','ALLO','INTE',ipIPSZ,2*n)
call GetMem('IWRK','ALLO','INTE',ipIWRK,liwrk)
call GetMem('RWRK','ALLO','REAL',ipWORK,lrwrk)
call GetMem('HDUP','ALLO','REAL',ipHDUP,lh)

call dcopy_(lh,H,1,Work(ipHDUP),1)

info = 0
call dsptrd_('U',n,Work(ipHDUP),Work(ipDIA),Work(ipOFF),Work(ipTAU),info)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(6,'(A,I4)') 'Failed to tridiagonalize matrix',info
# endif
  Go To 10
end if
#if defined(_ACML_) && defined(__PGI)
call ILAENVSET(10,'X','X',0,0,0,0,1,INFO)
call ILAENVSET(11,'X','X',0,0,0,0,1,INFO)
#endif
abstol = dlamch_('Safe minimum')
info = 0
call dstevr_('V','A',n,Work(ipDIA),Work(ipOFF),Work(ip_Dummy),Work(ip_Dummy),iWork(ip_iDummy),iWork(ip_iDummy),abstol,M, &
             Work(ipEVL),U,nv,iWork(ipIPSZ),Work(ipWORK),lrwrk,iWork(ipIWRK),liwrk,info)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(6,'(A,I4)') 'Failed to diagonalize matrix',info
# endif
  Go To 10
end if

call dopmtr_('Left','U','N',N,N,Work(ipHDUP),Work(ipTAU),U,nv,Work(ipWork),info)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(6,'(A,I4)') 'Failed to back transform vectors',info
# endif
  Go To 10
end if

call dcopy_(lh,Work(ipHDUP),1,H,1)

do I=1,N
  H((I*(I+1))/2) = Work(ipEVL+I-1)
end do

10 continue
call GetMem('DIA','FREE','REAL',ipDIA,n)
call GetMem('EVL','FREE','REAL',ipEVL,n)
call GetMem('OFF','FREE','REAL',ipOFF,(n-1))
call GetMem('TAU','FREE','REAL',ipTAU,(n-1))
call GetMem('IPSZ','FREE','INTE',ipIPSZ,2*n)
call GetMem('RWRK','FREE','REAL',ipWORK,lrwrk)
call GetMem('IWRK','FREE','INTE',ipIWRK,liwrk)
call GetMem('HDUP','FREE','REAL',ipHDUP,n*(n+1)/2)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(6,'(A)') 'Using the old Givens rot. based routine'
# endif
  call NIdiag(H,U,n,nv,iOpt)
end if

do i=1,n
  Tmp = OrbPhase(U(1,i),nv)
end do
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_real(Tmp)
#endif

end subroutine NIdiag_New
