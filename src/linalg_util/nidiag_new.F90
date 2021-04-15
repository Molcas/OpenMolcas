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

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, r8
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
! n    - Dimension of matrix                                           *
! nv   - Length of eigenvectors nv>=n                                  *
! H    - Matrix to be diagonalized                                     *
! U    - Eigenvectors                                                  *
! iOpt - Option flag, for future improvements.                         *
!----------------------------------------------------------------------*
integer(kind=iwp), intent(in) :: n, nv, iOpt
real(kind=wp), intent(inout) :: H(*)
real(kind=wp), intent(out) :: U(nv,n)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer(kind=iwp) :: lrwrk, liwrk, lh, info, I, M, idum(1)
real(kind=wp) :: abstol, dum(1), Tmp
integer(kind=iwp), allocatable :: IPSZ(:), IWRK(:)
real(kind=wp), allocatable :: DIA(:), EVL(:), HDUP(:), OFF(:), RWRK(:), TAU(:)
real(kind=r8), external :: dlamch_, OrbPhase
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if (n == 0) return
#ifdef _DEBUGPRINT_
write(u6,*) 'New nidiag'
#endif
call FZero(U,nv*n)

lh = n*(n+1)/2
liwrk = 10*n
lrwrk = 20*n

call mma_allocate(DIA,n,label='DIA')
call mma_allocate(EVL,n,label='EVL')
call mma_allocate(OFF,n-1,label='OFF')
call mma_allocate(TAU,n-1,label='TAU')
call mma_allocate(IPSZ,2*n,label='IPSZ')
call mma_allocate(IWRK,liwrk,label='IWRK')
call mma_allocate(RWRK,lrwrk,label='RWRK')
call mma_allocate(HDUP,lh,label='HDUP')

call dcopy_(lh,H,1,HDUP,1)

info = 0
call dsptrd_('U',n,HDUP,DIA,OFF,TAU,info)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(u6,'(A,I4)') 'Failed to tridiagonalize matrix',info
# endif
  Go To 10
end if
#if defined(_ACML_) && defined(__PGI)
call ILAENVSET(10,'X','X',0,0,0,0,1,INFO)
call ILAENVSET(11,'X','X',0,0,0,0,1,INFO)
#endif
abstol = dlamch_('Safe minimum')
info = 0
call dstevr_('V','A',n,DIA,OFF,dum,dum,idum,idum,abstol,M,EVL,U,nv,IPSZ,RWRK,lrwrk,IWRK,liwrk,info)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(u6,'(A,I4)') 'Failed to diagonalize matrix',info
# endif
  Go To 10
end if

call dopmtr_('Left','U','N',N,N,HDUP,TAU,U,nv,RWRK,info)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(u6,'(A,I4)') 'Failed to back transform vectors',info
# endif
  Go To 10
end if

call dcopy_(lh,HDUP,1,H,1)

do I=1,N
  H((I*(I+1))/2) = EVL(I)
end do

10 continue
call mma_deallocate(DIA)
call mma_deallocate(EVL)
call mma_deallocate(OFF)
call mma_deallocate(TAU)
call mma_deallocate(IPSZ)
call mma_deallocate(IWRK)
call mma_deallocate(RWRK)
call mma_deallocate(HDUP)

if (info /= 0) then
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'Using the old Givens rot. based routine'
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
