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

subroutine momloc2(N,NL,nneq,neq,neqv,r_rot,nsites,nexch,W,Z,dipexch,s_exch,dipso,s_so)

use Constants, only: Zero, cZero, cOne, gElectron
use Definitions, only: wp, u6

implicit none
integer :: N
integer :: NL
integer :: nneq
integer :: neq(nneq)
integer :: neqv  ! neqv = MAXVAL(neq(:))
integer :: nsites
integer :: nexch(nneq)
real(kind=8) :: W(N)
! assuming 10 equivalent magnetic sites, which is too much for many cases
real(kind=8) :: R_rot(NNEQ,neqv,3,3)
complex(kind=8) :: dipexch(3,N,N)
complex(kind=8) :: s_exch(3,N,N)
complex(kind=8) :: dipso(nneq,3,NL,NL)
complex(kind=8) :: s_so(nneq,3,NL,NL)
complex(kind=8) :: Z(N,N)
#include "stdalloc.fh"
! local variables:
integer :: L, i, j, m, k
integer :: nmult, isite
integer :: icod(nsites)
integer :: ib(N,nsites)
integer :: nind(nsites,2)
integer :: l_exch, jEnd
integer :: i1, j1, iss1, jss1, nb1, nb2, iss
integer :: icoord(nsites)
integer :: norder
real(kind=8) :: gtens(3)
real(kind=8) :: maxes(3,3)
real(kind=8) :: st(3)
real(kind=8) :: H
real(kind=8) :: E_thres
real(kind=8) :: zJ
character(len=60) :: fmtline
logical :: DBG
real(kind=8), allocatable :: WM(:)     ! WM(N)
real(kind=8), allocatable :: MM(:,:,:) ! MM(nsites,3,N)
real(kind=8), allocatable :: LM(:,:,:) ! LM(nsites,3,N)
real(kind=8), allocatable :: SM(:,:,:) ! SM(nsites,3,N)
real(kind=8), allocatable :: JM(:,:,:) ! JM(nsites,3,N)
complex(kind=8), allocatable :: ZM(:,:)  ! ZM(N,N)
complex(kind=8), allocatable :: VL(:,:)  ! VL(N,N)
complex(kind=8), allocatable :: TMP(:,:) ! TMP(N,N)
! temporary data for ZEEM:
real(kind=8), allocatable :: RWORK(:)
complex(kind=8), allocatable :: HZEE(:), WORK(:), W_c(:)
real(kind=8), parameter :: g_e = -gElectron

DBG = .false.

call mma_allocate(WM,N,'WM')
call mma_allocate(MM,nsites,3,N,'MM')
call mma_allocate(LM,nsites,3,N,'LM')
call mma_allocate(SM,nsites,3,N,'SM')
call mma_allocate(JM,nsites,3,N,'JM')
call mma_allocate(ZM,N,N,'ZM')
call mma_allocate(VL,N,N,'VL')
call mma_allocate(TMP,N,N,'TMP')
! temporary arrays used in ZEEM_SA:
call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')

! zero everything:
call dcopy_(N,[Zero],0,WM,1)
call dcopy_(nsites*3*N,[Zero],0,MM,1)
call dcopy_(nsites*3*N,[Zero],0,LM,1)
call dcopy_(nsites*3*N,[Zero],0,SM,1)
call dcopy_(nsites*3*N,[Zero],0,JM,1)
call zcopy_(N*N,[cZero],0,ZM,1)
call zcopy_(N*N,[cZero],0,VL,1)
call zcopy_(N*N,[cZero],0,TMP,1)

call dcopy_(3*N-2,[Zero],0,RWORK,1)
call zcopy_(N*(N+1)/2,[cZero],0,HZEE,1)
call zcopy_(2*N-1,[cZero],0,WORK,1)
call zcopy_(N,[cZero],0,W_c,1)

if (N == 1) goto 199
! initialisations:
isite = 0
nind = 0
do i=1,nneq
  do j=1,neq(i)
    isite = isite+1
    nind(isite,1) = i
    nind(isite,2) = j
  end do
end do
isite = 0
icod = 0
icod(1) = 1
do i=2,nsites
  isite = nind(i-1,1)
  icod(i) = icod(i-1)*nexch(isite)
end do
ib = 0
do i=1,N
  j = i-1
  do isite=1,nsites
    k = nsites-isite+1
    ib(i,k) = j/icod(k)
    j = j-ib(i,k)*icod(k)
  end do
end do
! find the multiplicity of the low-lying group of states,
! energy threshold E_thres = 1.0e-2 cm-1;
nmult = 0
E_thres = 1.0e-3_wp
do i=1,N
  if (abs(W(i)-W(1)) < E_thres) then
    nmult = nmult+1
  end if
end do
if (nmult < 2) nmult = 2 !minimum value needed to compute g-tensor
! find the main magnetic axes of this manifold:
call atens(dipexch(1:3,1:nmult,1:nmult),nmult,gtens,maxes,1)
if (DBG) then
  write(u6,'(A)') 'MOMLOC2:  g tensor of the ground manifold:'
  write(u6,*)
  write(u6,'((A,F12.6,A,3F12.7))') 'gX=',gtens(1),' axis X: ',(maxes(j,1),j=1,3)
  write(u6,'((A,F12.6,A,3F12.7))') 'gY=',gtens(2),' axis Y: ',(maxes(j,2),j=1,3)
  write(u6,'((A,F12.6,A,3F12.7))') 'gZ=',gtens(3),' axis Z: ',(maxes(j,3),j=1,3)
end if
! construct the Zeeman matrix in the lowest N exchange states  and diagonalize it
st(:) = Zero
H = 1.0e-4_wp ! tesla
zJ = Zero ! absence of intermolecular interaction
call zeem_sa(N,H,maxes(1,3),maxes(2,3),maxes(3,3),W(1:N),dipexch(1:3,1:N,1:N),s_exch(1:3,1:N,1:N),ST,zJ,WM(1:N),ZM(1:N,1:N),DBG, &
             RWORK,HZEE,WORK,W_c)
!cc  eigenvectors
if (DBG) then
  write(u6,*)
  write(u6,'(100a)') (('%'),i=1,96)
  write(u6,'(10x,a)') 'MOMLOC2:  EigenVectors of the Total Magnetic Interaction'
  write(u6,'(100a)') (('%'),i=1,96)

  write(u6,'(A,16A)') '-',('---',m=1,nsites),'|',('---------------------------|',i=1,4)
  write(fmtline,'(A,i2,A)') '(1x,',nsites,'A,39x,A,38x,A)'
  write(u6,fmtline) ('   ',m=1,nsites),'eigenvectors of the exchange matrix','|'
  do J=1,N,4
    jEnd = min(N,J+3)
    write(u6,'(A,16A)') '-',('---',m=1,nsites),'|',('---------------------------|',i=j,jEnd)

    write(u6,'(A,5A)') 'Exch.  |',('      exchange state       |',i=j,jEnd)

    write(u6,'(A,6A)') 'basis  |',('                           |',i=j,jEnd)
    write(u6,'(A,6(a,i5,a))') 'on site|',('         ',i,'             |',i=j,jEnd)
    write(fmtline,'(A,i2,A)') '(A,',nsites,'i3,a,6(a,f19.12,2x,a))'
    write(u6,fmtline) ' ',(i,i=1,nsites),'|',('  E =',wm(i)-wm(1),' |',i=j,jEnd)

    write(u6,'(A,16A)') '-',('---',m=1,nsites),'|',('----- Real ------- Imag ---|',i=j,jEnd)

    write(fmtline,'(A,i2,A)') '(A,',nsites,'i3,A,5(2F13.9,1x,a))'
    do iss=1,N
      write(u6,fmtline) '<',(ib(iss,m)+1,m=1,nsites),'|',(ZM(iss,i),'|',i=j,jEnd)
    end do  !iss
    write(u6,'(A,16A)') '-',('---',m=1,nsites),'|',('---------------------------|',i=j,jEnd)
    write(u6,*)
  end do ! j
end if !DBG

! generate the exchange basis:
do isite=1,nsites
  do L=1,3
    call zcopy_(N*N,[cZero],0,VL,1)
    do nb1=1,N
      do l_exch=1,nsites
        icoord(l_exch) = ib(nb1,l_exch)
      end do
      i1 = nind(isite,1)
      j1 = nind(isite,2)
      iss1 = ib(nb1,isite)+1

      do jss1=1,nexch(i1)
        icoord(isite) = jss1-1
        nb2 = norder(icoord,icod,nsites)
        VL(nb1,nb2) = VL(nb1,nb2)+r_rot(i1,j1,l,1)*dipso(i1,1,iss1,jss1)+r_rot(i1,j1,l,2)*dipso(i1,2,iss1,jss1)+ &
                      r_rot(i1,j1,l,3)*dipso(i1,3,iss1,jss1)
      end do  ! jss1
    end do  ! nb1
    ! rotate this matrix to exchange basis:
    call ZGEMM_('C','N',N,N,N,cOne,Z,N,VL,N,cZero,TMP,N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP,N,Z,N,cZero,VL,N)
    ! rotate this matrix to Zeeman basis:
    call ZGEMM_('C','N',N,N,N,cOne,ZM(1:N,1:N),N,VL(1:N,1:N),N,cZero,TMP(1:N,1:N),N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP(1:N,1:N),N,ZM(1:N,1:N),N,cZero,VL(1:N,1:N),N)
    do i=1,N
      MM(isite,L,i) = real(VL(i,i))
    end do
    ! spin moment
    call zcopy_(N*N,[cZero],0,VL,1)
    do nb1=1,N
      do l_exch=1,nsites
        icoord(l_exch) = ib(nb1,l_exch)
      end do
      i1 = nind(isite,1)
      j1 = nind(isite,2)
      iss1 = ib(nb1,isite)+1

      do jss1=1,nexch(i1)
        icoord(isite) = jss1-1
        nb2 = norder(icoord,icod,nsites)
        VL(nb1,nb2) = VL(nb1,nb2)+r_rot(i1,j1,l,1)*s_so(i1,1,iss1,jss1)+r_rot(i1,j1,l,2)*s_so(i1,2,iss1,jss1)+ &
                      r_rot(i1,j1,l,3)*s_so(i1,3,iss1,jss1)
      end do  ! jss1
    end do  ! nb1
    ! rotate this matrix to exchange basis and to Zeeman basis
    call ZGEMM_('C','N',N,N,N,cOne,Z(1:N,1:N),N,VL(1:N,1:N),N,cZero,TMP(1:N,1:N),N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP(1:N,1:N),N,Z(1:N,1:N),N,cZero,VL(1:N,1:N),N)
    call ZGEMM_('C','N',N,N,N,cOne,ZM(1:N,1:N),N,VL(1:N,1:N),N,cZero,TMP(1:N,1:N),N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP(1:N,1:N),N,ZM(1:N,1:N),N,cZero,VL(1:N,1:N),N)
    do i=1,N
      SM(isite,L,i) = real(VL(i,i))
    end do
  end do  ! L
end do  ! isite
! compute and print the calculated expectation values:
write(u6,*)
write(u6,'(A)') 'EXPECTATION VALUES'
write(u6,'(5A)') '--------|----|',('------------------------------|',i=1,4)
write(u6,'(A)') 'Exchange|Site|   MAGNETIC MOMENT (M=-L-2S)  |        SPIN MOMENT (S)       |      ORBITAL MOMENT (L)      |'// &
                '     TOTAL MOMENT (J=L+S)     |'
write(u6,'(5A)') ' state  | Nr.|',('     X         Y         Z    |',i=1,4)
write(u6,'(5A)') '--------|----|',('------------------------------|',i=1,4)
do i=1,N
 ! we proceed to compute expectation values for this nb1 exchange state
  do isite=1,nsites
    do L=1,3
      LM(isite,L,i) = -MM(isite,L,i)-g_e*SM(isite,L,i)
      JM(isite,L,i) = LM(isite,L,i)+SM(isite,L,i)
    end do
  end do

  do isite=1,nsites
    if (isite == int((nsites+1)/2)) then
      write(u6,'(i5,3x,A,1x,i2,1x,A,4(3(F9.5,1x),A))') i,'|',isite,'|',MM(isite,1,i),MM(isite,2,i),MM(isite,3,i),'|', &
                                                       SM(isite,1,i),SM(isite,2,i),SM(isite,3,i),'|',LM(isite,1,i),LM(isite,2,i), &
                                                       LM(isite,3,i),'|',JM(isite,1,i),JM(isite,2,i),JM(isite,3,i),'|'
    else
      write(u6,'(8x,A,1x,i2,1x,A,4(3(F9.5,1x),A))') '|',isite,'|',MM(isite,1,i),MM(isite,2,i),MM(isite,3,i),'|',SM(isite,1,i), &
                                                    SM(isite,2,i),SM(isite,3,i),'|',LM(isite,1,i),LM(isite,2,i),LM(isite,3,i), &
                                                    '|',JM(isite,1,i),JM(isite,2,i),JM(isite,3,i),'|'
    end if
  end do
  write(u6,'(5A)') '--------|----|',('------------------------------|',i1=1,4)
end do

199 continue

! deallocate temporary arrays:
call mma_deallocate(WM)
call mma_deallocate(MM)
call mma_deallocate(LM)
call mma_deallocate(SM)
call mma_deallocate(JM)
call mma_deallocate(ZM)
call mma_deallocate(VL)
call mma_deallocate(TMP)
call mma_deallocate(RWORK)
call mma_deallocate(HZEE)
call mma_deallocate(WORK)
call mma_deallocate(W_c)

return

end subroutine momloc2
