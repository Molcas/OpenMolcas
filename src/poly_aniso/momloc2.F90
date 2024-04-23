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

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne, gElectron
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, NL, nneq, neq(nneq), neqv, nsites, nexch(nneq)
real(kind=wp), intent(in) :: R_rot(NNEQ,neqv,3,3), W(N)
complex(kind=wp), intent(in) :: Z(N,N), dipexch(3,N,N), s_exch(3,N,N), dipso(nneq,3,NL,NL), s_so(nneq,3,NL,NL)
integer(kind=iwp) :: i, i1, ib(N,nsites), icod(nsites), icoord(nsites), isite, iss1, j, j1, jss1, k, L, nb1, nb2, nind(nsites,2), &
                     nmult, norder
real(kind=wp) :: E_thres, gtens(3), H, maxes(3,3), st(3), zJ
real(kind=wp), allocatable :: JM(:,:,:), LM(:,:,:), MM(:,:,:), RWORK(:), SM(:,:,:), WM(:)
complex(kind=wp), allocatable :: aux(:,:,:), HZEE(:), TMP(:,:), VL(:,:), W_c(:), WORK(:), ZM(:,:)
real(kind=wp), parameter :: g_e = -gElectron
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
integer(kind=iwp) :: iss, jEnd, m
character(len=60) :: fmtline
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: DBG = _DBG_

if (N > 1) return

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
call mma_allocate(HZEE,(nTri_Elem(N)),'ZEEM_HZEE')
call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')

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
  if (abs(W(i)-W(1)) < E_thres) nmult = nmult+1
end do
if (nmult < 2) nmult = 2 !minimum value needed to compute g-tensor
! find the main magnetic axes of this manifold:
if (N < nmult) then
  ! FIXME: this is a hack, atens should be made to work with size 1
  call mma_allocate(aux,3,nmult,nmult,label='aux')
  aux(:,:,:) = cZero
  aux(:,1:N,1:N) = dipexch(:,:,:)
  call atens(aux,nmult,gtens,maxes,1)
  call mma_deallocate(aux)
else
  call atens(dipexch,nmult,gtens,maxes,1)
end if
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'MOMLOC2:  g tensor of the ground manifold:'
write(u6,*)
write(u6,'((A,F12.6,A,3F12.7))') 'gX=',gtens(1),' axis X: ',(maxes(j,1),j=1,3)
write(u6,'((A,F12.6,A,3F12.7))') 'gY=',gtens(2),' axis Y: ',(maxes(j,2),j=1,3)
write(u6,'((A,F12.6,A,3F12.7))') 'gZ=',gtens(3),' axis Z: ',(maxes(j,3),j=1,3)
#endif
! construct the Zeeman matrix in the lowest N exchange states  and diagonalize it
st(:) = Zero
H = 1.0e-4_wp ! tesla
zJ = Zero ! absence of intermolecular interaction
call zeem_sa(N,H,maxes(1,3),maxes(2,3),maxes(3,3),W,dipexch,s_exch,ST,zJ,WM,ZM,DBG,RWORK,HZEE,WORK,W_c)
!cc  eigenvectors
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,'(a)') repeat('%',96)
write(u6,'(10x,a)') 'MOMLOC2:  EigenVectors of the Total Magnetic Interaction'
write(u6,'(a)') repeat('%',96)

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
#endif

! generate the exchange basis:
do isite=1,nsites
  do L=1,3
    VL(:,:) = Zero
    do nb1=1,N
      icoord(:) = ib(nb1,:)
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
    call ZGEMM_('C','N',N,N,N,cOne,ZM,N,VL,N,cZero,TMP,N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP,N,ZM,N,cZero,VL,N)
    do i=1,N
      MM(isite,L,i) = real(VL(i,i))
    end do
    ! spin moment
    VL(:,:) = Zero
    do nb1=1,N
      icoord(:) = ib(nb1,:)
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
    call ZGEMM_('C','N',N,N,N,cOne,Z,N,VL,N,cZero,TMP,N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP,N,Z,N,cZero,VL,N)
    call ZGEMM_('C','N',N,N,N,cOne,ZM,N,VL,N,cZero,TMP,N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP,N,ZM,N,cZero,VL,N)
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
  LM(:,:,i) = -MM(:,:,i)-g_e*SM(:,:,i)
  JM(:,:,i) = LM(:,:,i)+SM(:,:,i)

  do isite=1,nsites
    if (isite == (nsites+1)/2) then
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
