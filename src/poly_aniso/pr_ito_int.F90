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

subroutine pr_ito_int(npair,i_pair,lmax,nexch,nneq,neqv,itype,neq,nmax,soe,MM,SM,rot,Dipol,AnisoLines1,AnisoLines3,AnisoLines9, &
                      DM_exchange,JITO_exchange,HLIN1,HLIN3,HLIN9,HDIP,HDMO,HITO)
! this function prints the parameters of the exchange interaction in an accessible format

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: npair, i_pair(npair,2), lmax, nneq, nexch(nneq), neqv, neq(nneq), nmax
character, intent(in) :: itype(nneq)
real(kind=wp), intent(in) :: soe(nneq,nmax), rot(nneq,neqv,3,3)
complex(kind=wp), intent(in) :: MM(nneq,3,nmax,nmax), SM(nneq,3,nmax,nmax), HLIN1(npair,nmax,nmax,nmax,nmax), &
                                HLIN3(npair,nmax,nmax,nmax,nmax), HLIN9(npair,nmax,nmax,nmax,nmax), &
                                HDIP(npair,nmax,nmax,nmax,nmax), HDMO(npair,nmax,nmax,nmax,nmax), HITO(npair,nmax,nmax,nmax,nmax)
logical(kind=iwp), intent(in) :: Dipol, AnisoLines1, AnisoLines3, AnisoLines9, DM_exchange, JITO_exchange
! local variables
integer(kind=iwp) :: i, i1, i2, ibuf, j, k, k1, k2, l, lb1, lb2, lp, n1, n2, nind(lmax,2), q1, q2
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: is1, is2, js1, js2
#endif
real(kind=wp) :: g1(3), g2(3), J1C(3,3), J1Cr(3,3), mg1(3,3), mg2(3,3), RDI, RDM, RIT, RL1, RL3, RL9
complex(kind=wp) :: J_tmp(-1:1,-1:1)
complex(kind=wp), allocatable :: HTMP(:,:,:,:), JB(:,:,:,:), JN(:,:,:,:), JS(:,:,:,:), MM_tmp1(:,:,:), MM_tmp2(:,:,:)
!real(kind=wp), parameter :: cm_to_MHz = cm_s*1.0e-6_wp
real(kind=wp), external :: dznrm2_

#include "macros.fh"
unused_var(SM)
#ifndef _DEBUGPRINT_
unused_var(itype)
unused_var(soe)
unused_var(rot)
#endif

! some initializations:
l = 0
do i=1,nneq
  do j=1,neq(i)
    l = l+1
    nind(l,1) = i
    nind(l,2) = j
  end do
end do
nind(l+1:,:) = 0

ibuf = npair*nmax*nmax*nmax*nmax
if (ibuf == 0) then
  write(u6,'(A)') 'in UTMU:   ibuf=0 !!!'
  write(u6,*) 'npair= ',npair
  write(u6,*) 'nmax = ',nmax
  call xFlush(u6)
  call xquit(128)
end if
RL1 = dznrm2_(ibuf,HLIN1,1)
RL3 = dznrm2_(ibuf,HLIN3,1)
RL9 = dznrm2_(ibuf,HLIN9,1)
RDI = dznrm2_(ibuf,HDIP,1)
RDM = dznrm2_(ibuf,HDMO,1)
RIT = dznrm2_(ibuf,HITO,1)
!ccccccccccccccccccccccccccccccc
#ifdef _DEBUGPRINT_
write(u6,'(A,i6)') ' nmax =',nmax
write(u6,'(A,i6)') ' lmax =',lmax
write(u6,'(A,i6)') 'npair =',nmax
write(u6,'(A,i6)') ' nneq =',nneq
write(u6,'(20A)') ' itype =',(itype(i),i=1,nneq)
write(u6,'(A)') 'i_pair'
do i=1,npair
  write(u6,'(2I5)') i_pair(i,1),i_pair(i,2)
end do
write(u6,'(A)') 'equivalent sites'
do i=1,nneq
  write(u6,'(A,i2,A,i4)') 'neq(',i,') =',neq(i)
end do
write(u6,'(A)') 'exchange basis'
do i=1,nneq
  write(u6,'(A,i2,A,i4)') 'nexch(',i,') =',nexch(i)
end do
write(u6,'(A)') 'rotation matrix'
do i=1,nneq
  write(u6,'(A,i3)') 'site',i
  do j=1,neq(i)
    do l=1,3
      write(u6,'(3F12.6)') (rot(i,j,l,i1),i1=1,3)
    end do
  end do
end do
write(u6,'(A,i3)') 'site',i
write(u6,'(A)') 'magnetic moment, initial'
do i=1,nneq
  write(u6,'(A ,i3)') 'site',i
  call prMom('pr_ito_int:: magnetic moment, initial',mm(i,:,:,:),nmax)
end do
write(u6,'(A)') 'spin-orbit energies'
do i=1,nmax
  write(u6,'(1 5F14.5)') (soe(j,i),j=1,nneq)
end do

do lp=1,npair
  lb1 = i_pair(lp,1)
  lb2 = i_pair(lp,2)
  i1 = nind(lb1,1) ! indices of non-equivalent sites
  i2 = nind(lb2,1) ! indices of non-equivalent sites
  !j1 = nind(lb1,2) ! indices of equivalent sites
  !j2 = nind(lb2,2) ! indices of equivalent sites
  if (AnisoLines1 .and. (RL1 > Zero)) then
    write(u6,'(A,i5)') 'HLIN1,  interacting pair ',lp
    do is1=1,nexch(i1)
      do is2=1,nexch(i1)
        do js1=1,nexch(i2)
          write(u6,'(10(2F10.6,2x))') (HLIN1(lp,is1,is2,js1,js2),js2=1,nexch(i2))
        end do
      end do
    end do
  end if

  if (AnisoLines3 .and. (RL3 > Zero)) then
    write(u6,'(A,i5)') 'HLIN3,  interacting pair ',lp
    do is1=1,nexch(i1)
      do is2=1,nexch(i1)
        do js1=1,nexch(i2)
          write(u6,'(10(2F10.6,2x))') (HLIN3(lp,is1,is2,js1,js2),js2=1,nexch(i2))
        end do
      end do
    end do
  end if

  if (AnisoLines9 .and. (RL9 > Zero)) then
    write(u6,'(A,i5)') 'HLIN9,  interacting pair ',lp
    do is1=1,nexch(i1)
      do is2=1,nexch(i1)
        do js1=1,nexch(i2)
          write(u6,'(10(2F10.6,2x))') (HLIN9(lp,is1,is2,js1,js2),js2=1,nexch(i2))
        end do
      end do
    end do
  end if

  if (Dipol .and. (RDI > Zero)) then
    write(u6,'(A,i5)') 'HDIP,  interacting pair ',lp
    do is1=1,nexch(i1)
      do is2=1,nexch(i1)
        do js1=1,nexch(i2)
          write(u6,'(10(2F10.6,2x))') (HDIP(lp,is1,is2,js1,js2),js2=1,nexch(i2))
        end do
      end do
    end do
  end if

  if (DM_exchange .and. (RDM > Zero)) then
    write(u6,'(A,i5)') 'HDMO,  interacting pair ',lp
    do is1=1,nexch(i1)
      do is2=1,nexch(i1)
        do js1=1,nexch(i2)
          write(u6,'(10(2F10.6,2x))') (HDMO(lp,is1,is2,js1,js2),js2=1,nexch(i2))
        end do
      end do
    end do
  end if

  if (JITO_exchange .and. (RIT > Zero)) then
    write(u6,'(A,i5)') 'HITO,  interacting pair ',lp
    do is1=1,nexch(i1)
      do is2=1,nexch(i1)
        do js1=1,nexch(i2)
          write(u6,'(10(2F10.6,2x))') (HITO(lp,is1,is2,js1,js2),js2=1,nexch(i2))
        end do
      end do
    end do
  end if

end do !lp
!call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
!call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
!SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
!SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
!call prMom('SM(i1) bf Lines1',SM_tmp1,n1)
!call prMom('SM(i2) bf Lines1',SM_tmp2,n2)
!call mma_deallocate(SM_tmp1)
!call mma_deallocate(SM_tmp2)
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
write(u6,*)
write(u6,'(a)') repeat('%',100)
if ((.not. AnisoLines1) .and. (.not. AnisoLines3) .and. (.not. AnisoLines9) .and. (.not. Dipol) .and. (.not. DM_exchange) .and. &
    (.not. JITO_exchange)) then
  write(u6,'(20x,A)') 'ITO decomposition of exchange and/or dipolar couplings.'
  write(u6,'(A)') 'AnisoLines1 =  FALSE.'
  write(u6,'(A)') 'AnisoLines3 =  FALSE.'
  write(u6,'(A)') 'AnisoLines9 =  FALSE.'
  write(u6,'(A)') 'Dipol       =  FALSE.'
  write(u6,'(A)') 'JITO_exch   =  FALSE.'
  write(u6,'(A)') 'DM_exchange =  FALSE.'
  write(u6,'(A)') 'Nothing to Do.'
  return

else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. DM_exchange) .and. &
         (.not. JITO_exchange)) then
  write(u6,'(20x,A)') 'ITO decomposition of the Lines exchange interactions.'

else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. DM_exchange) .and. &
         (.not. JITO_exchange)) then
  write(u6,'(20x,A)') 'ITO decomposition of the dipolar interactions.'

else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. DM_exchange) .and. (.not. JITO_exchange)) then
  write(u6,'(20x,A)') 'ITO decomposition of exchange and/or dipolar couplings.'
else if (.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. DM_exchange) .and. &
         JITO_exchange) then
  write(u6,'(20x,A)') 'ITO decomposition of anisotropic exchange interaction. '
else if (.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. DM_exchange) .and. JITO_exchange) then
  write(u6,'(20x,A)') 'ITO decomposition of anisotropic exchange interaction and dipolar couplings.'

end if
write(u6,'(a)') repeat('%',100)
write(u6,*)

! decompose the exchange interaction in products of ITO
! first rotate the magnetic moments to the general coordinate system:
do lp=1,npair
  lb1 = i_pair(lp,1)
  lb2 = i_pair(lp,2)
  i1 = nind(lb1,1) ! indices of non-equivalent sites
  i2 = nind(lb2,1) ! indices of non-equivalent sites
  !j1 = nind(lb1,2) ! indices of equivalent sites
  !j2 = nind(lb2,2) ! indices of equivalent sites

  n1 = nexch(i1)
  n2 = nexch(i2)
  write(u6,'(A)') 'PART 1: Magnetic exchange is written in the coordinate systems of the LOCAL main magnetic axes of the '// &
                  'interacting sites.'
  write(u6,'(A)')
  write(u6,'(A)') repeat('-',100)
  write(u6,'(A,i2)') 'Interacting pair',lp

  call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
  call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
  MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
  MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

  ! JN= exch. parameters in Naoya's ITO operators
  ! JL= exch. parameters in Liviu's ITO operators
  ! JS= exch. parameters in Stevens ESO operators
  call mma_allocate(JN,[1,n1-1],[-(n1-1),n1-1],[1,n2-1],[-(n2-1),n2-1],'JN')
  call mma_allocate(JB,[1,n1-1],[-(n1-1),n1-1],[1,n2-1],[-(n2-1),n2-1],'JB')
  call mma_allocate(JS,[1,n1-1],[-(n1-1),n1-1],[1,n2-1],[-(n2-1),n2-1],'JS')

  !=====================================================================
  if (AnisoLines1 .and. (RL1 > Zero)) then

    call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
    HTMP(:,:,:,:) = HLIN1(lp,1:n1,1:n1,1:n2,1:n2)
    call newjkqpar(n1,n2,HTMP,JN,JB,JS)
    call mma_deallocate(HTMP)

    ! re-write the first rank tensor in cartesian representation:
    ! using Naoya's ITO parameters
    J1Cr(:,:) = Zero
    J_tmp(:,:) = JB(1,-1:1,1,-1:1)
    call tensor2cart(J_tmp,J1C)

    ! rotate cartesian J1C matrix by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM_tmp1,n1,g1,mg1,2)
    call atens(MM_tmp2,n2,g2,mg2,2)
    do j=1,3
      do l=1,3
        do k=1,3
          J1Cr(:,j) = J1Cr(:,j)+mg1(:,l)*mg2(j,k)*J1C(l,k)
        end do
      end do
    end do

    write(u6,'(A)')
    write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES-1'
    write(u6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(u6,'(A)') 'LOCAL AXES:::'
    write(u6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(u6,'(A)') ' H = -J * S1 * S2'
    write(u6,'(A)') '     (  xx   xy  xz  )  '
    write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(u6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(u6,'(3ES22.14)') (-J1Cr(i,j),j=1,3)
    end do

    ! print out the data:
    write(u6,'(A)')
    write(u6,'(10x,A)') 'Parameters of the ITOs: (Liviu ITO)'
    write(u6,'( 5x,A)') 'with absolute values larger than:  0.5e-14 '
    write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(u6,'(A)') ' rank | proj.| rank | proj.|           Lines  Exchange  Interaction        |'
    write(u6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'

    do k1=1,n1-1,2
      do q1=-k1,k1
        do k2=1,n2-1,2
          do q2=-k2,k2
            write(u6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(u6,'(A)') '------|------|------|------|-----------------------------------------------|'

    !-------------------------------------------------------------------
    ! verify the back transform for HAM!:
    !call mma_allocate(HAM,n1,n1,n2,n2,'HAM')
    !
    !S1a(:,:,:) = cZero
    !S2a(:,:,:) = cZero
    !
    !call mma_allocate(tmp1,n1,n1,label='tmp1')
    !call mma_allocate(tmp2,n1,n1,label='tmp2')
    !call ESO(n1,1,1,tmp1,tmp2,redME)
    !S1b(1,1:n1,1:n1) = tmp1(:,:)
    !S1b(2,1:n1,1:n1) = tmp2(:,:)
    !call ESO(n1,1,0,tmp1,tmp2,redME)
    !S1b(3,1:n1,1:n1) = tmp1(:,:)
    !W1b(1:n1,1:n1) = tmp2(:,:)
    !S2b(:,:,:) = S1b(:,:,:)
    !call mma_deallocate(tmp1)
    !call mma_deallocate(tmp2
    !
    !call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
    !call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
    !SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
    !SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
    !call prmom('SM(i1) bf recover HAM',SM_tmp,n1)
    !call prmom('SM(i2) bf recover HAM',SM_tmp,n2)
    !
    !J1C_trans(:,:) = -J1Cr(:,:)
    !write(u6,'(/)')
    !SM_tmp1(:,:,:) = S1b(:,1:n1,1:n1)
    !SM_tmp2(:,:,:) = S2b(:,1:n2,1:n2)
    !call Aniso_Lines_Exchange9(J1C_trans,n1,n2,SM_tmp1,SM_tmp2,HAM)
    !SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
    !SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
    !call Aniso_Lines_Exchange9(J1C_trans,n1,n2,SM_tmp1,SM_tmp2,HAM)
    !call mma_deallocate(SM_tmp1)
    !call mma_deallocate(SM_tmp2)
    !
    !write(u6,'(A,i5)') 'HLIN1: ORIG, REGEN, DIFF:'
    !do is1=1,n1
    !  do is2=1,n1
    !    do js1=1,n2
    !      do js2=1,n2
    !        write(u6,'(4I3,3(2ES22.13,3x))') is1,is2,js1,js2,HLIN1(lp,is1,is2,js1,js2),HAM(is1,is2,js1,js2), &
    !                                         HLIN1(lp,is1,is2,js1,js2)-HAM(is1,is2,js1,js2)
    !      end do
    !    end do
    !  end do
    !end do
    !
    !call mma_deallocate(HAM)
  end if

  !====================================================================
  if (AnisoLines3 .and. (RL3 > Zero)) then

    call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
    HTMP(:,:,:,:) = HLIN3(lp,1:n1,1:n1,1:n2,1:n2)
    call newjkqpar(n1,n2,HTMP,JN,JB,JS)
    call mma_deallocate(HTMP)

    J1Cr(:,:) = Zero
    ! using Liviu's ITO parameters
    J_tmp(:,:) = JB(1,-1:1,1,-1:1)
    call tensor2cart(J_tmp,J1C)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM_tmp1,n1,g1,mg1,2)
    call atens(MM_tmp2,n2,g2,mg2,2)
    do j=1,3
      do l=1,3
        do k=1,3
          J1Cr(:,j) = J1Cr(:,j)+mg1(:,l)*mg2(j,k)*J1C(l,k)
        end do
      end do
    end do

    write(u6,'(A)')
    write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES-3'
    write(u6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(u6,'(A)') 'LOCAL AXES:::'
    write(u6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(u6,'(A)') ' H = -J * S1 * S2'
    write(u6,'(A)') '     (  xx   xy  xz  )  '
    write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(u6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(u6,'(3ES22.14)') (-J1Cr(i,j),j=1,3)
    end do

    ! print out the data:
    write(u6,'(A)')
    write(u6,'(10x,A)') 'Parameters of the ITOs: (Liviu ITO)'
    write(u6,'( 5x,A)') 'with absolute values larger than:  1.0e-21 '
    write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(u6,'(A)') ' rank | proj.| rank | proj.|        Lines-3  Exchange  Interaction         |'
    write(u6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            !if (ABS(JLin3(lp,k1,q1,k2,q2)) > 1.0e-21_wp) then
            write(u6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
            !end if
          end do
        end do
      end do
    end do
    write(u6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================
  if (AnisoLines9 .and. (RL9 > Zero)) then

    call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
    HTMP(:,:,:,:) = HLIN9(lp,1:n1,1:n1,1:n2,1:n2)
    call newjkqpar(n1,n2,HTMP,JN,JB,JS)
    call mma_deallocate(HTMP)

    J1Cr(:,:) = Zero
    J_tmp(:,:) = JB(1,-1:1,1,-1:1)
    call tensor2cart(J_tmp,J1C)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM_tmp1,n1,g1,mg1,2)
    call atens(MM_tmp2,n2,g2,mg2,2)
    do j=1,3
      do l=1,3
        do k=1,3
          J1Cr(:,j) = J1Cr(:,j)+mg1(:,l)*mg2(j,k)*J1C(l,k)
        end do
      end do
    end do

    write(u6,'(A)')
    write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES-9'
    write(u6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(u6,'(A)') 'LOCAL AXES:::'
    write(u6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(u6,'(A)') ' H = -J * S1 * S2'
    write(u6,'(A)') '     (  xx   xy  xz  )  '
    write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(u6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(u6,'(3ES22.14)') (-J1Cr(i,j),j=1,3)
    end do

    ! print out the data:
    write(u6,'(A)')
    write(u6,'(10x,A)') 'Parameters of the ITOs:'
    write(u6,'( 5x,A)') 'with absolute values larger than:  1.0e-21 '
    write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(u6,'(A)') ' rank | proj.| rank | proj.|        Lines-9  Exchange  Interaction         |'
    write(u6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(u6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(u6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  if (Dipol .and. (RDI > Zero)) then

    call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
    HTMP(:,:,:,:) = HDIP(lp,1:n1,1:n1,1:n2,1:n2)
    call newjkqpar(n1,n2,HTMP,JN,JB,JS)
    call mma_deallocate(HTMP)

    J1Cr(:,:) = Zero
    J_tmp(:,:) = JB(1,-1:1,1,-1:1)
    call tensor2cart(J_tmp,J1C)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM_tmp1,n1,g1,mg1,2)
    call atens(MM_tmp2,n2,g2,mg2,2)
    do j=1,3
      do l=1,3
        do k=1,3
          J1Cr(:,j) = J1Cr(:,j)+mg1(:,l)*mg2(j,k)*J1C(l,k)
        end do
      end do
    end do

    write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: DIPOL'
    write(u6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(u6,'(A)') 'LOCAL AXES:::'
    write(u6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(u6,'(A)') ' H = -J * S1 * S2'
    write(u6,'(A)') '     (  xx   xy  xz  )  '
    write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(u6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(u6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
    end do
    write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(u6,'(A)') ' rank | proj.| rank | proj.|      Dipolar  Exchange  Interaction           |'
    write(u6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(u6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(u6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  if (DM_exchange .and. (RDM > Zero)) then

    call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
    HTMP(:,:,:,:) = HDMO(lp,1:n1,1:n1,1:n2,1:n2)
    call newjkqpar(n1,n2,HTMP,JN,JB,JS)
    call mma_deallocate(HTMP)

    J1Cr(:,:) = Zero
    J_tmp(:,:) = JB(1,-1:1,1,-1:1)
    call tensor2cart(J_tmp,J1C)

    call atens(MM_tmp1,n1,g1,mg1,2)
    call atens(MM_tmp2,n2,g2,mg2,2)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    do j=1,3
      do l=1,3
        do k=1,3
          J1Cr(:,j) = J1Cr(:,j)+mg1(:,l)*mg2(j,k)*J1C(l,k)
        end do
      end do
    end do

    write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: Dzyaloshinski-Morya'
    write(u6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(u6,'(A)') 'LOCAL AXES:::'
    write(u6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(u6,'(A)') ' H = -J * S1 * S2'
    write(u6,'(A)') '     (  xx   xy  xz  )  '
    write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(u6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(u6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
    end do
    write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(u6,'(A)') ' rank | proj.| rank | proj.|     Dzyaloshinsky - Morya Interaction         |'
    write(u6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(u6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(u6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  if (JITO_exchange .and. (RIT > Zero)) then

    call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
    HTMP(:,:,:,:) = HITO(lp,1:n1,1:n1,1:n2,1:n2)
    call newjkqpar(n1,n2,HTMP,JN,JB,JS)
    call mma_deallocate(HTMP)

    J_tmp(:,:) = JB(1,-1:1,1,-1:1)
    call tensor2cart(J_tmp,J1C)

    J1Cr(:,:) = Zero
    call atens(MM_tmp1,n1,g1,mg1,2)
    call atens(MM_tmp2,n2,g2,mg2,2)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    do j=1,3
      do l=1,3
        do k=1,3
          J1Cr(:,j) = J1Cr(:,j)+mg1(:,l)*mg2(j,k)*J1C(l,k)
        end do
      end do
    end do

    write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: Anisotropic ITO exchange'
    write(u6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
    write(u6,'(A)') 'LOCAL AXES:::'
    write(u6,'(A)') '     (  xx   xy  xz  )  '
    write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(u6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(u6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
    end do
    write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(u6,'(A)') ' rank | proj.| rank | proj.|     Anisotropic ITO Exchange Interaction      |'
    write(u6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(u6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(u6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !write(u6,'(A)') repeat('-',100)
  !write(u6,*)
  !write(u6,'(A)') 'PART 2: Magnetic exchange is written in the GENERAL coordinate systems of the computed system.'
  !write(u6,'(A)') 'The cartesian Z axis of the system is chosen the quantisation axis for the local pseudospins of the two '// &
  !                'interacting sites.'
  !write(u6,'(A)')
  !
  !iopt = 2
  !call unitmat(unity,3)
  !
  !if (Lines .or. AnisoLines) then
  !
  !  call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
  !  call mma_allocate(HTMP3,n1,n1,n2,n2,label='HTMP3')
  !  call mma_allocate(J_tmp2,[1,n1-1],[-(n1-1),n1-1],[1,n2-1],[-(n2-1),n2-1],label='J_tmp2')
  !  HTMP(:,:,:,:) = HLIN(lp,1:n1,1:n1,1:n2,1:n2)
  !  call transHam(n1,n2,unity,unity,MM_tmp1,MM_tmp2,itype(i1),itype(i2),HTMP,HTMP3),iopt)
  !  HLIN3(lp,1:n1,1:n1,1:n2,1:n2) = HTMP3(:,:,:,:)
  !  call JKQPar(n1,n2,HTMP3,J_tmp2)
  !  JLinG(lp,1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1) = J_tmp2(:,:,:,:)
  !  J_tmp(:,:) = JLinG(lp,1,-1:1,1,-1:1)
  !  call tensor2cart(1,1,J_tmp,J1C)
  !  JLinCG(lp,:,:) = J1C(:,:)
  !  call mma_deallocate(HTMP)
  !  call mma_deallocate(HTMP3)
  !  call mma_deallocate(J_tmp2)
  !  write(u6,'(A)')
  !  write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES'
  !  write(u6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !  write(u6,'(A)') 'GENERAL COORD:::'
  !  write(u6,'(A)') '     (  xx   xy  xz  )  '
  !  write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
  !  write(u6,'(A)') '     (  zx   zy  zz  )  '
  !  do i=1,3
  !    write(u6,'(3ES22.14)') (JLinCG(lp,i,j),j=1,3)
  !  end do
  !end if
  !
  !if (Dipol) then
  !  call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
  !  call mma_allocate(HTMP3,n1,n1,n2,n2,label='HTMP3')
  !  call mma_allocate(J_tmp2,[1,n1-1],[-(n1-1),n1-1],[1,n2-1],[-(n2-1),n2-1],label='J_tmp2')
  !  HDIP(:,:,:,:) = HLIN(lp,1:n1,1:n1,1:n2,1:n2)
  !  call transHam(n1,n2,unity,unity,MM_tmp1,MM_tmp2,itype(i1),itype(i2),HTMP,HTMP3),iopt)
  !  HDIP3(lp,1:n1,1:n1,1:n2,1:n2) = HTMP3(:,:,:,:)
  !  call JKQPar(n1,n2,HTMP3,J_tmp2)
  !  JDipG(lp,1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1) = J_tmp2(:,:,:,:)
  !  J_tmp(:,:) = JDipG(lp,1,-1:1,1,-1:1)
  !  call tensor2cart(1,1,J_tmp,JDipCG(lp,:,:))
  !  JDipCG(lp,:,:) = J1C(:,:)
  !  call mma_deallocate(HTMP)
  !  call mma_deallocate(HTMP3)
  !  call mma_deallocate(J_tmp2)
  !  write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: DIPOL'
  !  write(u6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !  write(u6,'(A)') 'GENERAL COORD:::'
  !  write(u6,'(A)') '     (  xx   xy  xz  )  '
  !  write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
  !  write(u6,'(A)') '     (  zx   zy  zz  )  '
  !  do i=1,3
  !    write(u6,'(3ES22.14)') (JDipCG(lp,i,j),j=1,3)
  !  end do
  !end if
  !
  !write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES+DIPOL'
  !write(u6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !write(u6,'(A)') 'GENERAL COORD:::'
  !write(u6,'(A)') '     (  xx   xy  xz  )  '
  !write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
  !write(u6,'(A)') '     (  zx   zy  zz  )  '
  !do i=1,3
  !  write(u6,'(3F24.14)') ((JLinCG(lp,i,j)+JDipCG(lp,i,j))*cm_to_MHz,j=1,3)
  !end do
  !
  ! print out the data:
  !write(u6,'(A)')
  !write(u6,'(10x,A)') 'Parameters of the ITOs:'
  !write(u6,'( 5x,A)') 'with absolute values larger than:  0.5e-14 '
  !if ((Lines .or. AnisoLines) .and. (.not. Dipol)) then
  !  write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-------------------------------------|'
  !  write(u6,'(A)') ' rank | proj.| rank | proj.|    Lines  Exchange  Interaction     |'
  !  write(u6,'(A)') '------|------|------|------|------ Real ----------- Imag --------|'
  !  do k1=1,N1-1,2
  !    do q1=-k1,k1
  !      do k2=1,N2-1,2
  !        do q2=-k2,k2
  !          if (abs(JLinG(lp,k1,q1,k2,q2)) > 0.5e-14_wp) &
  !            write(u6,'(4(i4,2x,A),2(1x,ES17.10),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JLinG(lp,k1,q1,k2,q2),'|'
  !        end do
  !      end do
  !    end do
  !  end do
  !  write(u6,'(A)') '------|------|------|------|-------------------------------------|'
  !else if ((.not. (Lines .or. AnisoLines)) .and. Dipol) then
  !  write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-------------------------------------|'
  !  write(u6,'(A)') ' rank | proj.| rank | proj.|  Dipolar  Exchange  Interaction     |'
  !  write(u6,'(A)') '------|------|------|------|------ Real ----------- Imag --------|'
  !  do k1=1,N1-1,2
  !    do q1=-k1,k1
  !      do k2=1,N2-1,2
  !        do q2=-k2,k2
  !          if (abs(JDipG(lp,k1,q1,k2,q2)) > 0.5e-14_wp) &
  !            write(u6,'(4(i4,2x,A),2(1x,ES17.10),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JDipG(lp,k1,q1,k2,q2),'|'
  !        end do
  !      end do
  !    end do
  !  end do
  !  write(u6,'(A)') '------|------|------|------|-------------------------------------|'
  !else if ((Lines .or. AnisoLines) .and. Dipol) then
  !  write(u6,'(A)') '--- SITE 1 --|--- SITE 2 --|-------------------------------------|-------------------------------------|'
  !  write(u6,'(A)') ' rank | proj.| rank | proj.|    Lines  Exchange  Interaction     |  Dipolar  Exchange  Interaction     |'
  !  write(u6,'(A)') '------|------|------|------|------ Real ----------- Imag --------|------ Real ----------- Imag --------|'
  !  do k1=1,N1-1,2
  !    do q1=-k1,k1
  !      do k2=1,N2-1,2
  !        do q2=-k2,k2
  !          if ((abs(JLinG(lp,k1,q1,k2,q2)) > 0.5e-14_wp) .or. (abs(JDipG(lp,k1,q1,k2,q2)) > 0.5e-14_wp)) &
  !            write(u6,'(4(i4,2x,A),2(1x,ES17.10),1x,A,2(1x,ES17.10),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JLinG(lp,k1,q1,k2,q2), &
  !                                                                           '|',JDipG(lp,k1,q1,k2,q2),'|'
  !        end do
  !      end do
  !    end do
  !  end do
  !  write(u6,'(A)') '------|------|------|------|-------------------------------------|-------------------------------------|'
  !end if
  !
  !write(u6,'(A)')
  !write(u6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction'
  !Write(u6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !Write(u6,'(A)') '     (  xx   xy  xz  )  '
  !Write(u6,'(A)') 'J =  (  yx   yy  yz  )  '
  !Write(u6,'(A)') '     (  zx   zy  zz  )  '
  !
  ! express the rank-1 tensors in the sum of
  ! Isotrop part    C * unit matrix
  ! Symmetric part
  ! Antisymmetric part
  ! print out the data:
  !
  !JDipCG(lp,:,:) = Zero
  !if(Lines .or. AnisoLines) then
  !  write(u6,'(A)') 'Lines matrix:'
  !  do i=1,3
  !    write(u6,'(3F12.6)') (JLinCG(lp,i,j),j=1,3)
  !  end do
  !
  !  ELin(lp,:,:) = Zero
  !  ALin(lp,:,:) = Zero
  !  SLin(lp,:,:) = Zero
  !  ELin(lp,1,1) = (JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/Three
  !  ELin(lp,2,2) = (JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/Three
  !  ELin(lp,3,3) = (JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/Three
  !  ALin(lp,:,:) = (JLinCG(lp,:,:)-JLinCG(lp,:,:))*Half
  !  tsum = Zero
  !  do l=1,3
  !    tsum = tsum+JLinCG(lp,l,l)
  !  end do
  !  test(:,:) = Zero
  !  do is1=1,3
  !    test(is1,is1) = tsum/OneHalf
  !  end do
  !
  !  SLin(lp,:,:) = (JLinCG(lp,:,:)+JLinCG(lp,:,:)-test(:,:))*Half
  !  write(u6,'(A)') 'Elin * unit'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (ELin(lp,is1,is2),is2=1,3)
  !  end do
  !  write(u6,'(A)') '-2/3 * test'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (-test(is1,is2)/OneHalf,is2=1,3)
  !  end do
  !
  !  write(u6,'(A)') 'ALin:'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (ALin(lp,is1,is2),is2=1,3)
  !  end do
  !  write(u6,'(A)') 'SLin:'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (SLin(lp,is1,is2),is2=1,3)
  !  end do
  !end if
  !
  !if (Dipol) then
  !  write(u6,'(A)') 'Dipolar exchange matrix:'
  !  do i=1,3
  !    write(u6,'(3F12.6)') (JDipCG(lp,i,j),j=1,3)
  !  end do
  !
  !  EDip(lp,:,:) = Zero
  !  ADip(lp,:,:) = Zero
  !  SDip(lp,:,:) = Zero
  !  EDip(lp,1,1) = (JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/Three
  !  EDip(lp,2,2) = (JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/Three
  !  EDip(lp,3,3) = (JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/Three
  !  ADip(lp,:,:) = (JDipCG(lp,:,:)-JDipCG(lp,:,:))*Half
  !  tsum = Zero
  !  do l=1,3
  !    tsum = tsum+JDipCG(lp,l,l)
  !  end do
  !  test(:,:) = Zero
  !  do is1=1,3
  !    test(is1,is1) = tsum/OneHalf
  !  end do
  !  SDip(lp,:,:) = (JDipCG(lp,:,:)+JDipCG(lp,:,:)-test(:,:))*Half
  !  write(u6,'(A)') 'EDip * unit'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (EDip(lp,is1,is2),is2=1,3)
  !  end do
  !  write(u6,'(A)') '-2/3 * test'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (-test(is1,is2)/OneHalf,is2=1,3)
  !  end do
  !  write(u6,'(A)') 'ADip:'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (ADip(lp,is1,is2),is2=1,3)
  !  end do
  !  write(u6,'(A)') 'SDip:'
  !  do is1=1,3
  !    write(u6,'(3F12.6)') (SDip(lp,is1,is2),is2=1,3)
  !  end do
  !end if

  call mma_deallocate(MM_tmp1)
  call mma_deallocate(MM_tmp2)

  call mma_deallocate(JN)
  call mma_deallocate(JB)
  call mma_deallocate(JS)
end do ! lp, interacting pairs

return

end subroutine pr_ito_int
