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
! this function prints the parameters of the exchange interaction in an
! accessible format

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: npair
integer, intent(in) :: nneq
integer, intent(in) :: neqv
integer, intent(in) :: lmax
integer, intent(in) :: nmax
integer, intent(in) :: nexch(nneq)
integer, intent(in) :: neq(nneq)
integer, intent(in) :: i_pair(npair,2)
real(kind=8), intent(in) :: rot(nneq,neqv,3,3)
real(kind=8), intent(in) :: soe(nneq,nmax)
complex(kind=8), intent(in) :: MM(nneq,3,nmax,nmax)
complex(kind=8), intent(in) :: SM(nneq,3,nmax,nmax)
complex(kind=8), intent(in) :: HLIN1(npair,nmax,nmax,nmax,nmax)
complex(kind=8), intent(in) :: HLIN3(npair,nmax,nmax,nmax,nmax)
complex(kind=8), intent(in) :: HLIN9(npair,nmax,nmax,nmax,nmax)
complex(kind=8), intent(in) :: HDIP(npair,nmax,nmax,nmax,nmax)
complex(kind=8), intent(in) :: HDMO(npair,nmax,nmax,nmax,nmax)
complex(kind=8), intent(in) :: HITO(npair,nmax,nmax,nmax,nmax)
logical, intent(in) :: Dipol
logical, intent(in) :: AnisoLines1
logical, intent(in) :: AnisoLines3
logical, intent(in) :: AnisoLines9
logical, intent(in) :: DM_exchange
logical, intent(in) :: JITO_exchange
character(len=1), intent(in) :: itype(nneq)
! local variables
!integer :: iopt
integer :: i, j, l, k, lp, i1, i2, lb1, lb2, ibuf, is1, is2, js1, js2, k1, k2, q1, q2, n1, n2, nsize
integer :: nind(lmax,2), l1(2), l2(2), l3(2), l4(2)
real(kind=8) :: J1C(3,3), J1Cr(3,3) !, J1C_trans(3,3)
complex(kind=8), allocatable :: JN(:,:,:,:)
complex(kind=8), allocatable :: JB(:,:,:,:)
complex(kind=8), allocatable :: JS(:,:,:,:)
real(kind=8) :: dznrm2_, RL1, RL3, RL9, RDI, RDM, RIT
real(kind=8) :: g1(3), g2(3), mg1(3,3), mg2(3,3)
external :: dznrm2_
!real(kind=8) :: cm_to_MHz
logical DBG

!cm_to_MHz = 29979.2458_wp
DBG = .false.
! some initializations:
nind(:,:) = 0
l = 0
do i=1,nneq
  do j=1,neq(i)
    l = l+1
    nind(l,1) = i
    nind(l,2) = j
  end do
end do

ibuf = npair*nmax*nmax*nmax*nmax
if (ibuf == 0) then
  write(6,'(A)') 'in UTMU:   ibuf=0 !!!'
  write(6,*) 'npair= ',npair
  write(6,*) 'nmax = ',nmax
  call xFlush(6)
  call xquit(128)
end if
RL1 = dznrm2_(ibuf,HLIN1,1)
RL3 = dznrm2_(ibuf,HLIN3,1)
RL9 = dznrm2_(ibuf,HLIN9,1)
RDI = dznrm2_(ibuf,HDIP,1)
RDM = dznrm2_(ibuf,HDMO,1)
RIT = dznrm2_(ibuf,HITO,1)
!ccccccccccccccccccccccccccccccc
if (DBG) then
  write(6,'(A,i6)') ' nmax =',nmax
  write(6,'(A,i6)') ' lmax =',lmax
  write(6,'(A,i6)') 'npair =',nmax
  write(6,'(A,i6)') ' nneq =',nneq
  write(6,'(20A)') ' itype =',(itype(i),i=1,nneq)
  write(6,'(A)') 'i_pair'
  do i=1,npair
    write(6,'(2I5)') i_pair(i,1),i_pair(i,2)
  end do
  write(6,'(A)') 'equivalent sites'
  do i=1,nneq
    write(6,'(A,i2,A,i4)') 'neq(',i,') =',neq(i)
  end do
  write(6,'(A)') 'exchange basis'
  do i=1,nneq
    write(6,'(A,i2,A,i4)') 'nexch(',i,') =',nexch(i)
  end do
  write(6,'(A)') 'rotation matrix'
  do i=1,nneq
    write(6,'(A,i3)') 'site',i
    do j=1,neq(i)
      do l=1,3
        write(6,'(3F12.6)') (rot(i,j,l,i1),i1=1,3)
      end do
    end do
  end do
  write(6,'(A,i3)') 'site',i
  write(6,'(A)') 'magnetic moment, initial'
  do i=1,nneq
    write(6,'(A ,i3)') 'site',i
    call prMom('pr_ito_int:: magnetic moment, initial',mm(i,:,:,:),nmax)
  end do
  write(6,'(A)') 'spin-orbit energies'
  do i=1,nmax
    write(6,'(1 5F14.5)') (soe(j,i),j=1,nneq)
  end do

  do lp=1,npair
    lb1 = i_pair(lp,1)
    lb2 = i_pair(lp,2)
    i1 = nind(lb1,1) ! indices of non-equivalent sites
    i2 = nind(lb2,1) ! indices of non-equivalent sites
    !j1 = nind(lb1,2) ! indices of equivalent sites
    !j2 = nind(lb2,2) ! indices of equivalent sites
    if (AnisoLines1 .and. (RL1 > 0.0_wp)) then
      write(6,'(A,i5)') 'HLIN1,  interacting pair ',lp
      do is1=1,nexch(i1)
        do is2=1,nexch(i1)
          do js1=1,nexch(i2)
            write(6,'(10(2F10.6,2x))') (HLIN1(lp,is1,is2,js1,js2),js2=1,nexch(i2))
          end do
        end do
      end do
    end if

    if (AnisoLines3 .and. (RL3 > 0.0_wp)) then
      write(6,'(A,i5)') 'HLIN3,  interacting pair ',lp
      do is1=1,nexch(i1)
        do is2=1,nexch(i1)
          do js1=1,nexch(i2)
            write(6,'(10(2F10.6,2x))') (HLIN3(lp,is1,is2,js1,js2),js2=1,nexch(i2))
          end do
        end do
      end do
    end if

    if (AnisoLines9 .and. (RL9 > 0.0_wp)) then
      write(6,'(A,i5)') 'HLIN9,  interacting pair ',lp
      do is1=1,nexch(i1)
        do is2=1,nexch(i1)
          do js1=1,nexch(i2)
            write(6,'(10(2F10.6,2x))') (HLIN9(lp,is1,is2,js1,js2),js2=1,nexch(i2))
          end do
        end do
      end do
    end if

    if (Dipol .and. (RDI > 0.0_wp)) then
      write(6,'(A,i5)') 'HDIP,  interacting pair ',lp
      do is1=1,nexch(i1)
        do is2=1,nexch(i1)
          do js1=1,nexch(i2)
            write(6,'(10(2F10.6,2x))') (HDIP(lp,is1,is2,js1,js2),js2=1,nexch(i2))
          end do
        end do
      end do
    end if

    if (DM_exchange .and. (RDM > 0.0_wp)) then
      write(6,'(A,i5)') 'HDMO,  interacting pair ',lp
      do is1=1,nexch(i1)
        do is2=1,nexch(i1)
          do js1=1,nexch(i2)
            write(6,'(10(2F10.6,2x))') (HDMO(lp,is1,is2,js1,js2),js2=1,nexch(i2))
          end do
        end do
      end do
    end if

    if (JITO_exchange .and. (RIT > 0.0_wp)) then
      write(6,'(A,i5)') 'HITO,  interacting pair ',lp
      do is1=1,nexch(i1)
        do is2=1,nexch(i1)
          do js1=1,nexch(i2)
            write(6,'(10(2F10.6,2x))') (HITO(lp,is1,is2,js1,js2),js2=1,nexch(i2))
          end do
        end do
      end do
    end if

  end do !lp
  call prMom('SM(i1) bf Lines1',SM(i1,1:3,1:n1,1:n1),n1)
  call prMom('SM(i2) bf Lines1',SM(i2,1:3,1:n2,1:n2),n2)
end if !DBG

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
write(6,*)
write(6,'(100a)') (('%'),j=1,100)
if ((.not. AnisoLines1) .and. (.not. AnisoLines3) .and. (.not. AnisoLines9) .and. (.not. Dipol) .and. (.not. DM_exchange) .and. &
    (.not. JITO_exchange)) then
  write(6,'(20x,A)') 'ITO decomposition of exchange and/or dipolar couplings.'
  write(6,'(A)') 'AnisoLines1 =  FALSE.'
  write(6,'(A)') 'AnisoLines3 =  FALSE.'
  write(6,'(A)') 'AnisoLines9 =  FALSE.'
  write(6,'(A)') 'Dipol       =  FALSE.'
  write(6,'(A)') 'JITO_exch   =  FALSE.'
  write(6,'(A)') 'DM_exchange =  FALSE.'
  write(6,'(A)') 'Nothing to Do.'
  return

else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. DM_exchange) .and. &
         (.not. JITO_exchange)) then
  write(6,'(20x,A)') 'ITO decomposition of the Lines exchange interactions.'

else if ((.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9)) .and. Dipol .and. (.not. DM_exchange) .and. &
         (.not. JITO_exchange)) then
  write(6,'(20x,A)') 'ITO decomposition of the dipolar interactions.'

else if ((AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. DM_exchange) .and. (.not. JITO_exchange)) then
  write(6,'(20x,A)') 'ITO decomposition of exchange and/or dipolar couplings.'
else if (.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. (.not. Dipol) .and. (.not. DM_exchange) .and. &
         JITO_exchange) then
  write(6,'(20x,A)') 'ITO decomposition of anisotropic exchange interaction. '
else if (.not. (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) .and. Dipol .and. (.not. DM_exchange) .and. JITO_exchange) then
  write(6,'(20x,A)') 'ITO decomposition of anisotropic exchange interaction and dipolar couplings.'

end if
write(6,'(100a)') (('%'),j=1,100)
write(6,*)

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
  write(6,'(A)') 'PART 1: Magnetic exchange is written in the coordinate systems of the LOCAL main magnetic axes of the '// &
                 'interacting sites.'
  !iopt = 1
  write(6,'(A)')
  write(6,'(100A)') ('-',i=1,100)
  write(6,'(A,i2)') 'Interacting pair',lp

  ! JN= exch. parameters in Naoya's ITO operators
  ! JL= exch. parameters in Liviu's ITO operators
  ! JS= exch. parameters in Stevens ESO operators
  l1(1) = 1
  l1(2) = n1-1
  l2(1) = -(n1-1)
  l2(2) = n1-1
  l3(1) = 1
  l3(2) = n2-1
  l4(1) = -(n2-1)
  l4(2) = n2-1
  nsize = (n1-1)*(n2-1)*(2*(n1-1)+1)*(2*(n2-1)+1)
  call mma_allocate(JN,l1,l2,l3,l4,'JN')
  call mma_allocate(JB,l1,l2,l3,l4,'JB')
  call mma_allocate(JS,l1,l2,l3,l4,'JS')

  !=====================================================================
  if (AnisoLines1 .and. (RL1 > 0.0_wp)) then

    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)

    call newjkqpar(n1,n2,HLIN1(lp,1:n1,1:n1,1:n2,1:n2),JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1), &
                   JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))

    ! re-write the first rank tensor in cartesian representation:
    ! using Naoya's ITO parameters
    J1C = 0.0_wp
    J1Cr = 0.0_wp
    call tensor2cart(JB(1,-1:1,1,-1:1),J1C)

    ! rotate cartesian J1C matrix by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM(i1,1:3,1:n1,1:n1),n1,g1,mg1,2)
    call atens(MM(i2,1:3,1:n2,1:n2),n2,g2,mg2,2)
    do i=1,3
      do j=1,3
        do l=1,3
          do k=1,3
            J1Cr(i,j) = J1Cr(i,j)+mg1(i,l)*mg2(j,k)*J1C(l,k)
          end do
        end do
      end do
    end do

    write(6,'(A)')
    write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES-1'
    write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(6,'(A)') 'LOCAL AXES:::'
    write(6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(6,'(A)') ' H = -J * S1 * S2'
    write(6,'(A)') '     (  xx   xy  xz  )  '
    write(6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(6,'(3ES22.14)') (-J1Cr(i,j),j=1,3)
    end do

    ! print out the data:
    write(6,'(A)')
    write(6,'(10x,A)') 'Parameters of the ITOs: (Liviu ITO)'
    write(6,'( 5x,A)') 'with absolute values larger than:  0.5d-14 '
    write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(6,'(A)') ' rank | proj.| rank | proj.|           Lines  Exchange  Interaction        |'
    write(6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'

    do k1=1,n1-1,2
      do q1=-k1,k1
        do k2=1,n2-1,2
          do q2=-k2,k2
            write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(6,'(A)') '------|------|------|------|-----------------------------------------------|'

    !-------------------------------------------------------------------
    ! verify the back transform for HAM!:
    !call mma_allocate(HAM,n1,n1,n2,n2,'HAM')
    !call zcopy_(n1*n1*n2*n2,[(0.0_wp,0.0_wp)],0,HAM,1)
    !
    !S1a = (0.0_wp,0.0_wp)
    !S2a = (0.0_wp,0.0_wp)
    !S1b = (0.0_wp,0.0_wp)
    !S2b = (0.0_wp,0.0_wp)
    !
    !call ESO(n1,1,1,S1b(1,1:n1,1:n1),S1b(2,1:n1,1:n1),redME)
    !call ESO(n1,1,0,S1b(3,1:n1,1:n1),W1b(1:n1,1:n1),redME)
    !call zcopy_(3*n1*n1,S1b,1,S2b,1)
    !
    !call prmom('SM(i1) bf recover HAM',SM(i1,1:3,1:n1,1:n1),n1)
    !call prmom('SM(i2) bf recover HAM',SM(i2,1:3,1:n2,1:n2),n2)
    !
    !J1C_trans = 0.0_wp
    !do i=1,3
    !  do j=1,3
    !    J1C_trans(i,j) = -J1Cr(i,j)
    !  end do
    !end do
    !write(6,'(/)')
    !call Aniso_Lines_Exchange9(J1C_trans,n1,n2,S1b(1:3,1:n1,1:n1),S2b(1:3,1:n2,1:n2),HAM)
    !call Aniso_Lines_Exchange9(J1C_trans,n1,n2,SM(i1,1:3,1:n1,1:n1),SM(i2,1:3,1:n2,1:n2),HAM)
    !
    !write(6,'(A,i5)') 'HLIN1: ORIG, REGEN, DIFF:'
    !do is1=1,n1
    !  do is2=1,n1
    !    do js1=1,n2
    !      do js2=1,n2
    !        write(6,'(4I3,3(2ES22.13,3x))') is1,is2,js1,js2,HLIN1(lp,is1,is2,js1,js2),HAM(is1,is2,js1,js2), &
    !                                        HLIN1(lp,is1,is2,js1,js2)-HAM(is1,is2,js1,js2)
    !      end do
    !    end do
    !  end do
    !end do
    !
    !call mma_deallocate(HAM)
  end if

  !====================================================================
  if (AnisoLines3 .and. (RL3 > 0.0_wp)) then

    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)

    call newjkqpar(n1,n2,HLIN3(lp,1:n1,1:n1,1:n2,1:n2),JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1), &
                   JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))

    J1C = 0.0_wp
    J1Cr = 0.0_wp
    ! using Liviu's ITO parameters
    call tensor2cart(JB(1,-1:1,1,-1:1),J1C)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM(i1,1:3,1:n1,1:n1),n1,g1,mg1,2)
    call atens(MM(i2,1:3,1:n2,1:n2),n2,g2,mg2,2)
    do i=1,3
      do j=1,3
        do l=1,3
          do k=1,3
            J1Cr(i,j) = J1Cr(i,j)+mg1(i,l)*mg2(j,k)*J1C(l,k)
          end do
        end do
      end do
    end do

    write(6,'(A)')
    write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES-3'
    write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(6,'(A)') 'LOCAL AXES:::'
    write(6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(6,'(A)') ' H = -J * S1 * S2'
    write(6,'(A)') '     (  xx   xy  xz  )  '
    write(6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(6,'(3ES22.14)') (-J1Cr(i,j),j=1,3)
    end do

    ! print out the data:
    write(6,'(A)')
    write(6,'(10x,A)') 'Parameters of the ITOs: (Liviu ITO)'
    write(6,'( 5x,A)') 'with absolute values larger than:  0.1d-20 '
    write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(6,'(A)') ' rank | proj.| rank | proj.|        Lines-3  Exchange  Interaction         |'
    write(6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            !if (ABS(JLin3(lp,k1,q1,k2,q2)) > 1.0e-21_wp) then
            write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
            !end if
          end do
        end do
      end do
    end do
    write(6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================
  if (AnisoLines9 .and. (RL9 > 0.0_wp)) then

    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)

    call newjkqpar(n1,n2,HLIN9(lp,1:n1,1:n1,1:n2,1:n2),JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1), &
                   JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))

    J1C = 0.0_wp
    J1Cr = 0.0_wp
    call tensor2cart(JB(1,-1:1,1,-1:1),J1C)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM(i1,1:3,1:n1,1:n1),n1,g1,mg1,2)
    call atens(MM(i2,1:3,1:n2,1:n2),n2,g2,mg2,2)
    do i=1,3
      do j=1,3
        do l=1,3
          do k=1,3
            J1Cr(i,j) = J1Cr(i,j)+mg1(i,l)*mg2(j,k)*J1C(l,k)
          end do
        end do
      end do
    end do

    write(6,'(A)')
    write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES-9'
    write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(6,'(A)') 'LOCAL AXES:::'
    write(6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(6,'(A)') ' H = -J * S1 * S2'
    write(6,'(A)') '     (  xx   xy  xz  )  '
    write(6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(6,'(3ES22.14)') (-J1Cr(i,j),j=1,3)
    end do

    ! print out the data:
    write(6,'(A)')
    write(6,'(10x,A)') 'Parameters of the ITOs:'
    write(6,'( 5x,A)') 'with absolute values larger than:  0.1e-20 '
    write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(6,'(A)') ' rank | proj.| rank | proj.|        Lines-9  Exchange  Interaction         |'
    write(6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  if (Dipol .and. (RDI > 0.0_wp)) then

    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)

    call newjkqpar(n1,n2,HDIP(lp,1:n1,1:n1,1:n2,1:n2),JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1), &
                   JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))

    J1C = 0.0_wp
    J1Cr = 0.0_wp
    call tensor2cart(JB(1,-1:1,1,-1:1),J1C(1:3,1:3))

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    call atens(MM(i1,1:3,1:n1,1:n1),n1,g1,mg1,2)
    call atens(MM(i2,1:3,1:n2,1:n2),n2,g2,mg2,2)
    do i=1,3
      do j=1,3
        do l=1,3
          do k=1,3
            J1Cr(i,j) = J1Cr(i,j)+mg1(i,l)*mg2(j,k)*J1C(l,k)
          end do
        end do
      end do
    end do

    write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: DIPOL'
    write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(6,'(A)') 'LOCAL AXES:::'
    write(6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(6,'(A)') ' H = -J * S1 * S2'
    write(6,'(A)') '     (  xx   xy  xz  )  '
    write(6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
    end do
    write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(6,'(A)') ' rank | proj.| rank | proj.|      Dipolar  Exchange  Interaction           |'
    write(6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  if (DM_exchange .and. (RDM > 0.0_wp)) then

    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)

    call newjkqpar(n1,n2,HDMO(lp,1:n1,1:n1,1:n2,1:n2),JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1), &
                   JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))

    J1C = 0.0_wp
    J1Cr = 0.0_wp
    call tensor2cart(JB(1,-1:1,1,-1:1),J1C)

    call atens(MM(i1,1:3,1:n1,1:n1),n1,g1,mg1,2)
    call atens(MM(i2,1:3,1:n2,1:n2),n2,g2,mg2,2)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    do i=1,3
      do j=1,3
        do l=1,3
          do k=1,3
            J1Cr(i,j) = J1Cr(i,j)+mg1(i,l)*mg2(j,k)*J1C(l,k)
          end do
        end do
      end do
    end do

    write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: Dzyaloshinski-Morya'
    write(6,'(A)') 'Anisotropic exchange interaction: -J matrix:'
    write(6,'(A)') 'LOCAL AXES:::'
    write(6,'(A)') 'To be used directly in exchange Hamiltonian of the kind:'
    write(6,'(A)') ' H = -J * S1 * S2'
    write(6,'(A)') '     (  xx   xy  xz  )  '
    write(6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
    end do
    write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(6,'(A)') ' rank | proj.| rank | proj.|     Dzyaloshinsky - Morya Interaction         |'
    write(6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  if (JITO_exchange .and. (RIT > 0.0_wp)) then

    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)
    call zcopy_(nsize,[(0.0_wp,0.0_wp)],0,JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),1)

    call newjkqpar(n1,n2,HITO(lp,1:n1,1:n1,1:n2,1:n2),JN(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1), &
                   JB(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1),JS(1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))

    J1C(1:3,1:3) = 0.0_wp
    call tensor2cart(JB(1,-1:1,1,-1:1),J1C)

    call atens(MM(i1,1:3,1:n1,1:n1),n1,g1,mg1,2)
    call atens(MM(i2,1:3,1:n2,1:n2),n2,g2,mg2,2)

    ! rotate cartesian JLinC1 by mg1 and mg2, in order to represent
    ! the interaction matrix in the original coordinate system:
    J1Cr(1:3,1:3) = 0.0_wp
    do i=1,3
      do j=1,3
        do l=1,3
          do k=1,3
            J1Cr(i,j) = J1Cr(i,j)+mg1(i,l)*mg2(j,k)*J1C(l,k)
          end do
        end do
      end do
    end do

    write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: Anisotropic ITO exchange'
    write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
    write(6,'(A)') 'LOCAL AXES:::'
    write(6,'(A)') '     (  xx   xy  xz  )  '
    write(6,'(A)') 'J =  (  yx   yy  yz  )  '
    write(6,'(A)') '     (  zx   zy  zz  )  '
    do i=1,3
      write(6,'(3ES24.14)') (-J1Cr(i,j),j=1,3)
    end do
    write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-----------------------------------------------|'
    write(6,'(A)') ' rank | proj.| rank | proj.|     Anisotropic ITO Exchange Interaction      |'
    write(6,'(A)') '------|------|------|------|---------- Real ----------------- Imag --------|'
    do k1=1,N1-1,2
      do q1=-k1,k1
        do k2=1,N2-1,2
          do q2=-k2,k2
            write(6,'(4(i4,2x,A),2(1x,ES22.14),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JN(k1,q1,k2,q2),'|'
          end do
        end do
      end do
    end do
    write(6,'(A)') '------|------|------|------|-----------------------------------------------|'
  end if

  !=====================================================================

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  !write(6,'(100A)') ('-',i=1,100)
  !write(6,*)
  !write(6,'(A)') 'PART 2: Magnetic exchange is written in the GENERAL coordinate systems of the computed system.'
  !write(6,'(A)') 'The cartesian Z axis of the system is chosen the quantisation axis for the local pseudospins of the two '// &
  !               'interacting sites.'
  !write(6,'(A)')
  !
  !iopt = 2
  !unity(1:3,1:3) = 0.0_wp
  !unity(1,1) = 1.0_wp
  !unity(2,2) = 1.0_wp
  !unity(3,3) = 1.0_wp
  !
  !if (Lines .or. AnisoLines) then
  !
  !  call transHam(n1,n2,unity(1:3,1:3),unity(1:3,1:3),MM(i1,1:3,1:n1,1:n1),MM(i2,1:3,1:n2,1:n2),itype(i1),itype(i2), &
  !                HLIN(lp,1:n1,1:n1,1:n2,1:n2),HLIN3(lp,1:n1,1:n1,1:n2,1:n2),iopt)
  !  JLinG(lp,1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1) = (0.0_wp,0.0_wp)
  !  call JKQPar(n1,n2,HLIN3(lp,1:n1,1:n1,1:n2,1:n2),JLinG(lp,1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))
  !  JLinCG(lp,1:3,1:3) = 0.0_wp
  !  call tensor2cart(1,1,JLinG(lp,1,-1:1,1,-1:1),JLinCG(lp,1:3,1:3))
  !  write(6,'(A)')
  !  write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES'
  !  write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !  write(6,'(A)') 'GENERAL COORD:::'
  !  write(6,'(A)') '     (  xx   xy  xz  )  '
  !  write(6,'(A)') 'J =  (  yx   yy  yz  )  '
  !  write(6,'(A)') '     (  zx   zy  zz  )  '
  !  do i=1,3
  !    write(6,'(3ES22.14)') (JLinCG(lp,i,j),j=1,3)
  !  end do
  !end if
  !
  !if (Dipol) then
  !  call transHam(n1,n2,unity(1:3,1:3),unity(1:3,1:3),MM(i1,1:3,1:n1,1:n1),MM(i2,1:3,1:n2,1:n2),itype(i1),itype(i2), &
  !                HDIP(lp,1:n1,1:n1,1:n2,1:n2),HDIP3(lp,1:n1,1:n1,1:n2,1:n2),iopt)
  !  JDipG(lp,1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1) = (0.0_wp,0.0_wp)
  !  call JKQPar(n1,n2,HDIP3(lp,1:n1,1:n1,1:n2,1:n2),JDipG(lp,1:n1-1,-(n1-1):n1-1,1:n2-1,-(n2-1):n2-1))
  !  JDipCG(lp,1:3,1:3) = 0.0_wp
  !  call tensor2cart(1,1,JDipG(lp,1,-1:1,1,-1:1),JDipCG(lp,1:3,1:3))
  !  write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: DIPOL'
  !  write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !  write(6,'(A)') 'GENERAL COORD:::'
  !  write(6,'(A)') '     (  xx   xy  xz  )  '
  !  write(6,'(A)') 'J =  (  yx   yy  yz  )  '
  !  write(6,'(A)') '     (  zx   zy  zz  )  '
  !  do i=1,3
  !    write(6,'(3ES22.14)') (JDipCG(lp,i,j),j=1,3)
  !  end do
  !end if
  !
  !write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction: LINES+DIPOL'
  !write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !write(6,'(A)') 'GENERAL COORD:::'
  !write(6,'(A)') '     (  xx   xy  xz  )  '
  !write(6,'(A)') 'J =  (  yx   yy  yz  )  '
  !write(6,'(A)') '     (  zx   zy  zz  )  '
  !do i=1,3
  !  write(6,'(3F24.14)') ((JLinCG(lp,i,j)+JDipCG(lp,i,j))*cm_to_MHz,j=1,3)
  !end do
  !
  ! print out the data:
  !write(6,'(A)')
  !write(6,'(10x,A)') 'Parameters of the ITOs:'
  !write(6,'( 5x,A)') 'with absolute values larger than:  0.5d-14 '
  !if ((Lines .or. AnisoLines) .and. (.not. Dipol)) then
  !  write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-------------------------------------|'
  !  write(6,'(A)') ' rank | proj.| rank | proj.|    Lines  Exchange  Interaction     |'
  !  write(6,'(A)') '------|------|------|------|------ Real ----------- Imag --------|'
  !  do k1=1,N1-1,2
  !    do q1=-k1,k1
  !      do k2=1,N2-1,2
  !        do q2=-k2,k2
  !          if (abs(JLinG(lp,k1,q1,k2,q2)) > 0.5d-14) write(6,'(4(i4,2x,A),2(1x,ES17.10),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|', &
  !                                                                                               JLinG(lp,k1,q1,k2,q2),'|'
  !        end do
  !      end do
  !    end do
  !  end do
  !  write(6,'(A)') '------|------|------|------|-------------------------------------|'
  !else if ((.not. (Lines .or. AnisoLines)) .and. Dipol) then
  !  write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-------------------------------------|'
  !  write(6,'(A)') ' rank | proj.| rank | proj.|  Dipolar  Exchange  Interaction     |'
  !  write(6,'(A)') '------|------|------|------|------ Real ----------- Imag --------|'
  !  do k1=1,N1-1,2
  !    do q1=-k1,k1
  !      do k2=1,N2-1,2
  !        do q2=-k2,k2
  !          if (abs(JDipG(lp,k1,q1,k2,q2)) > 0.5d-14) write(6,'(4(i4,2x,A),2(1x,ES17.10),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|', &
  !                                                                                               JDipG(lp,k1,q1,k2,q2),'|'
  !        end do
  !      end do
  !    end do
  !  end do
  !  write(6,'(A)') '------|------|------|------|-------------------------------------|'
  !else if ((Lines .or. AnisoLines) .and. Dipol) then
  !  write(6,'(A)') '--- SITE 1 --|--- SITE 2 --|-------------------------------------|-------------------------------------|'
  !  write(6,'(A)') ' rank | proj.| rank | proj.|    Lines  Exchange  Interaction     |  Dipolar  Exchange  Interaction     |'
  !  write(6,'(A)') '------|------|------|------|------ Real ----------- Imag --------|------ Real ----------- Imag --------|'
  !  do k1=1,N1-1,2
  !    do q1=-k1,k1
  !      do k2=1,N2-1,2
  !        do q2=-k2,k2
  !          if ((abs(JLinG(lp,k1,q1,k2,q2)) > 0.5d-14) .or. (abs(JDipG(lp,k1,q1,k2,q2)) > 0.5d-14)) &
  !            write(6,'(4(i4,2x,A),2(1x,ES17.10),1x,A,2(1x,ES17.10),1x,A)') k1,'|',q1,'|',k2,'|',q2,'|',JLinG(lp,k1,q1,k2,q2), &
  !                                                                          '|',JDipG(lp,k1,q1,k2,q2),'|'
  !        end do
  !      end do
  !    end do
  !  end do
  !  write(6,'(A)') '------|------|------|------|-------------------------------------|-------------------------------------|'
  !end if
  !
  !write(6,'(A)')
  !write(6,'(A)') 'Cartesian representation of the (rank-1)*(rank-1) exchange interaction'
  !Write(6,'(A)') 'Anisotropic exchange interaction:  J matrix:'
  !Write(6,'(A)') '     (  xx   xy  xz  )  '
  !Write(6,'(A)') 'J =  (  yx   yy  yz  )  '
  !Write(6,'(A)') '     (  zx   zy  zz  )  '
  !
  ! express the rank-1 tensors in the sum of
  ! Isotrop part    C * unit matrix
  ! Symmetric part
  ! Antisymmetric part
  ! print out the data:
  !
  !JDipCG(lp,:,:) = 0.0_wp
  !if(Lines .or. AnisoLines) then
  !  write(6,'(A)') 'Lines matrix:'
  !  do i=1,3
  !    write(6,'(3F12.6)') (JLinCG(lp,i,j),j=1,3)
  !  end do
  !
  !  ELin(lp,:,:) = 0.0_wp
  !  ALin(lp,:,:) = 0.0_wp
  !  SLin(lp,:,:) = 0.0_wp
  !  ELin(lp,1,1) = (JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/3.0_wp
  !  ELin(lp,2,2) = (JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/3.0_wp
  !  ELin(lp,3,3) = (JLinCG(lp,1,1)+JLinCG(lp,2,2)+JLinCG(lp,3,3))/3.0_wp
  !  do is1=1,3
  !    do is2=1,3
  !      ALin(lp,is1,is2) = (JLinCG(lp,is1,is2)-JLinCG(lp,is2,is1))/2.0_wp
  !    end do
  !  end do
  !  tsum = 0.0_wp
  !  do l=1,3
  !    tsum = tsum+JLinCG(lp,l,l)
  !  end do
  !  test = 0.0_wp
  !  do is1=1,3
  !    test(is1,is1) = 2.0_wp*tsum/3.0_wp
  !  end do
  !
  !  do is1=1,3
  !    do is2=1,3
  !      SLin(lp,is1,is2) = (JLinCG(lp,is1,is2)+JLinCG(lp,is2,is1)-test(is1,is2))/2.0_wp
  !    end do
  !  end do
  !  write(6,'(A)') 'Elin * unit'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (ELin(lp,is1,is2),is2=1,3)
  !  end do
  !  write(6,'(A)') '-2/3 * test'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (-2.0_wp*test(is1,is2)/3.0_wp,is2=1,3)
  !  end do
  !
  !  write(6,'(A)') 'ALin:'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (ALin(lp,is1,is2),is2=1,3)
  !  end do
  !  write(6,'(A)') 'SLin:'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (SLin(lp,is1,is2),is2=1,3)
  !  end do
  !end if
  !
  !if (Dipol) then
  !  write(6,'(A)') 'Dipolar exchange matrix:'
  !  do i=1,3
  !    write(6,'(3F12.6)') (JDipCG(lp,i,j),j=1,3)
  !  end do
  !
  !  EDip(lp,:,:) = 0.0_wp
  !  ADip(lp,:,:) = 0.0_wp
  !  SDip(lp,:,:) = 0.0_wp
  !  EDip(lp,1,1) = (JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/3.0_wp
  !  EDip(lp,2,2) = (JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/3.0_wp
  !  EDip(lp,3,3) = (JDipCG(lp,1,1)+JDipCG(lp,2,2)+JDipCG(lp,3,3))/3.0_wp
  !  do is1=1,3
  !    do is2=1,3
  !      ADip(lp,is1,is2) = (JDipCG(lp,is1,is2)-JDipCG(lp,is2,is1))/2.0_wp
  !    end Do
  !  end Do
  !  tsum = 0.0_wp
  !  do l=1,3
  !    tsum = tsum+JDipCG(lp,l,l)
  !  end do
  !  test = 0.0_wp
  !  do is1=1,3
  !    test(is1,is1) = 2.0_wp*tsum/3.0_wp
  !  end do
  !  do is1=1,3
  !    do is2=1,3
  !      SDip(lp,is1,is2) = (JDipCG(lp,is1,is2)+JDipCG(lp,is2,is1)-test(is1,is2))/2.0_wp
  !    end do
  !  end do
  !  write(6,'(A)') 'EDip * unit'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (EDip(lp,is1,is2),is2=1,3)
  !  end do
  !  write(6,'(A)') '-2/3 * test'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (-2.0_wp*test(is1,is2)/3.0_wp,is2=1,3)
  !  end do
  !  write(6,'(A)') 'ADip:'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (ADip(lp,is1,is2),is2=1,3)
  !  end do
  !  write(6,'(A)') 'SDip:'
  !  do is1=1,3
  !    write(6,'(3F12.6)') (SDip(lp,is1,is2),is2=1,3)
  !  end do
  !end if

  call mma_deallocate(JN)
  call mma_deallocate(JB)
  call mma_deallocate(JS)
end do ! lp, interacting pairs

return

end subroutine pr_ito_int
