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

subroutine Kinetic_Exchange(N1,N2,M1,S1,M2,S2,eso1,eso2,tpar,upar,lant,OPT,HKEX,MR1,SR1,MR2,SR2)
! compute KE, within various options :
! the Ln site
!  N1, eso1, M1, S1, MR1, SR1
! the radical
!  N2, eso2, M2, S2, MR2, SR2
! exchange Hamiltonian
!  HKEX

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N1, N2, lant, OPT
complex(kind=wp), intent(in) :: M1(3,N1,N1), S1(3,N1,N1), M2(3,N2,N2), S2(3,N2,N2)
real(kind=wp), intent(in) :: eso1(N1), eso2(N2), tpar, upar
complex(kind=wp), intent(out) :: HKEX(N1,N1,N2,N2), MR1(3,N1,N1), SR1(3,N1,N1), MR2(3,N2,N2), SR2(3,N2,N2)
integer(kind=iwp) :: i, i1, i2, info, iprint, is1, is2, j, l
real(kind=wp) :: gtens(3,4), maxes(3,3,4)
complex(kind=wp) :: MM(3,2,2)
real(kind=wp), allocatable :: eloc1(:), eloc2(:), wcr(:)
complex(kind=wp), allocatable :: ABIT(:,:,:,:), H1(:,:), H1T(:,:), H2(:,:), HCOV(:,:), HEXC(:,:,:,:), MM1(:,:,:), SM1(:,:,:), &
                                 TMP(:), TMP2(:,:), TMP3(:,:), Z1(:,:), Z2(:,:), ZCR(:,:), ZZ1(:,:), ZZ2(:,:)

! determine the pseudospin on each site (Z1 and Z2):
iprint = 1
gtens(:,:) = 0.0_wp
maxes(:,:,:) = 0.0_wp
call mma_allocate(Z1,N1,N1,label='Z1')
call mma_allocate(Z2,N2,N2,label='Z2')
call pseudospin(M1,N1,Z1,3,1,iprint)
call pseudospin(M2,N2,Z2,3,1,iprint)
#ifdef _DEBUGPRINT_
call pa_prMat('KE_Exchange:: Pseudospin site 1',Z1,N1)
call pa_prMat('KE_Exchange:: Pseudospin site 2',Z2,N2)
#endif
! get the "traced-to-zero" local energy states on both sites:
call mma_allocate(eloc1,N1,label='eloc1')
call mma_allocate(eloc2,N1,label='eloc2')
call rtrace(N1,eso1,eloc1)
call rtrace(N2,eso2,eloc2)
call mma_allocate(H1,N1,N1,label='H1')
call mma_allocate(H2,N2,N2,label='H2')
H1(:,:) = cZero
H2(:,:) = cZero
do I2=1,N1
  do i=1,N1
    H1(:,I2) = H1(:,I2)+ELOC1(i)*conjg(Z1(i,:))*Z1(i,I2)
  end do
end do
do I2=1,N2
  do i=1,N2
    H2(:,I2) = H2(:,I2)+ELOC2(i)*conjg(Z2(i,:))*Z2(i,I2)
  end do
end do
call mma_deallocate(eloc1)
call mma_deallocate(eloc2)
call mma_allocate(HCOV,N1,N1,label='HCOV')
call mma_allocate(HEXC,N1,N1,N2,N2,label='HEXC')
if ((OPT == 1) .or. (OPT == 3) .or. (OPT > 4)) then ! full
  call KE_Covalent(N1,lant,tpar,upar,1,HCOV)
  call KE_Exchange(N1,N2,lant,tpar,upar,1,HEXC)
else if ((OPT == 2) .or. (OPT == 4)) then ! 1/U model
  call KE_Covalent(N1,lant,tpar,upar,2,HCOV)
  call KE_Exchange(N1,N2,lant,tpar,upar,2,HEXC)
end if

call mma_allocate(H1T,N1,N1,label='H1T')
call mma_allocate(TMP,max(N1,N2)**2,label='TMP')
H1T(:,:) = H1(:,:)+HCOV(:,:)
! rewrite the HCOV in the initial ab initio basis:
call ZGEMM_('N','N',N1,N1,N1,cOne,Z1,N1,HCOV,N1,cZero,TMP,N1)
call ZGEMM_('N','C',N1,N1,N1,cOne,TMP,N1,Z1,N1,cZero,HCOV,N1)
! rewrite the H1 in the initial ab initio basis:
call ZGEMM_('N','N',N1,N1,N1,cOne,Z1,N1,H1,N1,cZero,TMP,N1)
call ZGEMM_('N','C',N1,N1,N1,cOne,TMP,N1,Z1,N1,cZero,H1,N1)
! rewrite the H2 in the initial ab initio basis:
call ZGEMM_('N','N',N2,N2,N2,cOne,Z2,N2,H2,N2,cZero,TMP,N2)
call ZGEMM_('N','C',N2,N2,N2,cOne,TMP,N2,Z2,N2,cZero,H2,N2)
! rewrite the H1T in the initial ab initio basis:
call ZGEMM_('N','N',N1,N1,N1,cOne,Z1,N1,H1T,N1,cZero,TMP,N1)
call ZGEMM_('N','C',N1,N1,N1,cOne,TMP,N1,Z1,N1,cZero,H1T,N1)
! rewrite the HEXC in the initial ab initio basis:
do is1=1,N2
  do is2=1,N2
    call ZGEMM_('N','N',N1,N1,N1,cOne,Z1,N1,HEXC(:,:,is1,is2),N1,cZero,TMP,N1)
    call ZGEMM_('N','C',N1,N1,N1,cOne,TMP,N1,Z1,N1,cZero,HEXC(:,:,is1,is2),N1)
  end do
end do
do is1=1,N1
  do is2=1,N1
    call ZGEMM_('N','N',N2,N2,N2,cOne,Z2,N2,HEXC(is1,is2,:,:),N2,cZero,TMP,N2)
    call ZGEMM_('N','C',N2,N2,N2,cOne,TMP,N2,Z2,N2,cZero,HEXC(is1,is2,:,:),N2)
  end do
end do

call mma_deallocate(Z1)
call mma_deallocate(Z2)

call mma_allocate(ABIT,N1,N1,N2,N2,label='ABIT')
HKEX(:,:,:,:) = cZero
do i=1,N1
  do i1=1,N2
    ABIT(i,i,i1,i1) = H1(i,i)+H2(i1,i1)
  end do
end do
HKEX(:,:,:,:) = HEXC(:,:,:,:)+ABIT(:,:,:,:)
call mma_deallocate(ABIT)
call mma_deallocate(H1)
call mma_deallocate(H2)
call mma_deallocate(HCOV)
call mma_deallocate(HEXC)
! H1T = H1 + HCOV
call mma_allocate(wcr,N1,label='wcr')
call mma_allocate(zcr,N1,N1,label='zcr')
call diag_c2(H1T,N1,info,wcr,zcr)
do i=1,N1
  write(u6,'(2(A,i2,A,F15.9))') 'ESO1(',i,')=',ESO1(i),'  ESO1+COV(',i,')=',wcr(i)-wcr(1)
end do
do i=1,N1
  do j=1,N1
    write(u6,'(A,i2,A,i2,A,2F20.14)') 'ZCR(',i,',',j,')=',ZCR(i,j)
  end do
end do
call mma_deallocate(wcr)
call mma_deallocate(H1T)

! rotate to COV basis:
if ((opt == 3) .or. (opt == 4)) then
  do is1=1,N2
    do is2=1,N2
      call ZGEMM_('C','N',N1,N1,N1,cOne,ZCR,N1,HKEX(:,:,is1,is2),N1,cZero,TMP,N1)
      call ZGEMM_('N','N',N1,N1,N1,cOne,TMP,N1,ZCR,N1,cZero,HKEX(:,:,is1,is2),N1)
    end do
  end do
end if
! compute the g tensors for initial and initial+covalence:
call mma_allocate(MM1,3,N1,N1,label='MM1')
call mma_allocate(SM1,3,N1,N1,label='SM1')
MM(:,:,:) = M1(:,1:2,1:2)
call atens(MM,2,gtens(:,1),maxes(:,:,1),1)
MM(:,:,:) = M1(:,3:4,3:4)
call atens(MM,2,gtens(:,2),maxes(:,:,2),1)
call mma_allocate(TMP2,N1,N1,label='TMP2')
do L=1,3
  TMP2(:,:) = M1(L,:,:)
  call ZGEMM_('C','N',N1,N1,N1,cOne,ZCR,N1,TMP2,N1,cZero,TMP,N1)
  call ZGEMM_('N','N',N1,N1,N1,cOne,TMP,N1,ZCR,N1,cZero,TMP2,N1)
  MM1(L,:,:) = TMP2(:,:)
  TMP2(:,:) = S1(L,:,:)
  call ZGEMM_('C','N',N1,N1,N1,cOne,ZCR,N1,TMP2,N1,cZero,TMP,N1)
  call ZGEMM_('N','N',N1,N1,N1,cOne,TMP,N1,ZCR,N1,cZero,TMP2,N1)
  SM1(L,:,:) = TMP2(:,:)
end do
call mma_deallocate(zcr)
call mma_deallocate(TMP2)
MM(:,:,:) = MM1(:,1:2,1:2)
call atens(MM,2,gtens(:,3),maxes(:,:,3),1)
MM(:,:,:) = MM1(:,3:4,3:4)
call atens(MM,2,gtens(:,4),maxes(:,:,4),1)
write(u6,'(A)') 'Initial g tensors of the ground and first excited KD'
do i=1,2
  write(u6,'((A,F12.6,A,3F12.7))') 'gX=',gtens(1,i),' axis X: ',(maxes(j,1,i),j=1,3)
  write(u6,'((A,F12.6,A,3F12.7))') 'gY=',gtens(2,i),' axis Y: ',(maxes(j,2,i),j=1,3)
  write(u6,'((A,F12.6,A,3F12.7))') 'gZ=',gtens(3,i),' axis Z: ',(maxes(j,3,i),j=1,3)
  write(u6,*)
end do
write(u6,'(A)') 'Initial+Covalence g tensors of the ground and first excited KD'
do i=3,4
  write(u6,'((A,F12.6,A,3F12.7))') 'gX=',gtens(1,i),' axis X: ',(maxes(j,1,i),j=1,3)
  write(u6,'((A,F12.6,A,3F12.7))') 'gY=',gtens(2,i),' axis Y: ',(maxes(j,2,i),j=1,3)
  write(u6,'((A,F12.6,A,3F12.7))') 'gZ=',gtens(3,i),' axis Z: ',(maxes(j,3,i),j=1,3)
  write(u6,*)
end do
!------
if ((opt == 3) .or. (opt == 4)) then
  ! rotate the magnetic moment to the coordinate
  ! system of main magnetic axes on Ln
  call rotmom2(SM1,N1,maxes(:,:,3),SR1)
  call rotmom2(MM1,N1,maxes(:,:,3),MR1)
  call rotmom2(S2,N2,maxes(:,:,3),SR2)
  call rotmom2(M2,N2,maxes(:,:,3),MR2)
  ! find the local pseudospins on both sites:
  call mma_allocate(ZZ1,N1,N1,label='ZZ1')
  call mma_allocate(ZZ2,N2,N2,label='ZZ2')
  call pseudospin(MR1,N1,ZZ1,3,1,iprint)
  call pseudospin(MR2,N2,ZZ2,3,1,iprint)
# ifdef _DEBUGPRINT_
  call pa_prMat('KE_Exchange:: Pseudospin site 1',ZZ1,N1)
  call pa_prMat('KE_Exchange:: Pseudospin site 2',ZZ2,N2)
# endif
  call mma_allocate(TMP2,N1,N1,label='TMP2')
  call mma_allocate(TMP3,N2,N2,label='TMP3')
  ! rewrite the magnetic moments and spin moments in new local bases:
  do L=1,3
    TMP2(:,:) = MR1(L,:,:)
    call ZGEMM_('C','N',N1,N1,N1,cOne,ZZ1,N1,TMP2,N1,cZero,TMP,N1)
    call ZGEMM_('N','N',N1,N1,N1,cOne,TMP,N1,ZZ1,N1,cZero,TMP2,N1)
    MR1(L,:,:) = TMP2(:,:)
    TMP2(:,:) = SR1(L,:,:)
    call ZGEMM_('C','N',N1,N1,N1,cOne,ZZ1,N1,TMP2,N1,cZero,TMP,N1)
    call ZGEMM_('N','N',N1,N1,N1,cOne,TMP,N1,ZZ1,N1,cZero,TMP2,N1)
    SR1(L,:,:) = TMP2(:,:)

    TMP3(:,:) = MR2(L,:,:)
    call ZGEMM_('C','N',N2,N2,N2,cOne,ZZ2,N2,TMP3,N2,cZero,TMP,N2)
    call ZGEMM_('N','N',N2,N2,N2,cOne,TMP,N2,ZZ2,N2,cZero,TMP3,N2)
    MR2(L,:,:) = TMP3(:,:)
    TMP3(:,:) = SR2(L,:,:)
    call ZGEMM_('C','N',N2,N2,N2,cOne,ZZ2,N2,TMP3,N2,cZero,TMP,N2)
    call ZGEMM_('N','N',N2,N2,N2,cOne,TMP,N2,ZZ2,N2,cZero,TMP3,N2)
    SR2(L,:,:) = TMP3(:,:)
  end do
  ! rewrite the exchnage matrix in the basis of local pseudospins:
  do is1=1,N2
    do is2=1,N2
      call ZGEMM_('C','N',N1,N1,N1,cOne,ZZ1,N1,HKEX(:,:,is1,is2),N1,cZero,TMP,N1)
      call ZGEMM_('N','N',N1,N1,N1,cOne,TMP,N1,ZZ1,N1,cZero,HKEX(:,:,is1,is2),N1)
    end do
  end do
  do is1=1,N1
    do is2=1,N1
      TMP3(:,:) = HKEX(is1,is2,:,:)
      call ZGEMM_('C','N',N2,N2,N2,cOne,ZZ2,N2,TMP3,N2,cZero,TMP,N2)
      call ZGEMM_('N','N',N2,N2,N2,cOne,TMP,N2,ZZ2,N2,cZero,TMP3,N2)
      HKEX(is1,is2,:,:) = TMP3(:,:)
    end do
  end do
  call mma_deallocate(ZZ1)
  call mma_deallocate(ZZ2)
  call mma_deallocate(TMP2)
  call mma_deallocate(TMP3)
else !opt=1, opt=2, and opt>4
  MR1(:,:,:) = M1(:,:,:)
  SR1(:,:,:) = S1(:,:,:)
  MR2(:,:,:) = M2(:,:,:)
  SR2(:,:,:) = S2(:,:,:)
end if

call mma_deallocate(MM1)
call mma_deallocate(SM1)
call mma_deallocate(TMP)

return

end subroutine Kinetic_Exchange
