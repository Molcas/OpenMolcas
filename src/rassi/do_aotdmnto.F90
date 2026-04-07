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
! Copyright (C) 2021, Rulin Feng                                       *
!***********************************************************************

!****************************************************
!           Do SVD for SO-TDM in AO basis
!****************************************************
! This routine is made to do single value decompositon
! to spin-orbit coupled transition density matrices in AO
! basis.
! Input: TDMZZ, transition densitry matrix
!        TSDMZZ,transition spin density matrix(x,y,z)
! Remember that TDMZZ and TSDMZZ contain six TDM's each
! See sonatorbm_full
! TDMZZ(3,:) and TSDMZZ(1-3,:), for the real part of TDM
! The transition density matrix TDMZZ does not depend on
! spin matrices, thus TDMZZ(1-3,:) are the same, so are the
! imaginary part TDMZZ(4-6,:).
! TDMZZ(6,:) and TSDMZZ(4-6,:), for the imaginary part of TDM
!
!                                                -RF 8/18,2021
subroutine DO_AOTDMNTO(TDMZZ,TSDMZZ,ANTSIN,ISTATE,JSTATE,nb,nb2)

use OneDat, only: sNoNuc, sNoOri, sOpSiz
use Cntrl, only: IfArgu
use rassi_data, only: NBTRI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Half, Pi, cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) ISTATE, JSTATE, nb, nb2
real(kind=wp) :: TDMZZ(6,nb2), TSDMZZ(6,nb2), ANTSIN(6,nb2)
integer(kind=iwp) :: di, i, icmp, iDummy(7,8), info, iopt, irc, isylab, j, LU, lwork, SIZ(1)
real(kind=wp) :: Dummy(1), NumofEc, phi, sd, Sumofeigen, SumofTDMZZLC, ttdi(3), ttdr(3)
complex(kind=wp) :: Transition_Dipole
integer(kind=iwp), allocatable :: PIV(:)
real(kind=wp), allocatable :: BFF(:), Dip(:), Dips(:), EIG(:), EIGM(:), RESIR(:), RESIX(:), SM(:), SMI(:), SVDS(:), SVDUI(:), &
                              SVDUR(:), SVDVHI(:), SVDVHR(:), SVDVI(:), SVDVR(:), SZZ(:), SZZs(:), TDMZZL(:,:), TMP(:), TMPI(:), &
                              TMPR(:), TSDMZZL(:,:)
complex(kind=wp), allocatable :: BUFF(:), BUFF1(:), BUFF2(:), DIPsC(:), RESI(:), SIZC(:), SumofYdiag(:), SVDU(:), SVDVH(:), &
                                 TDMZZC(:), TDMZZLC(:), YMAT(:,:)
character(len=128) :: FNAME
character(len=72) :: NOTE
character(len=8) :: LABEL
character(len=7) :: STATENAME, STATENAMETMP
real(kind=wp), parameter :: eigen_print_limit = 1.0e-8_wp
integer(kind=iwp), external :: isfreeunit

! trace of transition dipole real and imaginary (x,y,and z)
! ttdr, ttdi

! ANTISYMMETRIC matrix needs a little fixing
do i=1,nb
  do j=1,nb
    if (i < j) then
      ANTSIN(3,(i-1)*nb+j) = -ANTSIN(3,(i-1)*nb+j)
    else if (i == j) then
      ANTSIN(3,(i-1)*nb+j) = Zero
    end if
  end do
end do
! The imaginary part may need a negative sign
TDMZZ(4,:) = -TDMZZ(4,:)
TDMZZ(5,:) = -TDMZZ(5,:)
TDMZZ(6,:) = -TDMZZ(6,:)
! 'HERMSING' ITYPE=1
! 'ANTISING' ITYPE=2
! 'HERMTRIP' ITYPE=3
! 'ANTITRIP' ITYPE=4

! Thus we obtained the AO based transition density matrix TDMZZ
! and the transition spin density matrix TSDMZZ
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! we can do testing here we calculate the oscillator strength
! by dot with dipole moment AO matrix and trace
! do the test before the Lowdin orthogonalization
!if (TestPrint) then
call MMA_ALLOCATE(BUFF,nb2,LABEL='LBUFF')
call MMA_ALLOCATE(TDMZZC,nb2,LABEL='TDMZZC')
do di=1,3
  call mma_allocate(DIPs,nb2,Label='Dips')
  LABEL(1:8) = 'MLTPL  1'
  IRC = -1
  ICMP = di
  ISYLAB = 1
  IOPT = ibset(0,sOpSiz) ! Only read the size of the array
  call IRDONE(IRC,IOPT,LABEL,ICMP,SIZ,ISYLAB)
  !no nuclear contrib, no origin of operator
  IOPT = ibset(ibset(0,sNoOri),sNoNuc)
  call mma_allocate(DIP,SIZ(1),Label='DIP')
  call RDONE(IRC,IOPT,LABEL,ICMP,DIP,ISYLAB)
  call DESYM_SONTO(DIP,SIZ(1),DIPs,ISYLAB)
  call mma_deallocate(DIP)
  write(u6,*) '  For istate ',ISTATE,' jstate ',JSTATE
  write(u6,*) '  Component ',ICMP
  ! Get complex matrices
  call MMA_ALLOCATE(DIPsC,nb2,LABEL='DIPsC')
  TDMZZC(:) = cmplx(TDMZZ(di,:),TDMZZ(di+3,:),kind=wp)
  DIPsC(:) = DIPs(:)*cOne
  ! TDM
  call ZGEMM_('N','N',nb,nb,nb,cOne,DIPsC,nb,TDMZZC,nb,cZero,BUFF,nb)
  ! Trace the resulting matrix
  Transition_Dipole = cZero
  do i=1,nb
    Transition_Dipole = Transition_Dipole+BUFF((i-1)*nb+i)
  end do
  ttdr(di) = real(Transition_Dipole)
  ttdi(di) = aimag(Transition_Dipole)
  write(u6,*) '  Transition_Dipole :',Transition_Dipole
  write(u6,*)
  call MMA_DEALLOCATE(DIPs)
  call MMA_DEALLOCATE(DIPsC)
end do
call MMA_DEALLOCATE(BUFF)
call MMA_DEALLOCATE(TDMZZC)
!end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! phase factor = cos phi + i * sin phi
! phi = 1/2 * arc tan (2*(x_r*x_i + y_r*y_i + z_r*z_i)/
!                        (x_i**2 + y_i**2 + z_i**2 -
!                         x_r**2 - y_r**2 - z_r**2))
! to minimize the imaginary part of total transition dipole
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if (IFARGU) then
  ! check if minimum of maximum
  if ((ttdi(1)**2+ttdi(2)**2+ttdi(3)**2-ttdr(1)**2-ttdr(2)**2-ttdr(3)**2) == Zero) then
    phi = Zero
  else
    phi = Half*atan(Two*(ttdr(1)*ttdi(1)+ttdr(2)*ttdi(2)+ttdr(3)*ttdi(3))/ &
                    (ttdi(1)**2+ttdi(2)**2+ttdi(3)**2-ttdr(1)**2-ttdr(2)**2-ttdr(3)**2))
  end if
  sd = Two*cos(Two*phi)*(ttdr(1)**2+ttdr(2)**2+ttdr(3)**2-ttdi(1)**2-ttdi(2)**2-ttdi(3)**2)- &
       Four*sin(Two*phi)*(ttdr(1)*ttdi(1)+ttdr(2)*ttdi(2)+ttdr(3)*ttdi(3))
  ! make sure it's minimum
  if (sd < Zero) phi = phi+Pi/Two
  ! multipole phase factor with tdm as a whole
  write(u6,*) 'Phase factor turned on with calculated'
  write(u6,'(2X,A,F6.2)') 'argument Phi: ',phi
  call mma_allocate(TMPR,nb2,Label='TMPR')
  call mma_allocate(TMPI,nb2,Label='TMPI')
  TMPR(:) = TDMZZ(3,:)*cos(phi)-TDMZZ(6,:)*sin(phi)
  TMPI(:) = TDMZZ(6,:)*cos(phi)+TDMZZ(3,:)*sin(phi)
  TDMZZ(1,:) = TMPR(:)
  TDMZZ(2,:) = TMPR(:)
  TDMZZ(3,:) = TMPR(:)
  TDMZZ(4,:) = TMPI(:)
  TDMZZ(5,:) = TMPI(:)
  TDMZZ(6,:) = TMPI(:)
  call mma_deallocate(TMPI)
  call mma_deallocate(TMPR)
end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Make NTO output file names without spin
write(STATENAME,'(I3)') ISTATE
write(STATENAMETMP,'(I3,a1,a)') JSTATE,'_',trim(adjustl(STATENAME))
write(STATENAME,'(a)') trim(adjustl(STATENAMETMP))
! Everything is in C1 symmetry for now
! Do the Lowdin Orthogonalization assuming C1 symmetry
! SZZ  - AO Overlap integral
! SZZs - AO Overlap integral in square
! EIG  - AO Overlap eigenvalues
call mma_allocate(SZZ,NBTRI,Label='SZZ')
call mma_allocate(SZZs,nb2,Label='SZZs')
call mma_allocate(EIG,nb,Label='EIG')
SZZ(:) = Zero
SZZs(:) = Zero
EIG(:) = Zero
! AO OVERLAP MATRIX
IRC = -1
! IOPT=6, origin and nuclear contrib not read
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
ICMP = 1
ISYLAB = 1
LABEL = 'MLTPL  0'
call RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
if (IRC /= 0) then
  write(u6,*)
  write(u6,*) '      *** ERROR IN SUBROUTINE  SONATORB ***'
  write(u6,*) '      OVERLAP INTEGRALS ARE NOT AVAILABLE'
  write(u6,*)
  call ABEND()
end if
call DESYM_SONTO(SZZ,NBTRI,SZZs,ISYLAB)
call mma_deallocate(SZZ)
!*************
! For tests
!*************
call mma_allocate(TMP,nb2,Label='TMP')
TMP(:) = TDMZZ(3,:)
call mma_allocate(BFF,nb2,Label='BFF')
call DGEMM_('N','N',nb,nb,nb,One,SZZs,nb,TMP,nb,Zero,BFF,nb)
! Trace the resulting matrix
NumOfEc = Zero
do i=1,nb
  NumOfEc = NumOfEc+BFF((i-1)*nb+i)
end do
!write(u6,*) 'NumOfEc ',NumOfEc
call mma_deallocate(BFF)
!*************

! DIAGONALIZE AO OVERLAP MATRIX
! Set LWORK=-1 to get the optimal scratch space RESI
! then let LWORK equal to length of scratch space
! free and reallocate memory for RESI using that length
call mma_allocate(RESIX,1,Label='RESIX')
LWORK = -1
call DSYEV_('V','U',nb,SZZs,nb,EIG,RESIX,LWORK,INFO)
LWORK = int(RESIX(1))
call mma_deallocate(RESIX)
call mma_allocate(RESIX,LWORK,Label='RESIX')
! SZZs in as the AO overlap sqaure matrix
! out as the eigenvector matrix of SZZs
! with eigenvalues in EIG
call DSYEV_('V','U',nb,SZZs,nb,EIG,RESIX,LWORK,INFO)
! Put EIG in sqrt and in diagonal in EIGM
call mma_allocate(EIGM,nb2,Label='EIGM')
EIGM(:) = Zero
do i=1,nb
  EIGM((i-1)*nb+i) = sqrt(EIG(i))
end do
call mma_deallocate(EIG)
call mma_deallocate(RESIX)
! Get S^1/2 from S^1/2 = U S_diag^1/2 U^T
call mma_allocate(SM,nb2,Label='SM')
call DGEMM_('N','T',nb,nb,nb,One,EIGM,nb,SZZs,nb,Zero,TMP,nb)
call DGEMM_('N','N',nb,nb,nb,One,SZZs,nb,TMP,nb,Zero,SM,nb)
call mma_deallocate(SZZs)
call mma_deallocate(EIGM)
! Get inverse of S^1/2 -> S^-1/2
! Before calling DGETRI, call DGETRF to factorize SM
! Set LWORK=-1 to get the optimal scratch space RESI
! then let LWORK equal to length of scratch space
! free and reallocate memory for RESI using that length
call mma_allocate(SMI,nb2,Label='SMI')
SMI(:) = SM(:)
call mma_allocate(PIV,nb,Label='PIV')
call DGETRF_(nb,nb,SMI,nb,PIV,INFO)
call mma_allocate(RESIX,1,Label='RESIX')
LWORK = -1
call DGETRI_(nb,SMI,nb,PIV,RESIX,LWORK,INFO)
LWORK = int(RESIX(1))
call mma_deallocate(RESIX)
call mma_allocate(RESIX,LWORK,Label='RESIX')
call DGETRI_(nb,SMI,nb,PIV,RESIX,LWORK,INFO)
call mma_deallocate(PIV)
call mma_deallocate(RESIX)

! Note: The density matrix should transform as S^1/2 D S^1/2
! Transform TDMZZ and TSDMZZ as S^1/2 T S^1/2
call MMA_ALLOCATE(TDMZZL,6,nb2,LABEL='LTDMZZL')
call MMA_ALLOCATE(TSDMZZL,6,nb2,LABEL='LTSDMZZL')
call mma_allocate(TMPR,nb2,Label='TMPR')
! Real part of TDMZZ
TMPR(:) = TDMZZ(3,:)
call DGEMM_('N','N',nb,nb,nb,One,TMPR,nb,SM,nb,Zero,TMP,nb)
call DGEMM_('N','N',nb,nb,nb,One,SM,nb,TMP,nb,Zero,TMPR,nb)
TDMZZL(3,:) = TMPR(:)
! Imaginary part of TDMZZ
TMPR(:) = TDMZZ(6,:)
call DGEMM_('N','N',nb,nb,nb,One,TMPR,nb,SM,nb,Zero,TMP,nb)
call DGEMM_('N','N',nb,nb,nb,One,SM,nb,TMP,nb,Zero,TMPR,nb)
TDMZZL(6,:) = TMPR(:)
! Do for all components of TSDMZZ
do i=1,6
  TMPR(:) = TSDMZZ(i,:)
  call DGEMM_('N','N',nb,nb,nb,One,TMPR,nb,SM,nb,Zero,TMP,nb)
  call DGEMM_('N','N',nb,nb,nb,One,SM,nb,TMP,nb,Zero,TMPR,nb)
  TSDMZZL(i,:) = TMPR(:)
end do
call mma_deallocate(TMPR)
! End of the Lowdin Orthogonalization

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Do SVD for transition density matrix as a whole
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! In order to use zgesvd first we combine TDMZZL(3,:)
! and TDMZZL(6,:) as a complex matrix
call MMA_ALLOCATE(TDMZZLC,nb2,LABEL='LTDMZZLC')
SumofTDMZZLC = Zero
TDMZZLC(:) = cmplx(TDMZZL(3,:),TDMZZL(6,:),kind=wp)
SumofTDMZZLC = sum(abs(TDMZZLC(:)))
! Do SVD by using ZGESVD, see lapack for documentation
! A = U * SIGMA * V^dagger
! Get work space for U, SIGMA, and V^dagger, VH
call MMA_ALLOCATE(SVDU,nb2,LABEL='SVDU')
call MMA_ALLOCATE(SVDS,nb,LABEL='SVDS')
call MMA_ALLOCATE(SVDVH,nb2,LABEL='SVDVH')
call MMA_ALLOCATE(SIZC,1,LABEL='RESI')
call mma_allocate(RESIR,5*nb,Label='RESIR')
! Set LWORK=-1 to get the optimal scratch space in SIZC
! then let LWORK equal to length of scratch space
! free and reallocate memory for SIZC using that length
LWORK = -1
call ZGESVD_('A','A',NB,NB,TDMZZLC,NB,SVDS,SVDU,NB,SVDVH,NB,SIZC,LWORK,RESIR,INFO)
LWORK = max(1,int(SIZC(1)))
call MMA_DEALLOCATE(SIZC)
call MMA_ALLOCATE(RESI,LWORK,LABEL='RESI')
! Do SVD for TDMZZLC
call ZCOPY_(nb2,[cZero],0,SVDU,1)
call DCOPY_(nb,[Zero],0,SVDS,1)
call ZCOPY_(nb2,[cZero],0,SVDVH,1)
call ZCOPY_(LWORK,[cZero],0,RESI,1)
if (SumofTDMZZLC >= 1.0e-20_wp) call ZGESVD_('A','A',nb,nb,TDMZZLC,nb,SVDS,SVDU,nb,SVDVH,nb,RESI,LWORK,RESIR,INFO)
if (INFO /= Zero) write(u6,*) 'SVD convergence issue'
! End testing SVD
! For partitioning properties in the NTO basis
! Partition of the MLTPL 1, dipole moment intergals
call mma_allocate(DIPs,nb2,Label='DIPs')
call MMA_ALLOCATE(DIPsC,nb2,LABEL='LDIPsC')
call MMA_ALLOCATE(BUFF1,nb2,LABEL='BUFF1')
call MMA_ALLOCATE(BUFF2,nb2,LABEL='BUFF2')
call MMA_ALLOCATE(YMAT,3,nb2,LABEL='YMAT')
call MMA_ALLOCATE(SumofYdiag,3,LABEL='SumofYdiag')
! The three components of dipole
do di=1,3
  LABEL = 'MLTPL  1'
  IRC = -1
  ICMP = di
  ISYLAB = 1
  IOPT = ibset(0,sOpSiz)
  call IRDONE(IRC,IOPT,LABEL,ICMP,SIZ,ISYLAB)
  IOPT = 6
  call mma_allocate(DIP,SIZ(1),Label='DIP')
  call RDONE(IRC,IOPT,LABEL,ICMP,DIP,ISYLAB)
  call DESYM_SONTO(DIP,SIZ(1),DIPs,ISYLAB)
  call mma_deallocate(DIP)
  ! Perform Lowdin orthogonalization on operator matrix
  ! They transform as S^-1/2 P S^-1/2
  call DGEMM_('N','N',nb,nb,nb,One,DIPs,nb,SMI,nb,Zero,TMP,nb)
  call DGEMM_('N','N',nb,nb,nb,One,SMI,nb,TMP,nb,Zero,DIPs,nb)
  ! ZGESVD destroys TDMZZLC after it finishes
  ! reconstruct TDMZZLC and DIPsC
  do i=1,nb2
    TDMZZLC(i) = cmplx(TDMZZL(3,i),TDMZZL(6,i),kind=wp)
    DIPsC(i) = cmplx(DIPs(i),Zero,kind=wp)
  end do
  TDMZZLC(:) = cmplx(TDMZZL(3,:),TDMZZL(6,:),kind=wp)
  DIPsC(:) = DIPs(:)*cOne
  ! Do U^H TDMZZLC DIP U = Y, Diagonal of Y contains the partition
  call ZGEMM_('N','N',nb,nb,nb,cOne,TDMZZLC(:),nb,DIPsC,nb,cZero,BUFF1(:),nb)
  call ZGEMM_('N','N',nb,nb,nb,cOne,BUFF1(:),nb,SVDU(:),nb,cZero,BUFF2(:),nb)
  call ZGEMM_('C','N',nb,nb,nb,cOne,SVDU(:),nb,BUFF2(:),nb,cZero,BUFF1(:),nb)
  YMAT(di,:) = BUFF1(:)
  SumofYdiag(di) = cZero
  do i=1,nb
    SumofYdiag(di) = SumofYdiag(di)+YMAT(di,(i-1)*nb+i)
  end do
end do
call mma_deallocate(DIPs)
call MMA_DEALLOCATE(BUFF1)
call MMA_DEALLOCATE(BUFF2)
call MMA_DEALLOCATE(DIPsC)
call MMA_DEALLOCATE(RESI)
call MMA_DEALLOCATE(RESIR)

!call MMA_DEALLOCATE(YMAT)

!write(u6,'(F11.5,SP,F8.5,"i")') SumofYdiag(1)
!write(u6,'(F11.5,SP,F8.5,"i")') SumofYdiag(2)
!write(u6,'(F11.5,SP,F8.5,"i")') SumofYdiag(3)
! End of partitioning properties

! But we still need to transform U and V back to original AO
! basis using U=S^{-1/2) U' and V=S^{-1/2} V'
! for V^t it is V'^t S^{-1/2} = V^t
call mma_allocate(SVDUR,nb2,Label='SVDUR')
call mma_allocate(SVDUI,nb2,Label='SVDUI')
SVDUR(:) = Zero
SVDUI(:) = Zero
call mma_allocate(SVDVHR,nb2,Label='SVDVHR')
call mma_allocate(SVDVHI,nb2,Label='SVDVHR')
SVDVHR(:) = Zero
SVDVHI(:) = Zero

TMP(:) = Zero

SVDUR(:) = real(SVDU(:))
SVDUI(:) = aimag(SVDU(:))
SVDVHR(:) = real(SVDVH(:))
SVDVHI(:) = aimag(SVDVH(:))
! U
call DGEMM_('N','N',nb,nb,nb,One,SMI,nb,SVDUR,nb,Zero,TMP,nb)
SVDUR(:) = TMP(:)
call DGEMM_('N','N',nb,nb,nb,One,SMI,nb,SVDUI,nb,Zero,TMP,nb)
SVDUI(:) = TMP(:)
! V^H
call DGEMM_('N','N',nb,nb,nb,One,SVDVHR,nb,SMI,nb,Zero,TMP,nb)
SVDVHR(:) = TMP(:)
call DGEMM_('N','N',nb,nb,nb,One,SVDVHI,nb,SMI,nb,Zero,TMP,nb)
SVDVHI(:) = TMP(:)
! V
call mma_allocate(SVDVR,nb2,Label='SVDVR')
call mma_allocate(SVDVI,nb2,Label='SVDVI')
do i=1,nb
  do j=1,nb
    SVDVR((i-1)*nb+j) = SVDVHR((j-1)*nb+i)
    ! imaginary part takes a negative sign
    SVDVI((i-1)*nb+j) = -SVDVHI((j-1)*nb+i)
  end do
end do
call mma_deallocate(SVDVHR)
call mma_deallocate(SVDVHI)
! tests
call ADD_INFO('LAMBDA',SVDS,5,4)

! singular values
Sumofeigen = sum(SVDS(:)**2)

! Head of the output
write(u6,*)
write(u6,'(6X,A)') repeat('*',90)
write(u6,'(6X,A,88X,A)') '*','*'
write(u6,'(6X,A,29X,A31,28X,A)') '*','Natural transition orbitals','*'
write(u6,'(6X,A,88X,A)') '*','*'
write(u6,'(6X,A,27X,A25,I2,A12,I2,20X,A)') '*','Between spin-orbit state ',ISTATE,' and state ',JSTATE,'*'
write(u6,'(6X,A,88X,A)') '*','*'
write(u6,'(6X,A)') repeat('*',90)
write(u6,*)
! Start output singular value information for positive spin values
write(u6,'(6X,A)') repeat('=',90)
write(u6,'(5X,A12,A12,A16,A51)') 'EXCITATION','EIGENVALUE','EXCITATION','TRANSITION DIPOLE MOMENT'
write(u6,'(5X,A12,12X,A16,3A17)') 'AMPLITUDE','CONTRIBUTION(%)','(1)','(2)','(3)'
write(u6,'(6X,A)') repeat('-',90)
do i=1,nb
  if (SVDS(i)**2 < eigen_print_limit) exit
  write(u6,'(4X,3X,F8.5,4X,F8.5,8X,F8.2,2X,3(F9.4,SP,F7.4,"i",SS))') SVDS(i),SVDS(i)**2,SVDS(i)**2/Sumofeigen*100.0_wp, &
                                                                     YMAT(1,(i-1)*nb+i),YMAT(2,(i-1)*nb+i),YMAT(3,(i-1)*nb+i)
end do
write(u6,'(6X,A,F8.5)') 'SUM OF EIGENVALUES ',Sumofeigen
write(u6,'(6X,A24,15X,3(F9.4,SP,F7.4,"i",SS))') 'SUM OF TRANSITION DIPOLE',SumofYdiag(1),SumofYdiag(2),SumofYdiag(3)
write(u6,'(6X,A)') repeat('=',90)
write(u6,*)
write(u6,*)
! Write NTOs to file in C1 symmetry
SVDS(:) = SVDS(:)**2/Sumofeigen
LU = ISFREEUNIT(50)
Note = '*  Spin-orbit Natural Transition Orbitals'
! U real
write(FNAME,'(6(a))') 'NTORB.SO.',trim(adjustl(STATENAME)),'.','PART','.','Re'
write(u6,'(4(a))') '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',trim(STATENAME),' ARE WRITTEN ONTO FILE ',FNAME
call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDUR,SVDS,Dummy,iDummy,Note)
! U imaginary
write(FNAME,'(6(a))') 'NTORB.SO.',trim(adjustl(STATENAME)),'.','PART','.','Im'
write(u6,'(4(a))') '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',trim(STATENAME),' ARE WRITTEN ONTO FILE ',FNAME
call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDUI,SVDS,Dummy,iDummy,Note)
! V real
write(FNAME,'(6(a))') 'NTORB.SO.',trim(adjustl(STATENAME)),'.','HOLE','.','Re'
write(u6,'(4(a))') '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',trim(STATENAME),' ARE WRITTEN ONTO FILE ',FNAME
call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDVR,SVDS,Dummy,iDummy,Note)
! V imaginary
write(FNAME,'(6(a))') 'NTORB.SO.',trim(adjustl(STATENAME)),'.','HOLE','.','Im'
write(u6,'(4(a))') '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',trim(STATENAME),' ARE WRITTEN ONTO FILE ',FNAME
call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDVI,SVDS,Dummy,iDummy,Note)
! End of output
write(u6,*)
write(u6,'(6X,A)') repeat('*',90)
write(u6,'(6X,A,88X,A)') '*','*'
write(u6,'(6X,A,28X,A34,25X,A)') '*','End of natural transition orbitals','*'
write(u6,'(6X,A,88X,A)') '*','*'
write(u6,'(6X,A,27X,A25,I2,A12,I2,20X,A)') '*','Between spin-orbit state ',ISTATE,' and state ',JSTATE,'*'
write(u6,'(6X,A,88X,A)') '*','*'
write(u6,'(6X,A)') repeat('*',90)
write(u6,*)

! Free up workspace
call MMA_DEALLOCATE(SVDU)
call MMA_DEALLOCATE(SVDS)
call MMA_DEALLOCATE(SVDVH)
call MMA_DEALLOCATE(TDMZZL)
call MMA_DEALLOCATE(TSDMZZL)
call MMA_DEALLOCATE(TDMZZLC)
call MMA_DEALLOCATE(YMAT)
call MMA_DEALLOCATE(SumofYdiag)
call MMA_DEALLOCATE(SM)
call MMA_DEALLOCATE(SMI)
call MMA_DEALLOCATE(TMP)
call MMA_DEALLOCATE(SVDUR)
call MMA_DEALLOCATE(SVDUI)
call MMA_DEALLOCATE(SVDVR)
call MMA_DEALLOCATE(SVDVI)

end subroutine DO_AOTDMNTO
