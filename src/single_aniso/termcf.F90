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

subroutine termCF(ANGMOM,AMFI,ESFS,ldimcf,iDIM,maxes2,iopt,nlanth,iprint)
! this Subroutine calculates the parameters of term-specific crystal field
! for lanthanides:
!
!  nstate                      -- total number of spin free states
!                                 mixed in RASSI
!  lDIMcf                      -- number of orbital states in the
!                                 considered term (LS term)
!  ANGMOM ( 3, nstate, nstate) -- a real array holding the matrix
!                                 elements of L in this basis
!                                 (read from RASSI):
!  ESFS (nstate)               -- energy of the orbital states form the LS term
!  maxes2(3,3)                 -- rotation matrix needed to choose the
!                                 main quantization axis
!                                 ( determinant(maxes)=1.0_wp, orthogonal vectors )
!  iopt                        -- option for choosing the main quanization axis
!                   iopt = 1   -- axis of the ground orbital multiplet,
!                                 iDIM specifies the size of pseudo L
!                   iopt = 2   -- axis of the entire L manifold
!                   iopt = 3   -- maxes is defined by the user
!                                 ( maxes2 is given as input )
!                   iopt = 4   -- maxes is the unity matrix ( original Z
!                                 is the quantization axis )

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: ldimcf, iprint, iopt, nlanth, iDIM
real(kind=8), intent(in) :: esfs(ldimcf)
real(kind=8), intent(in) :: angmom(3,ldimcf,ldimcf)
real(kind=8), intent(in) :: amfi(3,ldimcf,ldimcf)
real(kind=8), intent(in) :: maxes2(3,3)
real(kind=8) :: finddetr
!real(kind=8) :: knm(12,0:12)
real(kind=8), allocatable :: maxes(:,:)
real(kind=8), allocatable :: gtens(:)
real(kind=8), allocatable :: eloc(:) ! lDIMcf
real(kind=8), allocatable :: Winit(:) ! lDIMcf
real(kind=8) :: BNC(lDIMcf,0:lDIMcf)
real(kind=8) :: BNS(lDIMcf,0:lDIMcf)
real(kind=8) :: Bstev(lDIMcf,-lDIMcf:lDIMcf)
complex(kind=8) :: Akq((lDIMcf-1),-(lDIMcf-1):(lDIMcf-1))
complex(kind=8), allocatable :: Angm(:,:,:) ! 3,ldimcf,ldimcf
complex(kind=8), allocatable :: dipso(:,:,:) ! 3,ldimcf,ldimcf
complex(kind=8), allocatable :: amfi_c(:,:,:) ! 3,ldimcf,ldimcf
complex(kind=8), allocatable :: amfi2(:,:,:) ! 3,ldimcf,ldimcf
complex(kind=8), allocatable :: amfi_l(:,:,:) ! 3,ldimcf,ldimcf
complex(kind=8), allocatable :: Z(:,:) !ldimcf,ldimcf
complex(kind=8), allocatable :: tmp(:,:) !ldimcf,ldimcf
complex(kind=8), allocatable :: HCF(:,:) !ldimcf,ldimcf
complex(kind=8), allocatable :: Zinit(:,:) !ldimcf,ldimcf
integer :: i, j, l, info, i1, i2
external :: finddetr
logical :: debug = .false.
real(kind=8) :: au2cm = 2.194746313705d5
real(kind=8) :: tS, tL, tJ, coeffCG, spinM, orbM, tJM, CF(100,100)
real(kind=8) :: det
integer :: MS, ML, MJ
integer :: ij, iLS, nLS, ibasS(100), ibasL(100), ibasJ(100)
integer :: irootL(100), ir, icas, k
complex(kind=8) :: CFC(100,100)

call mma_allocate(maxes,3,3,'maxes')
call mma_allocate(gtens,3,'gtens')
call mma_allocate(eloc,lDIMcf,'eloc')
call mma_allocate(Winit,lDIMcf,'Winit')

call mma_allocate(Angm,3,ldimcf,ldimcf,'angm')
call mma_allocate(dipso,3,ldimcf,ldimcf,'dipso')
call mma_allocate(amfi_c,3,ldimcf,ldimcf,'amfi_c')
call mma_allocate(amfi2,3,ldimcf,ldimcf,'amfi2')
call mma_allocate(amfi_l,3,ldimcf,ldimcf,'amfi_l')
call mma_allocate(Z,ldimcf,ldimcf,'Z')
call mma_allocate(Zinit,ldimcf,ldimcf,'Zinit')
call mma_allocate(tmp,ldimcf,ldimcf,'tmp')
call mma_allocate(HCF,ldimcf,ldimcf,'HCF')

call dcopy_(3,[0.0_wp],0,gtens,1)
call dcopy_(3*3,[0.0_wp],0,maxes,1)
call dcopy_(lDIMcf,[0.0_wp],0,eloc,1)
call dcopy_(lDIMcf,[0.0_wp],0,Winit,1)

call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,Angm,1)
call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,dipso,1)
call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi_c,1)
call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi2,1)
call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi_l,1)
call zcopy_(ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,Z,1)
call zcopy_(ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,Zinit,1)
call zcopy_(ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,tmp,1)
call zcopy_(ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,HCF,1)
!-----------------------------------------------------------------------
write(6,'(/)')
write(6,'(100A)') ('%',i=1,95)
if (mod(lDIMcf,2) == 1) then
  write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L = ',(lDIMcf-1)/2,'.'
else
  write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L = ',(lDIMcf-1),'/2.'
end if
write(6,'(100A)') ('%',i=1,95)
write(6,*)

do l=1,3
  do i=1,ldimcf
    do j=1,ldimcf
      dipso(l,i,j) = -cmplx(0.0_wp,angmom(l,i,j),wp)
      amfi_c(l,i,j) = cmplx(amfi(l,i,j),0.0_wp,wp)
    end do
  end do
end do
! find the main anisotropy direction of the
! diagonalize the angmom
if (iopt == 1) then
  ! iopt = 1   -- axis of the ground orbital Doublet
  call atens(dipso(1:3,1:iDIM,1:iDIM),iDIM,gtens,maxes,iprint)
  write(6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  if (mod(iDIM,2) == 0) then
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseudo-L = |',iDIM-1,'/2> orbital multiplet.'
  else
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseudo-L = |',(iDIM-1)/2,'> orbital multiplet.'
  end if

else if (iopt == 2) then
  ! iopt = 2   -- axis of the entire L manIfold
  call atens(dipso(1:3,1:ldimcf,1:ldimcf),lDIMcf,gtens,maxes,iprint)
  if (mod(lDIMCF,2) == 0) then
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic L = |',lDIMCF-1,'/2> multiplet'
  else
    write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic L = |',(lDIMCF-1)/2,'> multiplet'
  end if

else if (iopt == 3) then
  ! iopt = 3   -- maxes is defined by the user ( maxes is given as input )
  write(6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  write(6,'(a)') '(Xm, Ym, Zm) -- defined in the input file.'
  ! copy the maxes2 to maxes
  do i=1,3
    do j=1,3
      maxes(i,j) = maxes2(i,j)
    end do
  end do
  det = finddetr(maxes(1:3,1:3),3)
  if (det == 0.0_wp) then
    write(6,'(A)') 'TermCF:   iopt=3, while  DET(maxes)= 0.0'
    call AbEnd()
  end if
  ! FIXME: note finddetr modifies the matrix!
  if (det < 0.0_wp) then
    do i=1,3
      maxes(i,1) = -maxes2(i,1)
    end do
    if (iprint > 2) write(6,'(a)') 'The original coordinate system was LEFT-handed. It has been changed to the RIGHT-handed'
  end if

else if (iopt == 4) then
  ! iopt = 4   -- maxes is the unity matrix ( original Z is the quantization axis )
  do i=1,3
    maxes(i,i) = 1.0_wp
  end do
else
  call AbEnd()
end if

! print out the rotation matrix:
write(6,'(a)') 'Rotation matrix from the initial coordinate system to the employed coordinate system is:'

if ((iopt == 1) .or. (iopt == 2)) then
  write(6,'(70a)') ('-',i=1,67),'|'
  write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(6,'(A,35x,A)') 'Xm, Ym, Zm -- main magnetic axes','|'
  write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Xm |',(maxes(j,1),j=1,3),'|'
  write(6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(maxes(j,2),j=1,3),'|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Zm |',(maxes(j,3),j=1,3),'|'
  write(6,'(83a)') ('-',i=1,67),'|'
  write(6,'(A,I3)') 'Quantization axis is Zm.'

else if (iopt == 3) then

  write(6,'(70a)') ('-',i=1,67),'|'
  write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(6,'(A,11x,A)') 'Xm, Ym, Zm -- the coordinate system defined in the input','|'
  write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Xm |',(maxes(j,1),j=1,3),'|'
  write(6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(maxes(j,2),j=1,3),'|'
  write(6,'(6x,A,3F18.14,1x,A)') '| Zm |',(maxes(j,3),j=1,3),'|'
  write(6,'(83a)') ('-',i=1,67),'|'
  write(6,'(A,I3)') 'Quantization axis is Zm.'

else

  write(6,'(A)') 'IDENTITY matrix.'
  write(6,'(A,I3)') 'Quantization axis is the initial z axis.'
end if

! rotate the angular momentum to the new axes, using "maxes"
call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi2,1)
call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,angm,1)
call rotmom2(dipso(1:3,1:ldimcf,1:ldimcf),ldimcf,maxes,angm)
call rotmom2(amfi_c(1:3,1:ldimcf,1:ldimcf),ldimcf,maxes,amfi2)
if (debug) call prmom('TERMCF:: ANGM',angm,ldimcf)
if (debug) call prmom('TERMCF:: AMFI2',amfi2,ldimcf)

!-----------------------------------------------------------------------
! Find the pseudo-L basis of the LS manifold
call pseudospin(angm,ldimcf,Z,3,1,iprint)

if ((iprint >= 4) .or. debug) then
  write(6,*)
  write(6,'(5X,A)') 'PSEUDO-L EIGENFUNCTIONS:'
  write(6,*)
  if (mod(lDIMcf,2) == 1) then
    do I=1,lDIMcf
      write(6,'(A,I3,A,3X,20(2F9.6,1X))') '|',(lDIMcf-1)/2+(1-I),' > :',(Z(j,I),j=1,lDIMcf)
    end do
  else
    do I=1,lDIMcf
      write(6,'(A,I3,A,3X,20(2F9.6,1X))') '|',(lDIMcf-1)-2*(I-1),'/2 > :',(Z(j,I),j=1,lDIMcf)
    end do
  end if
end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  decompose the orbital moment AMSL in ITOs
!  transform AMFI integrals to pseudo-L basis:
! calculate the matrix elements of the spin and magnetic moment
! in the spin-orbit basis:
call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi_l,1)
do L=1,3
  call zcopy_(ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,tmp,1)
  ! amfi:
  call ZGEMM_('C','N',lDIMcf,lDIMcf,lDIMcf,(1.0_wp,0.0_wp),Z,lDIMcf,AMFI2(L,:,:),lDIMcf,(0.0_wp,0.0_wp),TMP,lDIMcf)
  call ZGEMM_('N','N',lDIMcf,lDIMcf,lDIMcf,(1.0_wp,0.0_wp),TMP,lDIMcf,Z,lDIMcf,(0.0_wp,0.0_wp),AMFI_L(L,:,:),lDIMcf)

end do !L

if (debug) call prMom_herm('TERMCF:: AMFI_L',amfi_l*au2cm,ldimcf)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the RASSI Crystal Field matrix
call dcopy_(lDIMcf,[0.0_wp],0,eloc,1)
call rtrace(lDIMcf,ESFS,ELOC)
call zcopy_(ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,HCF,1)
do i=1,lDIMcf
  do I1=1,lDIMcf
    do I2=1,lDIMcf
      HCF(I1,I2) = HCF(I1,I2)+ELOC(i)*conjg(Z(i,I1))*Z(i,I2)
    end do
  end do
end do

info = 0
call dcopy_(lDIMcf,[0.0_wp],0,Winit,1)
call zcopy_(lDIMcf*lDIMcf,[(0.0_wp,0.0_wp)],0,Zinit,1)
call DIAG_C2(HCF,lDIMcf,info,Winit,Zinit)
call print_ZFS('Ab Initio Calculated Crystal-Field Splitting Matrix written in the basis of Pseudo-L Eigenfunctions',HCF,lDIMCF)

call NEWCF(HCF,lDIMcf,Akq,BNC,BNS,BStev)
!ccccccccccccc print CF parameter ccccccccccccccccccccccccccccccc
if ((iprint >= 4) .or. debug) then
  call print_CFP_LCLU(lDIMCF,BNC,BNS,.true.)
  call print_CFP_stev(lDIMCF,Bstev,.true.)
  call print_CFP_naoya(lDIMcf,Akq,.true.)
else
  call print_CFP_LCLU(lDIMCF,BNC,BNS,.false.)
  call print_CFP_stev(lDIMCF,Bstev,.false.)
  call print_CFP_naoya(lDIMcf,Akq,.false.)
end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
write(6,'(/)')
if (mod(lDIMCF,2) == 1) then
  write(6,'(A,I0)') 'DECOMPOSITION OF THE RASSCF WAVE FUNCTIONS CORRESPONDING TO THE LOWEST ATOMIC MULTIPLET L =',(lDIMCF-1)/2
else
  write(6,'(A,I0,A)') 'DECOMPOSITION OF THE RASSCF WAVE FUNCTIONS CORRESPONDING TO THE LOWEST ATOMIC MULTIPLET L = ',(lDIMCF-1),'/2'
end if
write(6,'(A,I0)') 'IN WAVE FUNCTIONS WITH DEFINITE PROJECTION OF THE TOTAL MOMENT ON THE QUANTIZATION AXIS'
call print_ZFS_naoya('L',Zinit,lDIMcf)
call individual_ranks(lDIMCF,BNC,BNS,HCF,'L',iprint)

! saving some information for tests:
call Add_Info('CRYS_TERM_BNMC_20',BNC(2,0),1,4)
call Add_Info('CRYS_TERM_BNMC_40',BNC(4,0),1,4)
call Add_Info('CRYS_TERM_BNMC_60',BNC(6,0),1,4)
goto 999
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! generate |J,MJ> states using the spin |S,MS> and |L,ML> states

tS = 0.d0
tL = 0.d0
tJ = 0.d0
MS = 0
ML = 0
MJ = 0
if (nlanth == 1) then  ! Ce  J=5/2
  tS = 0.5d0
  tJ = 2.5d0
else if (nlanth == 2) then  ! Pr  J=4
  tS = 1.d0
  tJ = 4.d0
else if (nlanth == 3) then ! Nd  J=9/2
  tS = 1.5d0
  tJ = 4.5d0
else if (nlanth == 4) then ! Pm  J=4
  tS = 2.0d0
  tJ = 4.d0
else if (nlanth == 5) then ! Sm  J=5/2
  tS = 2.5d0
  tJ = 2.5d0
else if (nlanth == 6) then ! Eu  J=0
  tS = 3.0d0
  tJ = 0.0d0
else if (nlanth == 7) then ! Gd  J=7/2; S=7/2
  tS = 3.5d0
  tJ = 3.5d0
else if (nlanth == 8) then ! Tb  J=6
  tS = 3.0d0
  tJ = 6.0d0
else if (nlanth == 9) then ! Dy  J=15/2
  tS = 2.5d0
  tJ = 7.5d0
else if (nlanth == 10) then ! Ho  J=8
  tS = 2.0d0
  tJ = 8.0d0
else if (nlanth == 11) then ! Er  J=15/2
  tS = 1.5d0
  tJ = 7.5d0
else if (nlanth == 12) then ! Tm  J=6
  tS = 1.0d0
  tJ = 6.0d0
else if (nlanth == 13) then ! Yb  J=7/2
  tS = 0.5d0
  tJ = 3.5d0
else
  write(6,'(A)') 'not implemented yet'
end if
tL = (dble(lDIMcf)-1.d0)/2.d0
ML = lDIMcf
MS = nint(2.d0*tS+1.d0)
MJ = nint(2.d0*tJ+1.d0)

!write(6,*) 'nlanth=',nlanth,'  S,  L,  J =',tS,tL,tJ,' MS, ML, MJ =',MS,ML,MJ

!write(6,*) 'build coupled basis |L,ML>|S,MS>'
nLS = ML*MS
ij = 0
ir = 0
do i=-ML+1,ML-1,2
  ir = ir+1
  do j=-MS+1,MS-1,2
    ij = ij+1
    ibasL(ij) = i
    ibasS(ij) = j
    irootL(ij) = ir
  end do
end do
do iLS=1,nLS
  write(6,'(A,3I4,2x,A,2F6.1)') 'nLS,  ML,  MS =',iLS,ibasL(iLS),ibasS(iLS),'ML, MS =',dble(ibasL(iLS))/2.d0,dble(ibasS(iLS))/2.d0
end do

write(6,*) 'proceed to build |J,MJ> states'
ij = 0
do i=-MJ+1,MJ-1,2
  ij = ij+1
  ibasJ(ij) = i
end do

do i=1,MJ
  write(6,'(i3,i6,F6.1)') i,ibasJ(i),dble(ibasJ(i))/2.d0
end do

Cf = 0.0_wp
CfC = (0.0_wp,0.0_wp)

do ij=1,MJ
  tJM = dble(ibasJ(ij))/2.d0
  do iLS=1,nLS
    ! set projections
    spinM = dble(ibasS(iLS))/2.d0
    orbM = dble(ibasL(iLS))/2.d0

    call Clebsch_Gordan(tL,orbM,tS,spinM,tJ,tJM,coeffCG)
    Cf(iJ,iLS) = coeffCG

    if (abs(coeffCG) > 1.d-20) write(6,*) 'ij,iLS,coeffCG',ij,iLS,coeffCG
  end do
end do

write(6,*) 'MJ ->  (1,16), (2,15), (3,14), (4,13), (5,12), (6,11), (7,10), (8,9)'
do iLS=1,nLS
  write(6,'(A,2F6.1,16F11.8)') 'ML,MS: ',dble(ibasL(iLS))/2.d0,dble(ibasS(iLS))/2.d0,Cf(1,iLS),Cf(MJ,iLS),Cf(2,iLS),Cf(MJ-1,iLS), &
                               Cf(3,iLS),Cf(MJ-2,iLS),Cf(4,iLS),Cf(MJ-3,iLS),Cf(5,iLS),Cf(MJ-4,iLS),Cf(6,iLS),Cf(MJ-5,iLS), &
                               Cf(7,iLS),Cf(MJ-6,iLS),Cf(8,iLS),Cf(MJ-7,iLS)
end do

write(6,*) 're-write initial CASSCF states into |J,MJ>, using ( Z(j,I),j=1,lDIMcf)  coefficients'

do ij=1,MJ
  do iLS=1,nLS

    k = ibasL(iLS)/2

    do iCAS=1,ML
      CFC(iJ,iLS) = CFC(iJ,iLS)+Z(k,iCAS)*Cf(iJ,iLS)
    end do

  end do
end do

write(6,*) 'MJ ->  (1,16), (2,15), (3,14), (4,13), (5,12), (6,11), (7,10), (8,9)'
do iLS=1,nLS
  write(6,'(A,i2,F6.1,16(2F8.4,2x))') 'iCAS,MS: ',irootL(iLS),dble(ibasS(iLS))/2.d0,CfC(1,iLS),CfC(MJ,iLS),CfC(2,iLS), &
                                      CfC(MJ-1,iLS),CfC(3,iLS),CfC(MJ-2,iLS),CfC(4,iLS),CfC(MJ-3,iLS) !,CfC(5,iLS),CfC(MJ-4,iLS), &
                                      !CfC(6,iLS),CfC(MJ-5,iLS),CfC(7,iLS),CfC(MJ-6,iLS),CfC(8,iLS),CfC(MJ-7,iLS)
end do

999 continue
!-----------------------------------------------------------------------
call mma_deallocate(maxes)
call mma_deallocate(gtens)
call mma_deallocate(eloc)
call mma_deallocate(Winit)

call mma_deallocate(Angm)
call mma_deallocate(dipso)
call mma_deallocate(amfi_c)
call mma_deallocate(amfi2)
call mma_deallocate(amfi_l)
call mma_deallocate(Z)
call mma_deallocate(Zinit)
call mma_deallocate(tmp)
call mma_deallocate(HCF)

return

end subroutine termCF
