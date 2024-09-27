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

subroutine termCF(ANGMOM,AMFI,ESFS,ldimcf,d,maxes2,iopt,nlanth,iprint)
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
!                                 ( determinant(maxes)=1.0, orthogonal vectors )
!  iopt                        -- option for choosing the main quanization axis
!                   iopt = 1   -- axis of the ground orbital multiplet,
!                                 d specifies the size of pseudo L
!                   iopt = 2   -- axis of the entire L manifold
!                   iopt = 3   -- maxes is defined by the user
!                                 ( maxes2 is given as input )
!                   iopt = 4   -- maxes is the unity matrix ( original Z
!                                 is the quantization axis )

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne, Onei, auTocm
#ifdef _DISABLED_
use Constants, only: One, Two, Half
#endif
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ldimcf, d, iopt, nlanth, iprint
real(kind=wp), intent(in) :: angmom(3,ldimcf,ldimcf), amfi(3,ldimcf,ldimcf), esfs(ldimcf), maxes2(3,3)
integer(kind=iwp) :: i, i2, info, j, l
real(kind=wp) :: det, gtens(3), maxes(3,3)
real(kind=wp), allocatable :: BNC(:,:), BNS(:,:), Bstev(:,:), eloc(:), Winit(:)
complex(kind=wp), allocatable :: Akq(:,:), amfi2(:,:,:), amfi_c(:,:,:), amfi_l(:,:,:), amfi_tmp(:,:), Angm(:,:,:), dipso(:,:,:), &
                                 dipso_tmp(:,:,:), HCF(:,:), tmp(:,:), Z(:,:), Zinit(:,:)
logical(kind=iwp), parameter :: debug = .false.
real(kind=wp), external :: finddetr
#ifdef _DISABLED_
integer(kind=iwp) :: ibasJ(100), ibasL(100), ibasS(100), icas, ij, iLS, ir, irootL(100), k, MJ, ML, MS, nLS
real(kind=wp) :: CF(100,100), coeffCG, orbM, spinM, tJ, tJM, tL, tS
complex(kind=wp) :: CFC(100,100)
#else
#include "macros.fh"
unused_var(nlanth)
#endif

call mma_allocate(eloc,lDIMcf,label='eloc')
call mma_allocate(Winit,lDIMcf,label='Winit')

call mma_allocate(Angm,3,ldimcf,ldimcf,label='angm')
call mma_allocate(dipso,3,ldimcf,ldimcf,label='dipso')
call mma_allocate(amfi_c,3,ldimcf,ldimcf,label='amfi_c')
call mma_allocate(amfi2,3,ldimcf,ldimcf,label='amfi2')
call mma_allocate(amfi_l,3,ldimcf,ldimcf,label='amfi_l')
call mma_allocate(Z,ldimcf,ldimcf,label='Z')
call mma_allocate(Zinit,ldimcf,ldimcf,label='Zinit')
call mma_allocate(tmp,ldimcf,ldimcf,label='tmp')
call mma_allocate(HCF,ldimcf,ldimcf,label='HCF')
call mma_allocate(BNC,[1,lDIMcf],[0,lDIMcf],label='BNC')
call mma_allocate(BNS,[1,lDIMcf],[0,lDIMcf],label='BNS')
call mma_allocate(Bstev,[1,lDIMcf],[-lDIMcf,lDIMcf],label='Bstev')
call mma_allocate(Akq,[1,lDIMcf-1],[-(lDIMcf-1),lDIMcf-1],label='Akq')
!-----------------------------------------------------------------------
write(u6,'(/)')
write(u6,'(A)') repeat('%',95)
if (mod(lDIMcf,2) == 1) then
  write(u6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L = ',(lDIMcf-1)/2,'.'
else
  write(u6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC TERM, L = ',lDIMcf-1,'/2.'
end if
write(u6,'(A)') repeat('%',95)
write(u6,*)

dipso(:,:,:) = -angmom(:,:,:)*Onei
amfi_c(:,:,:) = amfi(:,:,:)*cOne
! find the main anisotropy direction of the
! diagonalize the angmom
if (iopt == 1) then
  ! iopt = 1   -- axis of the ground orbital doublet
  call mma_allocate(dipso_tmp,3,d,d,label='dipso_tmp')
  dipso_tmp(:,:,:) = dipso(:,1:d,1:d)
  call atens(dipso_tmp,d,gtens,maxes,iprint)
  call mma_deallocate(dipso_tmp)
  write(u6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  if (mod(d,2) == 0) then
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseudo-L = |',d-1,'/2> orbital multiplet.'
  else
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground pseudo-L = |',(d-1)/2,'> orbital multiplet.'
  end if

else if (iopt == 2) then
  ! iopt = 2   -- axis of the entire L manifold
  call atens(dipso,lDIMcf,gtens,maxes,iprint)
  if (mod(lDIMCF,2) == 0) then
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic L = |',lDIMCF-1,'/2> multiplet'
  else
    write(u6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic axes of the ground atomic L = |',(lDIMCF-1)/2,'> multiplet'
  end if

else if (iopt == 3) then
  ! iopt = 3   -- maxes is defined by the user ( maxes is given as input )
  write(u6,'(a)') 'The parameters of the Crystal Field matrix are written in the coordinate system:'
  write(u6,'(a)') '(Xm, Ym, Zm) -- defined in the input file.'
  ! copy the maxes2 to maxes
  maxes(:,:) = maxes2(:,:)
  det = finddetr(maxes,3)
  if (det == Zero) then
    write(u6,'(A)') 'TermCF:   iopt=3, while  DET(maxes)= 0.0'
    call AbEnd()
  end if
  ! FIXME: note finddetr modifies the matrix!
  if (det < Zero) then
    maxes(:,1) = -maxes2(:,1)
    if (iprint > 2) write(u6,'(a)') 'The original coordinate system was LEFT-handed. It has been changed to the RIGHT-handed'
  end if

else if (iopt == 4) then
  ! iopt = 4   -- maxes is the unity matrix ( original Z is the quantization axis )
  call unitmat(maxes,3)
else
  call AbEnd()
end if

! print out the rotation matrix:
write(u6,'(a)') 'Rotation matrix from the initial coordinate system to the employed coordinate system is:'

if ((iopt == 1) .or. (iopt == 2)) then
  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(u6,'(A,35x,A)') 'Xm, Ym, Zm -- main magnetic axes','|'
  write(u6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Xm |',(maxes(j,1),j=1,3),'|'
  write(u6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(maxes(j,2),j=1,3),'|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Zm |',(maxes(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,I3)') 'Quantization axis is Zm.'

else if (iopt == 3) then

  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
  write(u6,'(A,11x,A)') 'Xm, Ym, Zm -- the coordinate system defined in the input','|'
  write(u6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Xm |',(maxes(j,1),j=1,3),'|'
  write(u6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',(maxes(j,2),j=1,3),'|'
  write(u6,'(6x,A,3F18.14,1x,A)') '| Zm |',(maxes(j,3),j=1,3),'|'
  write(u6,'(2a)') repeat('-',67),'|'
  write(u6,'(A,I3)') 'Quantization axis is Zm.'

else

  write(u6,'(A)') 'IDENTITY matrix.'
  write(u6,'(A,I3)') 'Quantization axis is the initial z axis.'
end if

! rotate the angular momentum to the new axes, using "maxes"
call rotmom2(dipso,ldimcf,maxes,angm)
call rotmom2(amfi_c,ldimcf,maxes,amfi2)
if (debug) call prmom('TERMCF:: ANGM',angm,ldimcf)
if (debug) call prmom('TERMCF:: AMFI2',amfi2,ldimcf)

!-----------------------------------------------------------------------
! Find the pseudo-L basis of the LS manifold
call pseudospin(angm,ldimcf,Z,3,1,iprint)

if ((iprint >= 4) .or. debug) then
  write(u6,*)
  write(u6,'(5X,A)') 'PSEUDO-L EIGENFUNCTIONS:'
  write(u6,*)
  if (mod(lDIMcf,2) == 1) then
    do I=1,lDIMcf
      write(u6,'(A,I3,A,3X,20(2F9.6,1X))') '|',(lDIMcf-1)/2+(1-I),' > :',(Z(j,I),j=1,lDIMcf)
    end do
  else
    do I=1,lDIMcf
      write(u6,'(A,I3,A,3X,20(2F9.6,1X))') '|',(lDIMcf-1)-2*(I-1),'/2 > :',(Z(j,I),j=1,lDIMcf)
    end do
  end if
end if
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  decompose the orbital moment AMSL in ITOs
!  transform AMFI integrals to pseudo-L basis:
! calculate the matrix elements of the spin and magnetic moment
! in the spin-orbit basis:
call mma_allocate(amfi_tmp,ldimcf,ldimcf,label='amfi_tmp')
do L=1,3
  ! amfi:
  amfi_tmp(:,:) = AMFI2(L,:,:)
  call ZGEMM_('C','N',lDIMcf,lDIMcf,lDIMcf,cOne,Z,lDIMcf,amfi_tmp,lDIMcf,cZero,TMP,lDIMcf)
  call ZGEMM_('N','N',lDIMcf,lDIMcf,lDIMcf,cOne,TMP,lDIMcf,Z,lDIMcf,cZero,amfi_tmp,lDIMcf)
  AMFI_L(L,:,:) = amfi_tmp(:,:)
end do !L
call mma_deallocate(amfi_tmp)

if (debug) call prMom_herm('TERMCF:: AMFI_L',amfi_l*auTocm,ldimcf)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate the RASSI Crystal Field matrix
call rtrace(lDIMcf,ESFS,ELOC)
HCF(:,:) = cZero
do i=1,lDIMcf
  do I2=1,lDIMcf
    HCF(:,I2) = HCF(:,I2)+ELOC(i)*conjg(Z(i,:))*Z(i,I2)
  end do
end do

info = 0
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
write(u6,'(/)')
if (mod(lDIMCF,2) == 1) then
  write(u6,'(A,I0)') 'DECOMPOSITION OF THE RASSCF WAVE FUNCTIONS CORRESPONDING TO THE LOWEST ATOMIC MULTIPLET L =',(lDIMCF-1)/2
else
  write(u6,'(A,I0,A)') 'DECOMPOSITION OF THE RASSCF WAVE FUNCTIONS CORRESPONDING TO THE LOWEST ATOMIC MULTIPLET L = ',lDIMCF-1,'/2'
end if
write(u6,'(A,I0)') 'IN WAVE FUNCTIONS WITH DEFINITE PROJECTION OF THE TOTAL MOMENT ON THE QUANTIZATION AXIS'
call print_ZFS_naoya('L',Zinit,lDIMcf)
call individual_ranks(lDIMCF,BNC,BNS,HCF,'L',iprint)

! saving some information for tests:
call Add_Info('CRYS_TERM_BNMC_20',BNC(2,0),1,4)
call Add_Info('CRYS_TERM_BNMC_40',BNC(4,0),1,4)
call Add_Info('CRYS_TERM_BNMC_60',BNC(6,0),1,4)
#ifdef _DISABLED_
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! generate |J,MJ> states using the spin |S,MS> and |L,ML> states

select case (nlanth)
  case (1)  ! Ce  J=5/2
    tS = 0.5_wp
    tJ = 2.5_wp
  case (2)  ! Pr  J=4
    tS = 1.0_wp
    tJ = 4.0_wp
  case (3)  ! Nd  J=9/2
    tS = 1.5_wp
    tJ = 4.5_wp
  case (4)  ! Pm  J=4
    tS = 2.0_wp
    tJ = 4.0_wp
  case (5)  ! Sm  J=5/2
    tS = 2.5_wp
    tJ = 2.5_wp
  case (6)  ! Eu  J=0
    tS = 3.0_wp
    tJ = 0.0_wp
  case (7)  ! Gd  J=7/2; S=7/2
    tS = 3.5_wp
    tJ = 3.5_wp
  case (8)  ! Tb  J=6
    tS = 3.0_wp
    tJ = 6.0_wp
  case (9)  ! Dy  J=15/2
    tS = 2.5_wp
    tJ = 7.5_wp
  case (10) ! Ho  J=8
    tS = 2.0_wp
    tJ = 8.0_wp
  case (11) ! Er  J=15/2
    tS = 1.5_wp
    tJ = 7.5_wp
  case (12) ! Tm  J=6
    tS = 1.0_wp
    tJ = 6.0_wp
  case (13) ! Yb  J=7/2
    tS = 0.5_wp
    tJ = 3.5_wp
  case default
    tS = Zero
    tJ = Zero
    write(u6,'(A)') 'not implemented yet'
end select
tL = real(lDIMcf-1,kind=wp)*Half
ML = lDIMcf
MS = nint(Two*tS+One)
MJ = nint(Two*tJ+One)

!write(u6,*) 'nlanth=',nlanth,'  S,  L,  J =',tS,tL,tJ,' MS, ML, MJ =',MS,ML,MJ

!write(u6,*) 'build coupled basis |L,ML>|S,MS>'
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
  write(u6,'(A,3I4,2x,A,2F6.1)') 'nLS,  ML,  MS =',iLS,ibasL(iLS),ibasS(iLS),'ML, MS =',real(ibasL(iLS),kind=wp)*Half, &
                                 real(ibasS(iLS),kind=wp)*Half
end do

write(u6,*) 'proceed to build |J,MJ> states'
ij = 0
do i=-MJ+1,MJ-1,2
  ij = ij+1
  ibasJ(ij) = i
end do

do i=1,MJ
  write(u6,'(i3,i6,F6.1)') i,ibasJ(i),real(ibasJ(i),kind=wp)*Half
end do

Cf(:,:) = Zero
CfC(:,:) = cZero

do ij=1,MJ
  tJM = real(ibasJ(ij),kind=wp)*Half
  do iLS=1,nLS
    ! set projections
    spinM = real(ibasS(iLS),kind=wp)*Half
    orbM = real(ibasL(iLS),kind=wp)*Half

    call Clebsch_Gordan(tL,orbM,tS,spinM,tJ,tJM,coeffCG)
    Cf(iJ,iLS) = coeffCG

    if (abs(coeffCG) > 1.0e-20_wp) write(u6,*) 'ij,iLS,coeffCG',ij,iLS,coeffCG
  end do
end do

write(u6,*) 'MJ ->  (1,16), (2,15), (3,14), (4,13), (5,12), (6,11), (7,10), (8,9)'
do iLS=1,nLS
  write(u6,'(A,2F6.1,16F11.8)') 'ML,MS: ',real(ibasL(iLS),kind=wp)*Half,real(ibasS(iLS),kind=wp)*Half,Cf(1,iLS),Cf(MJ,iLS), &
                                Cf(2,iLS),Cf(MJ-1,iLS),Cf(3,iLS),Cf(MJ-2,iLS),Cf(4,iLS),Cf(MJ-3,iLS),Cf(5,iLS),Cf(MJ-4,iLS), &
                                Cf(6,iLS),Cf(MJ-5,iLS),Cf(7,iLS),Cf(MJ-6,iLS),Cf(8,iLS),Cf(MJ-7,iLS)
end do

write(u6,*) 're-write initial CASSCF states into |J,MJ>, using ( Z(j,I),j=1,lDIMcf)  coefficients'

do ij=1,MJ
  do iLS=1,nLS

    k = ibasL(iLS)/2

    do iCAS=1,ML
      CFC(iJ,iLS) = CFC(iJ,iLS)+Z(k,iCAS)*Cf(iJ,iLS)
    end do

  end do
end do

write(u6,*) 'MJ ->  (1,16), (2,15), (3,14), (4,13), (5,12), (6,11), (7,10), (8,9)'
do iLS=1,nLS
  write(u6,'(A,i2,F6.1,16(2F8.4,2x))') 'iCAS,MS: ',irootL(iLS),real(ibasS(iLS),kind=wp)*Half,CfC(1,iLS),CfC(MJ,iLS),CfC(2,iLS), &
                                       CfC(MJ-1,iLS),CfC(3,iLS),CfC(MJ-2,iLS),CfC(4,iLS),CfC(MJ-3,iLS) !,CfC(5,iLS),CfC(MJ-4,iLS), &
                                       !CfC(6,iLS),CfC(MJ-5,iLS),CfC(7,iLS),CfC(MJ-6,iLS),CfC(8,iLS),CfC(MJ-7,iLS)
end do
#endif

!-----------------------------------------------------------------------
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
call mma_deallocate(BNC)
call mma_deallocate(BNS)
call mma_deallocate(Bstev)
call mma_deallocate(Akq)

return

end subroutine termCF
