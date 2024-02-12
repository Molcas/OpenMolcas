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

subroutine torque(Nss,NM,AngPoints,EM,eso,dipm,sm,zJ,thrs,mem,m_paranoid,smagn,H_torq,T_torq,ma,dbg)

implicit none
integer, parameter :: wp = kind(0.d0)
#include "cntrl_sa.fh"
#include "mgrid.fh"
#include "stdalloc.fh"
integer, intent(in) :: nss, nm
integer, intent(in) :: AngPoints
integer, intent(in) :: mem
logical, intent(in) :: smagn
logical, intent(in) :: m_paranoid
! ab initio data:
! exchange energies printed out in the previous part
real(kind=8), intent(in) :: ESO(nss)
real(kind=8), intent(in) :: EM, zJ, thrs
real(kind=8), intent(in) :: H_torq, T_torq
real(kind=8), intent(in) :: ma(3,3) ! main magnetic axes
complex(kind=8), intent(in) :: DIPM(3,nss,nss)
complex(kind=8), intent(in) :: SM(3,nss,nss)
! local data:
! magnetic field strength and orientation data:
integer :: nPlanes
parameter(nPlanes=1)
!real(kind=8) :: dlth
real(kind=8), allocatable :: W(:) ! W(NM) ! Zeeman exchange energies
real(kind=8) :: ZT(1)             ! ZT ! total statistical sum, Boltzmann distribution
real(kind=8), allocatable :: ST(:) ! ST(3) ! total spin magnetisation,
real(kind=8), allocatable :: MT(:) ! MT(3) ! total magnetisation
real(kind=8), allocatable :: dX(:) ! dX(AngPoints)
real(kind=8), allocatable :: dY(:) ! dY(AngPoints)
real(kind=8), allocatable :: dZ(:) ! dZ(AngPoints)
real(kind=8), allocatable :: Ang(:)  ! Ang(AngPoints)
complex(kind=8), allocatable :: M(:,:,:)
complex(kind=8), allocatable :: S(:,:,:)
! magnetic torque
!real(kind=8), allocatable :: tx(:,:) ! tx(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque
real(kind=8), allocatable :: ty(:) ! ty(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque
!real(kind=8), allocatable :: tz(:,:) ! tz(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque
! magnetic and spin moments (i.e. the BIG matrices):
character(len=99) :: STLNE1, STLNE2
real(kind=8) :: g(3), mg(3,3) !,ma_inv(3,3)!,det
real(kind=8) :: AngStep, AngRad, pi
logical :: DBG
integer :: IM, I, J, mem_local, RtoB, CtoB, nT_torq

!Boltz_k = 0.6950356000_wp ! in cm^-1*K-1
!mu_Bohr = 0.4668643740_wp ! in cm-1*T-1
!cm3tomB = 0.5584938904_wp ! in cm3 * mol-1 * T
pi = 3.1415926535897932384626433832795028841971_wp

write(6,*)
write(6,'(100A)') (('%'),J=1,96)
write(6,'(20X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(6,'(100A)') (('%'),J=1,96)
write(6,*)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!write(6,'(2X,A,i3,A)') 'Magnetization torque is calculated for the ',NH,' field points, in the field domain:'
!write(6,'(2X,F4.1,1x,a,1X,F4.1,a,10(F6.3,a))') HMIN,'--',HMAX,' T., at the following temperatures:', &
!                                               (TempMagn(l),' K.;',l=1,nTempMagn)
!do i=11,nTempMagn,10
!  j = MIN(nTempMagn,i+9)
!  write(6,'(17x,10(F9.3,A))') (TempMagn(k),' K.;',k=i,j)
!end do
write(6,'(2X,A,F10.5,A,F10.5)') 'Magnetization torque is calculated for one field point,',H_torq, &
                                ' T., at the following temperature:',T_torq

write(6,'(2x,A,i3,A)') 'Angular dependence of the magnetization torque is computed for ',AngPoints,' angular points distributed'
write(6,'(2x,A)') 'in the domain 0-180 deg.'
write(6,'(2X,A)') 'The cut-off energy for the exact diagonalization of the Zeeman Hamiltonian is:'
write(6,'(2x,a,F15.9,A)') 'E = ',EM,' cm(-1).'
if (NM < 10) then
  write(6,'(2X,A,i2,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 10) .and. (NM < 100)) then
  write(6,'(2X,A,i3,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 100) .and. (NM < 1000)) then
  write(6,'(2X,A,i4,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 1000) .and. (NM < 10000)) then
  write(6,'(2X,A,i5,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
end if
nT_torq = 1

!if (dbg) then
do i=1,3
  write(6,*) (ma(j,i),j=1,3)
end do
!end if
!-----------------------------------------------------------------------
! Allocate memory for this computation:
mem_local = 0
RtoB = 8
CtoB = 16
! Zeeman exchange energy spectrum
call mma_allocate(W,nM,'W')
call dcopy_(nM,[0.0_wp],0,W,1)
mem_local = mem_local+nM*RtoB

call mma_allocate(ST,3,'ST')
call dcopy_(3,[0.0_wp],0,ST,1)
mem_local = mem_local+3*RtoB

call mma_allocate(MT,3,'MT')
call dcopy_(3,[0.0_wp],0,MT,1)
mem_local = mem_local+3*RtoB

call mma_allocate(dX,AngPoints,'dX')
call mma_allocate(dY,AngPoints,'dY')
call mma_allocate(dZ,AngPoints,'dZ')
call dcopy_(AngPoints,[0.0_wp],0,dX,1)
call dcopy_(AngPoints,[0.0_wp],0,dY,1)
call dcopy_(AngPoints,[0.0_wp],0,dZ,1)
mem_local = mem_local+3*AngPoints*RtoB

call mma_allocate(Ang,AngPoints,'Ang')
call dcopy_(AngPoints,[0.0_wp],0,Ang,1)
mem_local = mem_local+AngPoints*RtoB

call mma_allocate(ty,AngPoints,'ty')
call dcopy_(AngPoints,[0.0_wp],0,ty,1)
mem_local = mem_local+AngPoints*RtoB

call mma_allocate(M,3,nss,nss,'Mrot')
call mma_allocate(S,3,nss,nss,'Srot')
call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,M,1)
call zcopy_(3*nss*nss,[(0.0_wp,0.0_wp)],0,S,1)
mem_local = mem_local+2*3*nss*nss*CtoB

if (dbg) write(6,*) 'TORQ:  memory allocated (local):'
if (dbg) write(6,*) 'mem_local=',mem_local
if (dbg) write(6,*) 'TORQ:  memory allocated (total):'
if (dbg) write(6,*) 'mem_total=',mem+mem_local
!-----------------------------------------------------------------------
! rotate the moments to the coordiante system of
! the ground state
!ma_inv = 0.0_wp
!call REVERSE(ma,ma_inv,DET)
call rotmom2(DIPM,nss,ma,M)
call rotmom2(SM,nss,ma,S)

g = 0.0_wp
mg = 0.0_wp
call atens(M(1:3,1:2,1:2),2,g,mg,2)
!-----------------------------------------------------------------------
call dcopy_(AngPoints,[0.0_wp],0,dX,1)
call dcopy_(AngPoints,[0.0_wp],0,dY,1)
call dcopy_(AngPoints,[0.0_wp],0,dZ,1)
call dcopy_(AngPoints,[0.0_wp],0,Ang,1)
!call hdir2(AngPoints,2,dX,dY,dZ,Ang,2)
AngStep = 0.0_wp
AngRad = 0.0_wp
AngStep = 360.0_wp/dble(AngPoints-1)
dX(1) = 1.0_wp
dZ(1) = 0.0_wp
do i=1,AngPoints
  AngRad = dble(i-1)*AngStep*Pi/180.0_wp !+122.625_wp*Pi/180.0_wp
  Ang(i) = dble(i-1)*AngStep
  dX(i) = cos(AngRad)
  dZ(i) = sin(AngRad)
end do
if (dbg) then
  !write(6,'(A,I5)') 'Angular grid for Magnetization Torque, Cartesian Component =',L
  write(6,'(2x,A,4x,A,5x,3(10X,A,10x))') 'Nr.','Angle','X','Y','Z'
  do i=1,AngPoints
    write(6,'(I4,F10.3,3x,3F21.14)') i,Ang(i),dX(i),dY(i),dZ(i)
  end do
end if
!-----------------------------------------------------------------------

do IM=1,AngPoints
  write(STLNE1,'(A   )') 'SINGLE_ANISO:  torque:'
  write(STLNE2,'(A,I3)') ' Magnetization at point ',IM
  call StatusLine(trim(STLNE1),trim(STLNE2))
  ZT = 0.0_wp
  call dcopy_(nM,[0.0_wp],0,W,1)
  call dcopy_(3,[0.0_wp],0,MT,1)
  call dcopy_(3,[0.0_wp],0,ST,1)
  ! exchange magnetization:
  call MAGN(nss,NM,dX(iM),dY(iM),dZ(iM),H_torq,eso,zJ,THRS,M,S,nT_torq,[T_torq],smagn,W,ZT,ST,MT,m_paranoid,DBG)
  if (dbg) write(6,'(A,3F18.10)') 'TORQ: MT=',MT(1),MT(2),MT(3)

  !tx(iPl,iM) = (MT(2)*dZ(iM)-MT(3)*dY(iM))*H_torq
  ty(iM) = (MT(3)*dX(iM)-MT(1)*dZ(iM))*H_torq
  !tz(iPl,iM) = (MT(1)*dY(iM)-MT(2)*dX(iM))*H_torq
end do ! iM
!-----------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
!-----------------------------------------------------------------------
write(6,*)
write(6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(6,'(30X,A)') '(Units of torque: [energy, cm-1])'
write(6,*)

write(6,'(5x,A)') 'Orientation of the applied magnetic field employed:'
write(6,'(10A)') '--------|','---------------------------|'

write(6,'(2x,A,10x,A)') 'Angle |','rotation in the XZ plane          |'
write(6,'(10A)') '--------|','--- proj X ---|','--- proj Y ---|','--- proj Z ---|'
do iM=1,AngPoints
  write(6,'(F7.3,1x,A,F13.10,1x,A,F13.10,1x,A,F13.10,1x,A,F20.14,1x,A)') Ang(iM),'|',dX(iM),' ',dY(iM),' ',dZ(iM),'|',ty(iM),'|'
end do
write(6,'(10A)') '--------|','---------------------------|'

write(6,*)
write(6,'(10A)') '--------|','---------------------------|'
write(6,'(A,F9.4,A)') 'Magnetic field strength = ',H_torq,' Tesla'
write(6,'(12x,A,F9.4,A)') 'Temperature = ',T_torq,' Kelvin'
write(6,'(10A)') '--------|','---------------------------|'

write(6,'(2x,A,3(10x,A))') 'Angle |','rotation in the XZ plane          |'
write(6,'(10A)') '--------|','-- torque along Y --|'
do iM=1,AngPoints
  write(6,'(F7.3,1x,A,ES18.10)') Ang(iM),'|',ty(iM)
end do
write(6,'(10A)') '--------|','---------------------------|'

!-----------------------------------------------------------------------
! deallocate memory for this computation:
call mma_deallocate(W)
call mma_deallocate(ST)
call mma_deallocate(MT)
call mma_deallocate(dX)
call mma_deallocate(dY)
call mma_deallocate(dZ)
call mma_deallocate(Ang)
call mma_deallocate(ty)
call mma_deallocate(m)
call mma_deallocate(s)

if (dbg) write(6,*) 'TORQ: allocated memory was sucessfully deallocated'

return

end subroutine torque
