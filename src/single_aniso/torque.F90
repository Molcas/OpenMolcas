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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, deg2rad
use Definitions, only: wp, iwp, u6, CtoB, RtoB

implicit none
integer(kind=iwp), intent(in) :: nss, nm, AngPoints, mem
real(kind=wp), intent(in) :: EM, ESO(nss), zJ, thrs, H_torq, T_torq, ma(3,3) ! main magnetic axes
complex(kind=wp), intent(in) :: DIPM(3,nss,nss), SM(3,nss,nss)
logical(kind=iwp), intent(in) :: m_paranoid, smagn, DBG
integer(kind=iwp) :: I, IM, J, mem_local, nT_torq
real(kind=wp) :: AngRad, AngStep, g(3), mg(3,3), MT(3), ST(3), ZT(1) !, det, dlth, ma_inv(3,3)
complex(kind=wp) :: MM(3,2,2)
character(len=99) :: STLNE1, STLNE2
real(kind=wp), allocatable :: Ang(:), dX(:), dY(:), dZ(:), ty(:), W(:) !, tx(:,:), tz(:,:)
complex(kind=wp), allocatable :: M(:,:,:), S(:,:,:)
integer(kind=iwp), parameter :: nPlanes = 1

write(u6,*)
write(u6,'(A)') repeat('%',96)
write(u6,'(20X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(u6,'(A)') repeat('%',96)
write(u6,*)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!write(u6,'(2X,A,i3,A)') 'Magnetization torque is calculated for the ',NH,' field points, in the field domain:'
!write(u6,'(2X,F4.1,1x,a,1X,F4.1,a,10(F6.3,a))') HMIN,'--',HMAX,' T., at the following temperatures:', &
!                                                (TempMagn(l),' K.;',l=1,nTempMagn)
!do i=11,nTempMagn,10
!  j = MIN(nTempMagn,i+9)
!  write(u6,'(17x,10(F9.3,A))') (TempMagn(k),' K.;',k=i,j)
!end do
write(u6,'(2X,A,F10.5,A,F10.5)') 'Magnetization torque is calculated for one field point,',H_torq, &
                                 ' T., at the following temperature:',T_torq

write(u6,'(2x,A,i3,A)') 'Angular dependence of the magnetization torque is computed for ',AngPoints,' angular points distributed'
write(u6,'(2x,A)') 'in the domain 0-180 deg.'
write(u6,'(2X,A)') 'The cut-off energy for the exact diagonalization of the Zeeman Hamiltonian is:'
write(u6,'(2x,a,F15.9,A)') 'E = ',EM,' cm(-1).'
if (NM < 10) then
  write(u6,'(2X,A,i2,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 10) .and. (NM < 100)) then
  write(u6,'(2X,A,i3,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 100) .and. (NM < 1000)) then
  write(u6,'(2X,A,i4,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 1000) .and. (NM < 10000)) then
  write(u6,'(2X,A,i5,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
end if
nT_torq = 1

!if (dbg) then
do i=1,3
  write(u6,*) (ma(j,i),j=1,3)
end do
!end if
!-----------------------------------------------------------------------
! Allocate memory for this computation:
mem_local = 0
! Zeeman exchange energy spectrum
call mma_allocate(W,nM,'W')
mem_local = mem_local+size(W)*RtoB

call mma_allocate(dX,AngPoints,'dX')
call mma_allocate(dY,AngPoints,'dY')
call mma_allocate(dZ,AngPoints,'dZ')
mem_local = mem_local+size(dX)*RtoB
mem_local = mem_local+size(dY)*RtoB
mem_local = mem_local+size(dZ)*RtoB

call mma_allocate(Ang,AngPoints,'Ang')
mem_local = mem_local+size(Ang)*RtoB

call mma_allocate(ty,AngPoints,'ty')
mem_local = mem_local+size(ty)*RtoB

call mma_allocate(M,3,nss,nss,'Mrot')
call mma_allocate(S,3,nss,nss,'Srot')
mem_local = mem_local+size(M)*CtoB
mem_local = mem_local+size(S)*CtoB

if (dbg) then
  write(u6,*) 'TORQ:  memory allocated (local):'
  write(u6,*) 'mem_local=',mem_local
  write(u6,*) 'TORQ:  memory allocated (total):'
  write(u6,*) 'mem_total=',mem+mem_local
end if
!-----------------------------------------------------------------------
! rotate the moments to the coordiante system of
! the ground state
!call REVERSE(ma,ma_inv,DET)
call rotmom2(DIPM,nss,ma,M)
call rotmom2(SM,nss,ma,S)

MM(:,:,:) = M(:,1:2,1:2)
call atens(MM,2,g,mg,2)
!-----------------------------------------------------------------------
!call hdir2(AngPoints,2,dX,dY,dZ,Ang,2)
dY(:) = Zero
AngStep = 360.0_wp/real(AngPoints-1,kind=wp)
do i=1,AngPoints
  AngRad = real(i-1,kind=wp)*AngStep*deg2rad !+122.625_wp*deg2rad
  Ang(i) = real(i-1,kind=wp)*AngStep
  dX(i) = cos(AngRad)
  dZ(i) = sin(AngRad)
end do
if (dbg) then
  !write(u6,'(A,I5)') 'Angular grid for Magnetization Torque, Cartesian Component =',L
  write(u6,'(2x,A,4x,A,5x,3(10X,A,10x))') 'Nr.','Angle','X','Y','Z'
  do i=1,AngPoints
    write(u6,'(I4,F10.3,3x,3F21.14)') i,Ang(i),dX(i),dY(i),dZ(i)
  end do
end if
!-----------------------------------------------------------------------

do IM=1,AngPoints
  write(STLNE1,'(A   )') 'SINGLE_ANISO:  torque:'
  write(STLNE2,'(A,I3)') ' Magnetization at point ',IM
  call StatusLine(trim(STLNE1),trim(STLNE2))
  ! exchange magnetization:
  call MAGN(nss,NM,dX(iM),dY(iM),dZ(iM),H_torq,eso,zJ,THRS,M,S,nT_torq,[T_torq],smagn,W,ZT,ST,MT,m_paranoid,DBG)
  if (dbg) write(u6,'(A,3F18.10)') 'TORQ: MT=',MT(1),MT(2),MT(3)

  !tx(iPl,iM) = (MT(2)*dZ(iM)-MT(3)*dY(iM))*H_torq
  ty(iM) = (MT(3)*dX(iM)-MT(1)*dZ(iM))*H_torq
  !tz(iPl,iM) = (MT(1)*dY(iM)-MT(2)*dX(iM))*H_torq
end do ! iM
!-----------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(u6,'(30X,A)') '(Units of torque: [energy, cm-1])'
write(u6,*)

write(u6,'(5x,A)') 'Orientation of the applied magnetic field employed:'
write(u6,'(10A)') '--------|','---------------------------|'

write(u6,'(2x,A,10x,A)') 'Angle |','rotation in the XZ plane          |'
write(u6,'(10A)') '--------|','--- proj X ---|','--- proj Y ---|','--- proj Z ---|'
do iM=1,AngPoints
  write(u6,'(F7.3,1x,A,F13.10,1x,A,F13.10,1x,A,F13.10,1x,A,F20.14,1x,A)') Ang(iM),'|',dX(iM),' ',dY(iM),' ',dZ(iM),'|',ty(iM),'|'
end do
write(u6,'(10A)') '--------|','---------------------------|'

write(u6,*)
write(u6,'(10A)') '--------|','---------------------------|'
write(u6,'(A,F9.4,A)') 'Magnetic field strength = ',H_torq,' tesla'
write(u6,'(12x,A,F9.4,A)') 'Temperature = ',T_torq,' kelvin'
write(u6,'(10A)') '--------|','---------------------------|'

write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the XZ plane          |'
write(u6,'(10A)') '--------|','-- torque along Y --|'
do iM=1,AngPoints
  write(u6,'(F7.3,1x,A,ES18.10)') Ang(iM),'|',ty(iM)
end do
write(u6,'(10A)') '--------|','---------------------------|'

!-----------------------------------------------------------------------
! deallocate memory for this computation:
call mma_deallocate(W)
call mma_deallocate(dX)
call mma_deallocate(dY)
call mma_deallocate(dZ)
call mma_deallocate(Ang)
call mma_deallocate(ty)
call mma_deallocate(m)
call mma_deallocate(s)

if (dbg) write(u6,*) 'TORQ: allocated memory was sucessfully deallocated'

return

end subroutine torque
