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

subroutine magnetization(nss,nM,nTempMagn,nDirTot,nDir,nDirZee,nH,iPrint,LUZee,mem,compute_Mdir_vector,zeeman_energy,hinput, &
                         m_paranoid,smagn,doplot,TempMagn,eso,dirX,dirY,dirZ,dir_weight,hexp,magn_exp,zJ,hmin,hmax,EM,thrs,dipm, &
                         sm,dbg)
!***********************************************************************
!                                                                      *
!     MAGNETIZATION control section                                    *
!                                                                      *
!     calling arguments:                                               *
!     NSS     : number of spin-orbit states (total)                    *
!               scalar integer                                         *
!     NM      : size of the Zeeman Hamiltonian matrix                  *
!               scalar integer                                         *
!     EM      : cut-off energy (energy of the last s-o state which is  *
!               included in the Zeeman matrix                          *
!               scalar real*8                                          *
!     EM      : cut-off energy (energy of the last s-o state which is  *
!               included in the Zeeman matrix                          *
!               scalar real*8                                          *
!     IFINAL  : integer                                                *
!               termination flag                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     Liviu Ungur                                                      *
!     University of Leuven, Belgium, 2008-2017                         *
!                                                                      *
!----------------------------------------------------------------------*
!     Hystory:                                                         *
!     Liviu Ungur, 2008-2017 various modifications                     *
!***********************************************************************

use Constants, only: Zero
use Definitions, only: wp, u6

implicit none
#include "mgrid.fh"
#include "stdalloc.fh"
!----------------------------------------------------------------
! input data (30)
integer, intent(in) :: nDir
integer, intent(in) :: nDirZee
integer, intent(in) :: nDirTot
integer, intent(in) :: nss
integer, intent(in) :: nM
integer, intent(in) :: nH
integer, intent(in) :: iprint
integer, intent(in) :: nTempMagn
integer, intent(in) :: mem
integer, intent(in) :: LUZee(nDirZee)
logical, intent(in) :: compute_Mdir_vector
logical, intent(in) :: zeeman_energy
logical, intent(in) :: DoPlot
logical, intent(in) :: hinput
logical, intent(in) :: smagn
logical, intent(in) :: m_paranoid
logical, intent(in) :: dbg
real(kind=8), intent(in) :: dirX(nDir), dirY(nDir), dirZ(nDir)
real(kind=8), intent(in) :: dir_weight(nDirZee,3)
real(kind=8), intent(in) :: hmin, hmax
real(kind=8), intent(in) :: zj, thrs
real(kind=8), intent(in) :: eso(nss)
real(kind=8), intent(in) :: EM
real(kind=8), intent(in) :: TempMagn(nTempMagn)
real(kind=8), intent(in) :: hexp(nH)
real(kind=8), intent(in) :: magn_exp(nH,nTempMagn)
complex(kind=8), intent(in) :: sm(3,nss,nss)
complex(kind=8), intent(in) :: dipm(3,nss,nss)
!-----------------------------------------------------------------------
! local variables
integer :: nP, iTEnd, iT, IM, I, L, J, IC, IDIR, IH, iTemp
real(kind=8) :: DLTH, mv, sv, dev
character(len=99) :: STLNE1, STLNE2
real(kind=8), allocatable :: WM(:)         ! WM(nm)
real(kind=8), allocatable :: MT(:,:,:)     ! MT(3,nH,nTempMagn)
real(kind=8), allocatable :: ST(:,:,:)     ! ST(3,nH,nTempMagn)
! magnetization and spin vectors
real(kind=8), allocatable :: MVEC(:,:,:,:) ! MVEC(nDirTot,nH,nTempMagn,3)
real(kind=8), allocatable :: SVEC(:,:,:,:) ! SVEC(nDirTot,nH,nTempMagn,3)
real(kind=8), allocatable :: H(:)          ! H(nH)
! average powder M and S:
real(kind=8), allocatable :: MAV(:,:) ! MAV(nH,nTempMagn)
real(kind=8), allocatable :: SAV(:,:) ! SAV(nH,nTempMagn)
real(kind=8), allocatable :: ZT(:,:)  ! ZT(nH,nTempMagn)
real(kind=8), allocatable :: STDEV(:) ! STDEV(nTempMagn)
real(kind=8), allocatable :: dHX(:)   ! dHX(nDirTot)
real(kind=8), allocatable :: dHY(:)   ! dHY(nDirTot)
real(kind=8), allocatable :: dHZ(:)   ! dHZ(nDirTot)
real(kind=8), allocatable :: dHW(:)   ! dHW(nDirTot)
integer :: mem_local, RtoB
external :: dev

!-----------------------------------------------------------------------
! Allocate necessary memory
mem_local = 0
RtoB = 8

! Zeeman exchange energy spectrum
call mma_allocate(WM,nM,'W')
call dcopy_(nM,[Zero],0,WM,1)
mem_local = mem_local+nM*RtoB

call mma_allocate(MT,3,nH,nTempMagn,'MT')
call dcopy_(3*nH*nTempMagn,[Zero],0,MT,1)
mem_local = mem_local+3*nH*nTempMagn*RtoB

call mma_allocate(ST,3,nH,nTempMagn,'ST')
call dcopy_(3*nH*nTempMagn,[Zero],0,ST,1)
mem_local = mem_local+3*nH*nTempMagn*RtoB

call mma_allocate(MAV,nH,nTempMagn,'MAV')
call dcopy_(nH*nTempMagn,[Zero],0,MAV,1)
mem_local = mem_local+nH*nTempMagn*RtoB

call mma_allocate(SAV,nH,nTempMagn,'SAV')
call dcopy_(nH*nTempMagn,[Zero],0,SAV,1)
mem_local = mem_local+nH*nTempMagn*RtoB

call mma_allocate(ZT,nH,nTempMagn,'ZT')
call dcopy_(nH*nTempMagn,[Zero],0,ZT,1)
mem_local = mem_local+nH*nTempMagn*RtoB

call mma_allocate(MVEC,nDirTot,nH,nTempMagn,3,'MVEC')
call mma_allocate(SVEC,nDirTot,nH,nTempMagn,3,'SVEC')
call dcopy_(3*nDirTot*nH*nTempMagn,[Zero],0,MVEC,1)
call dcopy_(3*nDirTot*nH*nTempMagn,[Zero],0,SVEC,1)
mem_local = mem_local+6*nDirTot*nH*nTempMagn*RtoB

call mma_allocate(H,nH,'H')
call dcopy_(nH,[Zero],0,H,1)
mem_local = mem_local+nH*RtoB

call mma_allocate(STDEV,nTempMagn,'H')
call dcopy_(nTempMagn,[Zero],0,STDEV,1)
mem_local = mem_local+nTempMagn*RtoB

call mma_allocate(dHX,nDirTot,'dHX')
call mma_allocate(dHY,nDirTot,'dHY')
call mma_allocate(dHZ,nDirTot,'dHZ')
call mma_allocate(dHW,nDirTot,'dHW')
call dcopy_(nDirTot,[Zero],0,dHX,1)
call dcopy_(nDirTot,[Zero],0,dHY,1)
call dcopy_(nDirTot,[Zero],0,dHZ,1)
call dcopy_(nDirTot,[Zero],0,dHW,1)
mem_local = mem_local+4*nDirTot*RtoB
if (dbg) write(u6,*) 'MAGNETIZATION:  memory allocated (local):'
if (dbg) write(u6,*) 'mem_local=',mem_local
if (dbg) write(u6,*) 'MAGNETIZATION:  memory allocated (total):'
if (dbg) write(u6,*) 'mem_total=',mem+mem_local
!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(100A)') (('%'),J=1,96)
write(u6,'(40X,A)') 'CALCULATION OF THE MOLAR MAGNETIZATION'
write(u6,'(100A)') (('%'),J=1,96)
write(u6,*)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
if (DBG .or. (iprint > 3)) then
  write(u6,'(A, I5)') 'nH                  = ',nH
  write(u6,'(A, I5)') 'nM                  = ',nM
  write(u6,'(A, I5)') 'nsymm               = ',nsymm
  write(u6,'(A, I5)') 'ngrid               = ',ngrid
  write(u6,'(A, I5)') 'nDir                = ',nDir
  write(u6,'(A, I5)') 'nDirZee             = ',nDirZee
  write(u6,'(A, I5)') 'nDirTot             = ',nDirTot
  write(u6,'(A, F9.5)') 'HMIN                = ',hmin
  write(u6,'(A, F9.5)') 'HMAX                = ',hmax
  write(u6,'(A, F9.5)') 'zJ                  = ',zJ
  write(u6,*) 'compute_Mdir_vector = ',compute_Mdir_vector
  write(u6,*) 'hinput              = ',hinput
  write(u6,*) 'zeeman_energy       = ',zeeman_energy
  write(u6,*) 'hinput              = ',hinput
  write(u6,*) 'smagn               = ',smagn
  write(u6,'(A)') 'dir_weight'
  do i=1,nDirZee
    write(u6,'(3F10.6)') (dir_weight(i,j),j=1,3)
  end do
  write(u6,'(A)') 'nDir'
  do i=1,nDir
    write(u6,'(3F10.6)') dirX(i),dirY(i),dirZ(i)
  end do
  write(u6,'(30(F6.3,a))') (TempMagn(iTemp),' K.;',iTemp=1,nTempMagn)
  if (zeeman_energy) then
    write(u6,'(A)') 'dir_weight'
    do i=1,nDirZee
      write(u6,'(3F10.6)') (dir_weight(i,j),j=1,3)
    end do
  end if
end if

nP = get_nP(nsymm,ngrid)

call hdir(nDir,nDirZee,dirX,dirY,dirZ,dir_weight,nP,nsymm,ngrid,nDirTot,dHX,dHY,dHZ,dHW)

write(u6,'(2X,A,i3,A)') 'Molar magnetization will be calculated in ',NH,' points, equally distributed in magnetic field range'
write(u6,'(2X,F4.1,1x,a,1X,F4.1,a,10(F6.3,a))') HMIN,'--',HMAX,' T., at the following temperatures:'
do i=1,nTempMagn,10
  j = min(nTempMagn,i+9)
  write(u6,'(10(F8.4,A))') (TempMagn(l),' K.;',l=i,j)
end do
write(u6,'(2X,A,I4,A)') 'Powder molar magnetization will be averaged on ',nP,' directions of the applied magnetic field.'
write(u6,'(2x,10A)') ('--------',i=1,10)
if (nsymm == 1) then
  write(u6,'(23x,A)') 'Lebedev-Laikov grid on a hemisphere:'
  write(u6,'(38x,A)') 'z >= 0;'
else if (nsymm == 2) then
  write(u6,'(23x,A)') 'Lebedev-Laikov grid on a 4th-of-a-sphere:'
  write(u6,'(34x,A)') 'x >= 0; z >= 0;'
else if (nsymm == 3) then
  write(u6,'(23x,A)') 'Lebedev-Laikov grid on a 8th-of-a-sphere:'
  write(u6,'(30x,A)') 'x >= 0; y >= 0; z >= 0;'
end if
write(u6,'(2x,10A)') ('--------',i=1,10)
write(u6,'(2x,A,12x,A,2(18x,A),16x,A)') 'Nr.','x','y','z','weight'
do i=1,nP
  write(u6,'(i4,2x,4(F18.12,1x))') i,dHX(i+nDir+nDirZee),dHY(i+nDir+nDirZee),dHZ(i+nDir+nDirZee),dHW(i+nDir+nDirZee)
end do

if (nDir > 0) then
  write(u6,'(2x,10A)') ('--------',i=1,10)
  write(u6,'(23x,A)') ' Magnetization vector will be computed'
  write(u6,'(2x,10A)') ('--------',i=1,10)
  write(u6,'(2x,A,12x,A,2(18x,A))') 'Nr.','x','y','z'
  do i=1,nDir
    write(u6,'(i4,2x,4(F18.12,1x))') i,dHX(i),dHY(i),dHZ(i)
  end do
end if

if (zeeman_energy) then
  write(u6,'(2x,10A)') ('--------',i=1,10)
  write(u6,'(23x,A)') 'Zeeman Energy Splitting will be computed'
  write(u6,'(2x,10A)') ('--------',i=1,10)
  write(u6,'(2x,A,12x,A,2(18x,A))') 'Nr.','x','y','z'
  do i=1,nDirZee
    j = i+nDir
    write(u6,'(i4,2x,4(F18.12,1x))') j,dHX(j),dHY(j),dHZ(j)
  end do
end if

write(u6,'(2x,10A)') ('--------',i=1,10)
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
if (compute_Mdir_vector) then
  write(u6,'(2x,A,i2,a)') 'The magnetization vector for ',nDir,' directions of the applied magnetic field will be calculated.'
else
  write(u6,'(2X,A)') 'The magnetization vector was not calculated.'
end if
if (zeeman_energy) then
  write(u6,'(2x,A,i2,a)') 'The Zeeman splitting for ',nDirZee,' directions of the applied magnetic field will be calculated.'
  write(u6,'(2x,a     )') 'The Zeeman energies for each direction of the applied magnetic field are written in files '// &
                         '"zeeman_energy_xxx.txt".'
else
  write(u6,'(2X,A)') 'Computation of the Zeeman splitting was not requested.'
end if
!smagn = .false.
!m_paranoid = .true.
!THRS = 1.0e-10_wp ! threshold for convergence of average spin, in case (zJ /= 0)
! /// opening the loop over the field points
do iH=1,nH
  !/// -----------------------------------------------------------------
  if (HINPUT) then
    H(iH) = HEXP(iH)
    if (H(iH) == Zero) H(iH) = 0.0001_wp
  else
    DLTH = (HMAX-HMIN)/real(NH-1,kind=wp)
    if (iH == 1) then
      H(IH) = HMIN+0.0001_wp
    else
      H(IH) = HMIN+DLTH*real(IH-1,kind=wp)
    end if
    if (H(iH) == Zero) H(iH) = 0.0001_wp
  end if

  if (DBG) write(u6,'(A,i0,A,F10.5,A,L2,A,L2)') 'MAGNETIZATION::  H(',iH,') = ',H(iH),'smagn=',smagn,' m_paranoid=',m_paranoid

  !/// opening the loop over dIfferent directions of the magnetic field
  do iM=1,NDIRTOT
    ! Entry into monitor: Status line
    write(STLNE1,'(A)') 'SINGLE_ANISO:  powder magnetization:'
    write(STLNE2,'(A,I4,A,I4,A,I4,A,I4)') ' Field: ',IH,' from ',nH,' at direction ',IM,' from ',NDIRTOT
    call StatusLine(trim(STLNE1),trim(STLNE2))
    ! actual calculation of the MT and ST, ZT
    call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),H(iH),ESO,zJ,THRS,DIPM,SM,nTempMagn,TempMagn,smagn,WM,ZT(iH,:),ST(:,iH,:), &
              MT(:,iH,:),m_paranoid,DBG)
    if (DBG .and. (iH == nH) .and. (iM == 23)) then
      write(u6,'(A,3ES16.8)') 'iM:',dHX(iM),dHY(iM),dHZ(iM)
      write(u6,'(2(A,3ES16.8,1x),A,ES16.8)') 'MT:',(MT(l,iH,1),l=1,3),'ST:',(ST(l,iH,1),l=1,3),'ZSTAT:',ZT(iH,1)
      write(u6,'(A,3ES16.8)') 'WM:',(WM(l),l=1,nM)
    end if
    !-------------------------------------------------------------------
    if (zeeman_energy) then
      if ((iH == 1) .and. (iM == nDir+1)) write(u6,'(A)') 'Energies of the Zeeman Hamiltonian for the following directions of '// &
                                                         'the applied field:'
      if ((iH == 1) .and. (iM > nDir) .and. (iM <= nDir+nDirZee)) then
        write(u6,'(A,I3,A,3F10.6,3x,5A)') 'direction Nr.',iM-nDir,' : ',dHX(iM),dHY(iM),dHZ(iM),'written in file "zeeman_energy_', &
                                          char(48+mod(int((iM-nDir)/100),10)),char(48+mod(int((iM-nDir)/10),10)), &
                                          char(48+mod((iM-nDir),10)),'.txt".'

        write(LUZee(iM-nDir),'(A,3F24.15)') '# direction of the applied magnetic field:',dHX(iM),dHY(iM),dHZ(iM)
        write(LUZee(iM-nDir),'(A,6x,A,1000(I4,6x) )') '# H(T)',' State =>',(i,i=1,nm)
      end if

      if ((iM > nDir) .and. (iM <= nDir+nDirZee)) write(LUZee(iM-nDir),'(F8.4,1000F10.3)') H(IH),(WM(I),I=1,NM)
    end if !zeeman_energy
    !-------------------------------------------------------------------
    ! computing the AVERAGE MOMENTS calculated at different temperatures
    ! (TempMagn(i))
    do iT=1,nTempMagn
      do ic=1,3
        MVEC(iM,iH,iT,ic) = MT(ic,iH,iT)
        SVEC(iM,iH,iT,ic) = ST(ic,iH,iT)
      end do  ! ic

      if (iM > nDir+nDirZee) then
        ! MAV(iH,iTemp) = MAV(iH,iTemp)+MT(1,iH,iT)*dHX(iM)*dHW(iM)+MT(2,iH,iT)*dHY(iM)*dHW(iM)+MT(3,iH,iT)*dHZ(iM)*dHW(iM)
        ! SAV(iH,iTemp) = SAV(iH,iTemp)+ST(1,iH,iT)*dHX(iM)*dHW(iM)+ST(2,iH,iT)*dHY(iM)*dHW(iM)+ST(3,iH,iT)*dHZ(iM)*dHW(iM)
        ! accumulate contributions:
        call daxpy_(1,dHX(iM)*dHW(iM),MT(1,iH,iT),1,MAV(iH,iT),1)
        call daxpy_(1,dHY(iM)*dHW(iM),MT(2,iH,iT),1,MAV(iH,iT),1)
        call daxpy_(1,dHZ(iM)*dHW(iM),MT(3,iH,iT),1,MAV(iH,iT),1)
        call daxpy_(1,dHX(iM)*dHW(iM),ST(1,iH,iT),1,SAV(iH,iT),1)
        call daxpy_(1,dHY(iM)*dHW(iM),ST(2,iH,iT),1,SAV(iH,iT),1)
        call daxpy_(1,dHZ(iM)*dHW(iM),ST(3,iH,iT),1,SAV(iH,iT),1)
      end if
      if (iprint > 2) then
        if ((iM == 1) .and. (iH == 1)) write(u6,'(2x,A,1x,A,4x,A,7x,A,7x,A)') 'iH','iM','iT','moment(iM,iT)','spin(iM,iT)'
        write(u6,'(3i4,3(F21.15,1x))') iH,iM,iT,(MT(1,iH,iT)*dHX(iM)+MT(2,iH,iT)*dHY(iM)+MT(3,iH,iT)*dHZ(iM)), &
                                       (ST(1,iH,iT)*dHX(iM)+ST(2,iH,iT)*dHY(iM)+ST(3,iH,iT)*dHZ(iM))
      end if
    end do ! iT
    !-------------------------------------------------------------------
    !/// closing the loop over directions of magnetic field
  end do   ! iM
  !/// -----------------------------------------------------------------
  !/// closing the loop over the points of magnetic field:
end do ! IH
!/// -------------------------------------------------------------------
! Close Zeeman files, if opened
if (Zeeman_Energy) then
  do i=1,nDirZee
    close(LUZee(i))
  end do
end if

!/// -------------------------------------------------------------------
if (nDir > 0) then

  if (smagn) then
    do iT=1,nTempMagn
      write(u6,*)
      do iDir=1,nDir
        write(u6,*)
        write(u6,'(A,A,1x,A)') '--------|','------------------------------------------------------------|', &
                               '|------------------------------------------------------------|'
        write(u6,'(A,i3,26x,A,1x,A,60x,A)') 'Direction of the applied magnetic field:',iDir,'|','|','|'
        write(u6,'(A,F18.14,44x,A,1x,A,60x,A)') 'proj X=',dHX(iDIR),'|','|','|'
        write(u6,'(A,F18.14,44x,A,1x,A,60x,A)') 'proj Y=',dHY(iDir),'|','|','|'
        write(u6,'(A,F18.14,44x,A,1x,A,60x,A)') 'proj Z=',dHZ(iDir),'|','|','|'
        write(u6,'(A,F7.4,A,41x,A,1x,A,60x,A)') 'Temperature = ',TempMagn(iT),' kelvin','|','|','|'
        write(u6,'(A,A,1x,A)') '--------|','------------------------------------------------------------|', &
                               '|------------------------------------------------------------|'
        write(u6,'(2x,A,12x,2A,1x,A,10x,2A)') 'Field |','Magnetization Vector            |','   Total Magn. |','|', &
                                              'Spin Magnetization Vector         |','   Total Magn. |'
        write(u6,'(5A,1x,5A)') '--------|','--- proj X ---|','--- proj Y ---|','--- proj Z ---|','- in this dir.-|', &
                               '|--- proj X ---|','--- proj Y ---|','--- proj Z ---|','- in this dir.-|'
        do iH=1,nH
          mv = MVEC(iDir,iH,iT,1)*dHX(iDir)+MVEC(iDir,iH,iT,2)*dHY(iDir)+MVEC(iDir,iH,iT,3)*dHZ(iDir)
          sv = SVEC(iDir,iH,iT,1)*dHX(iDir)+SVEC(iDir,iH,iT,2)*dHY(iDir)+SVEC(iDir,iH,iT,3)*dHZ(iDir)

          write(u6,'(F7.3,1x,A, 3(ES13.6,1x,A),ES14.7,1x,A,1x,A,3(ES13.6,1x,A),ES14.7,1x,A)') H(iH),'|',MVEC(iDir,iH,iT,1),' ', &
                                                                                              MVEC(iDir,iH,iT,2),' ', &
                                                                                              MVEC(iDir,iH,iT,3),'|',mv,'|','|', &
                                                                                              SVEC(iDir,iH,iT,1),' ', &
                                                                                              SVEC(iDir,iH,iT,2),' ', &
                                                                                              SVEC(iDir,iH,iT,3),'|',sv,'|'
        end do
        write(u6,'(A,A,1x,A)') '--------|','------------------------------------------------------------|', &
                               '|------------------------------------------------------------|'
      end do ! iDir
    end do ! iT

  else ! smagn == .false.

    do iT=1,nTempMagn
      write(u6,*)
      do iDir=1,nDir
        write(u6,*)
        write(u6,'(2A)') '--------|','------------------------------------------------------------|'
        write(u6,'(A,i3,26x,A)') 'Direction of the applied magnetic field:',iDir,'|'
        write(u6,'(A,F18.14,44x,A)') 'proj X=',dHX(iDIR),'|'
        write(u6,'(A,F18.14,44x,A)') 'proj Y=',dHY(iDir),'|'
        write(u6,'(A,F18.14,44x,A)') 'proj Z=',dHZ(iDir),'|'
        write(u6,'(A,F7.4,A,41x,A)') 'Temperature = ',TempMagn(iT),' kelvin','|'
        write(u6,'(2A)') '--------|','------------------------------------------------------------|'
        write(u6,'(2x,A,12x,2A)') 'Field |','Magnetization Vector            |','   Total Magn. |'
        write(u6,'(5A)') '--------|','--- proj X ---|','--- proj Y ---|','--- proj Z ---|','- in this dir.-|'
        do iH=1,nH
          mv = MVEC(iDir,iH,iT,1)*dHX(iDir)+MVEC(iDir,iH,iT,2)*dHY(iDir)+MVEC(iDir,iH,iT,3)*dHZ(iDir)
          write(u6,'(F7.3,1x,A,3(ES13.6,1x,A),ES14.7,1x,A)') H(iH),'|',MVEC(iDir,iH,iT,1),' ',MVEC(iDir,iH,iT,2),' ', &
                                                             MVEC(iDir,iH,iT,3),'|',mv,'|'
        end do
        write(u6,'(2A)') '--------|','------------------------------------------------------------|'
      end do !iDir
    end do !iT
  end if !(smagn)
end if !(nDir>0)
!/// -------------------------------------------------------------------
! COMPUTING THE STANDARD DEVIATION OF THE MAGNETIZATION
if (HINPUT) then
  do iT=1,nTempMagn
    STDEV(iT) = dev(nH,MAV(:,iT),magn_exp(:,iT))
  end do
end if
!/// -------------------------------------------------------------------

write(u6,*)
write(u6,'(25X,A)') 'HIGH-FIELD POWDER MAGNETIZATION'
write(u6,'(30X,A)') '(Units: Bohr magneton)'
write(u6,*)
do iT=1,nTempMagn,5
  iTEnd = min(nTempMagn,iT+4)

  write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  write(u6,'(A,10(F10.3,A))') '   H(T)    |STATISTICAL SUM|',(TempMagn(i),' K.  |',i=iT,iTEnd)
  write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)

  do iH=1,nH
    write(u6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))') H(iH),'|',ZT(iH,1),'|',(MAV(iH,i),'|',i=iT,iTEnd)
  end do

  write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  if (HINPUT) then
    write(u6,'(A,15x,11(f14.10,1x,A) )') 'ST.DEV.M   |',(STDEV(i),'|',i=1,nTempMagn)
    write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  end if
end do

if (smagn) then
  write(u6,*)
  write(u6,'(15X,A)') 'HIGH-FIELD POWDER SPIN MAGNETIZATION'
  write(u6,'(20X,A)') '(Units: Bohr magneton)'
  write(u6,*)
  do iT=1,nTempMagn,5
    iTEnd = min(nTempMagn,iT+4)
    write(u6,'(A,11A)') '---------|',('---------------|',i=iT,iTEnd+1)
    write(u6,'(A,11(F10.3,A))') '   H(T)  |STATISTICAL SUM|',(TempMagn(i),' K.  |',i=iT,iTEnd)
    write(u6,'(A,11A)') '---------|',('---------------|',i=iT,iTEnd+1)
    do iH=1,nH
      write(u6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))') H(iH),'|',ZT(iH,1),'|',(SAV(iH,i),'|',i=iT,iTEnd)
    end do
    write(u6,'(A,11A)') '---------|',('---------------|',i=iT,iTEnd+1)
  end do
end if!smagn

if (DoPlot) then
  if (hinput) then
    call plot_MH_with_Exp(nH,H,nTempMagn,TempMagn,MAV,magn_exp)
  else
    call plot_MH_no_Exp(nH,H,nTempMagn,TempMagn,MAV)
  end if
  !if (zeeman_energy) then
  !  call plot_zeeman(nH,nM,nDirZee,H,LuZee)
  !end if
end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

call Add_Info('MAGN_AVERAGED',MAV(1:nH,1:nTempMagn),nH*nTempMagn,5)
if (compute_Mdir_vector) then
  call Add_Info('MAGN_VECT_X(2)     ',MVEC(1,2,1,1),1,4)
  call Add_Info('MAGN_VECT_X(nH/2)  ',MVEC(1,(NH-1)/2,1,1),1,4)
  call Add_Info('MAGN_VECT_X(nH)    ',MVEC(1,NH,1,1),1,4)
  call Add_Info('MAGN_VECT_Y(2)     ',MVEC(1,2,1,2),1,4)
  call Add_Info('MAGN_VECT_Y(nH/2)  ',MVEC(1,(NH-1)/2,1,2),1,4)
  call Add_Info('MAGN_VECT_Y(nH)    ',MVEC(1,NH,1,2),1,4)
  call Add_Info('MAGN_VECT_Z(2)     ',MVEC(1,2,1,3),1,4)
  call Add_Info('MAGN_VECT_Z(nH/2)  ',MVEC(1,(NH-1)/2,1,3),1,4)
  call Add_Info('MAGN_VECT_Z(nH)    ',MVEC(1,NH,1,3),1,4)
end if

!-----------------------------------------------------------------------
! Deallocate necessary memory
call mma_deallocate(WM)
call mma_deallocate(MT)
call mma_deallocate(ST)
call mma_deallocate(MAV)
call mma_deallocate(SAV)
call mma_deallocate(ZT)
call mma_deallocate(MVEC)
call mma_deallocate(SVEC)
call mma_deallocate(H)
call mma_deallocate(STDEV)
call mma_deallocate(dHX)
call mma_deallocate(dHY)
call mma_deallocate(dHZ)
call mma_deallocate(dHW)

return

end subroutine magnetization
