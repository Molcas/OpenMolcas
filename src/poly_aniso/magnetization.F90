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

subroutine magnetization_pa(exch,nLoc,nM,nH,nneq,neq,neqv,nCenter,nTempMagn,nDir,nDirZee,nDirTot,nss,nexch,iopt,LUZee,TempMagn, &
                            hexp,mexp,hmin,hmax,em,zJ,thrs,dirX,dirY,dirZ,dir_weight,w,dipexch,s_exch,dipso,s_so,eso,hinput,r_rot, &
                            XLM,ZLM,XRM,ZRM,zeeman_energy,compute_Mdir_vector,m_paranoid,m_accurate,smagn,mem,doplot)

implicit none
integer, parameter :: wp = kind(0.d0)
#include "mgrid.fh"
#include "stdalloc.fh"
! constants defining the sizes
integer, intent(in) :: exch, nLoc, nM
integer, intent(in) :: nH, nCenter, nTempMagn
integer, intent(in) :: nDir, nDirZee, nDirTot, nneq, neqv
integer, intent(in) :: neq(nneq), nss(nneq), nexch(nneq)
integer, intent(in) :: iopt, mem
integer, intent(in) :: LUZee(nDirZee)
logical, intent(in) :: hinput
logical, intent(in) :: zeeman_energy
logical, intent(in) :: compute_Mdir_vector
logical, intent(in) :: m_paranoid
logical, intent(in) :: m_accurate
logical, intent(in) :: smagn
logical, intent(in) :: doplot
real(kind=8), intent(in) :: R_ROT(nneq,neqv,3,3)
real(kind=8), intent(in) :: W(exch)
! exchange energies printed out in the previous part
real(kind=8), intent(in) :: ESO(nneq,nLoc)
! spin-orbit energies from ANISO files
real(kind=8), intent(in) :: Hexp(nH), Mexp(nH,nTempMagn)
real(kind=8), intent(in) :: thrs
real(kind=8), intent(in) :: XLM(nCenter,nTempMagn,3,3)
real(kind=8), intent(in) :: ZLM(nCenter,nTempMagn)
real(kind=8), intent(in) :: XRM(nCenter,nTempMagn,3,3)
real(kind=8), intent(in) :: ZRM(nCenter,nTempMagn)
real(kind=8), intent(in) :: dirX(nDir), dirY(nDir), dirZ(nDir)
real(kind=8), intent(in) :: dir_weight(nDirZee,3)
real(kind=8), intent(in) :: TempMagn(nTempMagn)
real(kind=8), intent(in) :: zJ, hmin, hmax, em
complex(kind=8), intent(in) :: DIPEXCH(3,exch,exch)
complex(kind=8), intent(in) :: S_EXCH(3,exch,exch)
complex(kind=8), intent(in) :: dipso(nneq,3,nLoc,nLoc)
complex(kind=8), intent(in) :: s_so(nneq,3,nLoc,nLoc)
! exchange data:
!integer :: NM ! number of states included in the exchange Zeeman matrix, ( Nex <= exch)
real(kind=8), allocatable :: Wex(:) !WEX(NM) ! Zeeman exchange energies
real(kind=8), allocatable :: Zex(:) !ZEX(nTempMagn) ! exchange statistical sum, Boltzmann distribution
real(kind=8), allocatable :: Sex(:,:) !SEX(3,nTempMagn) ! spin magnetisation, from the exchange block;
real(kind=8), allocatable :: Mex(:,:) !MEX(3,nTempMagn) ! magnetisation, from the exchange block
! data for individual sites (all states):
real(kind=8), allocatable :: ZL(:,:) !ZL(nneq,nTempMagn) ! local statistical sum, Boltzmann distribution
real(kind=8), allocatable :: WL(:,:) !WL(nneq,nLoc) ! Zeeman local energis
real(kind=8), allocatable :: SL(:,:,:) !SL(nneq,3,nTempMagn) ! spin magnetisation, from the local sites, using ALL states ;
real(kind=8), allocatable :: ML(:,:,:) !ML(nneq,3,nTempMagn) ! magnetisation, from local sites, using ALL states;
! data for individual sites (only states that enter exchange):
real(kind=8), allocatable :: ZR(:,:) !ZR(nneq,nTempMagn) ! local statistical sum, Boltzmann distribution, using only Nexch states
real(kind=8), allocatable :: WR(:,:) !WR(nneq,nLoc) ! Zeeman local reduced energies, using only Nexch states;
real(kind=8), allocatable :: SR(:,:,:) !SR(nneq,3,nTempMagn) ! spin magnetisation, from the local sites, using only Nexch states ;
real(kind=8), allocatable :: MR(:,:,:) !MR(nneq,3,nTempMagn) ! magnetisation, from local sites, using only Nexch states;
! total vectors in general coordinate system:
real(kind=8), allocatable :: ZRT(:,:) !ZRT(nCenter,nTempMagn)
real(kind=8), allocatable :: ZLT(:,:) !ZLT(nCenter,nTempMagn)
real(kind=8), allocatable :: MRT(:,:,:) !MRT(nCenter,3,nTempMagn)
real(kind=8), allocatable :: MLT(:,:,:) !MLT(nCenter,3,nTempMagn)
real(kind=8), allocatable :: SRT(:,:,:) !SRT(nCenter,3,nTempMagn)
real(kind=8), allocatable :: SLT(:,:,:) !SLT(nCenter,3,nTempMagn)
! data for total system:
real(kind=8), allocatable :: ZT(:,:) !ZT(nH,nTempMagn) ! total statistical sum, Boltzmann distribution
real(kind=8), allocatable :: ST(:,:,:) !ST(3,nH,nTempMagn) ! total spin magnetisation,
real(kind=8), allocatable :: MT(:,:,:) !MT(3,nH,nTempMagn) ! total magnetisation
! magnetic field strength and orientation data:
real(kind=8) :: dltH
real(kind=8), allocatable :: H(:) !H(nH)
real(kind=8), allocatable :: dHX(:) !dHX(nDirTot)
real(kind=8), allocatable :: dHY(:) !dHY(nDirTot)
real(kind=8), allocatable :: dHZ(:) !dHZ(nDirTot)
real(kind=8), allocatable :: dHW(:) !dHW(nDirTot)
! total average M and average S data:
real(kind=8), allocatable :: MAV(:,:) !MAV(nH,nTempMagn)
real(kind=8), allocatable :: SAV(:,:) !SAV(nH,nTempMagn)
real(kind=8), allocatable :: MVEC(:,:,:,:) !MVEC(nDirTot,nH,nTempMagn,3)
real(kind=8), allocatable :: SVEC(:,:,:,:) !SVEC(nDirTot,nH,nTempMagn,3)

integer :: IM, I, it, itEnd, J, iH, k, isite, l, n, nP
integer :: iDir, rtob, ibuf, mem_local
real(kind=8) :: cm3tomB
real(kind=8) :: dev, dnrm2_
external :: dev, dnrm2_
logical :: DBG
character(len=15) :: lbl_X, lbl_Y, lbl_Z

DBG = .false.
!Boltz_k = 0.6950356000_wp  ! in cm^-1*K-1
!mu_Bohr = 0.4668643740_wp  ! in cm-1*T-1
cm3tomB = 0.5584938904_wp   ! in cm3 * mol-1 * T

write(6,*)
write(6,'(100A)') (('%'),J=1,96)
write(6,'(40X,A)') 'CALCULATION OF THE MOLAR MAGNETIZATION'
write(6,'(100A)') (('%'),J=1,96)
write(6,*)
!----------------------------------------------------------------------
mem_local = 0
RtoB = 8
if (dbg) write(6,*) 'MAGN:        nM=',nM
if (dbg) write(6,*) 'MAGN:      exch=',exch
if (dbg) write(6,*) 'MAGN:      nLoc=',nLoc
if (dbg) write(6,*) 'MAGN: nTempMagn=',nTempMagn

! Zeeman exchange energy spectrum
call mma_allocate(Wex,nM,'Wex')
call dcopy_(nM,[0.0_wp],0,Wex,1)
mem_local = mem_local+nM*RtoB

! exchange statistical sum, Boltzmann distribution
call mma_allocate(Zex,nTempMagn,'Zex')
call dcopy_(nTempMagn,[0.0_wp],0,Zex,1)
mem_local = mem_local+nTempMagn*RtoB
! spin magnetisation, from the exchange block
call mma_allocate(Sex,3,nTempMagn,'Sex')
call dcopy_(3*nTempMagn,[0.0_wp],0,Sex,1)
mem_local = mem_local+3*nTempMagn*RtoB
! magnetisation, from the exchange block
call mma_allocate(Mex,3,nTempMagn,'Mex')
call dcopy_(3*nTempMagn,[0.0_wp],0,Mex,1)
mem_local = mem_local+3*nTempMagn*RtoB

! local statistical sum, Boltzmann distribution
call mma_allocate(ZL,nneq,nTempMagn,'ZL')
call dcopy_(nneq*nTempMagn,[0.0_wp],0,ZL,1)
mem_local = mem_local+nneq*nTempMagn*RtoB
! spin magnetisation, from the local sites, using ALL states
call mma_allocate(SL,nneq,3,nTempMagn,'SL')
call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,SL,1)
mem_local = mem_local+3*nneq*nTempMagn*RtoB
! magnetisation, from local sites, using ALL states
call mma_allocate(ML,nneq,3,nTempMagn,'ML')
call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,ML,1)
mem_local = mem_local+3*nneq*nTempMagn*RtoB

! local statistical sum, Boltzmann distribution, using only Nexch states
call mma_allocate(ZR,nneq,nTempMagn,'ZR')
call dcopy_(nneq*nTempMagn,[0.0_wp],0,ZR,1)
mem_local = mem_local+nneq*nTempMagn*RtoB
! spin magnetisation, from the local sites, using only Nexch states
call mma_allocate(SR,nneq,3,nTempMagn,'SR')
call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,SR,1)
mem_local = mem_local+3*nneq*nTempMagn*RtoB
! magnetisation, from local sites, using only Nexch states
call mma_allocate(MR,nneq,3,nTempMagn,'MR')
call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,MR,1)
mem_local = mem_local+3*nneq*nTempMagn*RtoB

! Zeeman local energies
call mma_allocate(WL,nneq,nLoc,'WL')
call dcopy_(nneq*nLoc,[0.0_wp],0,WL,1)
mem_local = mem_local+nneq*nLoc*RtoB
! Zeeman local reduced energies, using only Nexch states
call mma_allocate(WR,nneq,nLoc,'WR')
call dcopy_(nneq*nLoc,[0.0_wp],0,WR,1)
mem_local = mem_local+nneq*nLoc*RtoB

! ZRT(nCenter,nTempMagn)
call mma_allocate(ZRT,nCenter,nTempMagn,'ZRT')
call dcopy_(nCenter*nTempMagn,[0.0_wp],0,ZRT,1)
mem_local = mem_local+nCenter*nTempMagn*RtoB
! ZLT(nCenter,nTempMagn)
call mma_allocate(ZLT,nCenter,nTempMagn,'ZLT')
call dcopy_(nCenter*nTempMagn,[0.0_wp],0,ZLT,1)
mem_local = mem_local+nCenter*nTempMagn*RtoB
! MRT(nCenter,3,nTempMagn)
call mma_allocate(MRT,nCenter,3,nTempMagn,'MRT')
call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,MRT,1)
mem_local = mem_local+3*nCenter*nTempMagn*RtoB
! MLT(nCenter,3,nTempMagn)
call mma_allocate(MLT,nCenter,3,nTempMagn,'MLT')
call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,MLT,1)
mem_local = mem_local+3*nCenter*nTempMagn*RtoB
! SRT(nCenter,3,nTempMagn)
call mma_allocate(SRT,nCenter,3,nTempMagn,'SRT')
call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,SRT,1)
mem_local = mem_local+3*nCenter*nTempMagn*RtoB
! SLT(nCenter,3,nTempMagn)
call mma_allocate(SLT,nCenter,3,nTempMagn,'SLT')
call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,SLT,1)
mem_local = mem_local+3*nCenter*nTempMagn*RtoB

! total statistical sum, Boltzmann distribution
call mma_allocate(ZT,nH,nTempMagn,'ZT')
call dcopy_(nH*nTempMagn,[0.0_wp],0,ZT,1)
mem_local = mem_local+nH*nTempMagn*RtoB
! total spin magnetisation
call mma_allocate(ST,3,nH,nTempMagn,'ST')
call dcopy_(3*nH*nTempMagn,[0.0_wp],0,ST,1)
mem_local = mem_local+3*nH*nTempMagn*RtoB
! total magnetisation
call mma_allocate(MT,3,nH,nTempMagn,'MT')
call dcopy_(3*nH*nTempMagn,[0.0_wp],0,MT,1)
mem_local = mem_local+3*nH*nTempMagn*RtoB
! total spin magnetisation
call mma_allocate(SAV,nH,nTempMagn,'SAV')
call dcopy_(nH*nTempMagn,[0.0_wp],0,SAV,1)
mem_local = mem_local+nH*nTempMagn*RtoB
! total magnetisation
call mma_allocate(MAV,nH,nTempMagn,'MAV')
call dcopy_(nH*nTempMagn,[0.0_wp],0,MAV,1)
mem_local = mem_local+nH*nTempMagn*RtoB
! total spin magnetisation vector
call mma_allocate(SVEC,nDirTot,nH,nTempMagn,3,'SVEC')
call dcopy_(nDirTot*nH*nTempMagn*3,[0.0_wp],0,SVEC,1)
mem_local = mem_local+nDirTot*nH*nTempMagn*3*RtoB
! total magnetisation vector
call mma_allocate(MVEC,nDirTot,nH,nTempMagn,3,'MVEC')
call dcopy_(nDirTot*nH*nTempMagn*3,[0.0_wp],0,MVEC,1)
mem_local = mem_local+nDirTot*nH*nTempMagn*3*RtoB

! orientation of the field
call mma_allocate(dHX,nDirTot,'dHX')
call mma_allocate(dHY,nDirTot,'dHY')
call mma_allocate(dHZ,nDirTot,'dHZ')
call mma_allocate(dHW,nDirTot,'dHW')
call dcopy_(nDirTot,[0.0_wp],0,dHX,1)
call dcopy_(nDirTot,[0.0_wp],0,dHY,1)
call dcopy_(nDirTot,[0.0_wp],0,dHZ,1)
call dcopy_(nDirTot,[0.0_wp],0,dHW,1)
mem_local = mem_local+4*nDirTot*RtoB

call mma_allocate(H,nH,'H  field')
call dcopy_(nH,[0.0_wp],0,H,1)
mem_local = mem_local+nH*RtoB
if (dbg) write(6,*) 'MAGN:  memory allocated (local):'
if (dbg) write(6,*) 'mem_local=',mem_local
if (dbg) write(6,*) 'MAGN:  memory allocated (total):'
if (dbg) write(6,*) 'mem_total=',mem+mem_local
!-----------------------------------------------------------------------
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
nP = get_nP(nsymm,ngrid)

call hdir(nDir,nDirZee,dirX,dirY,dirZ,dir_weight,nP,nsymm,ngrid,nDirTot,dHX,dHY,dHZ,dHW)
call Add_Info('MR_MAGN  dHX',dHX(1:5),5,10)
call Add_Info('MR_MAGN  dHY',dHY(1:5),5,10)
call Add_Info('MR_MAGN  dHZ',dHZ(1:5),5,10)
call Add_Info('MR_MAGN  dWX',dHW(1:5),5,10)
if (DBG) then
  write(6,'(A,F10.6)') '        zJ = ',zJ
  write(6,'(A,F10.6)') '      thrs = ',thrs
  write(6,'(A,   I6)') '      iopt = ',iopt
  write(6,'(A, 10I6)') '     nss(i)=',(nss(i),i=1,nneq)
  write(6,'(A, 10I6)') '   nexch(i)=',(nexch(i),i=1,nneq)
  write(6,'(A,   I6)') '  nTempMagn= ',nTempMagn
  write(6,'(A      )') 'TempMagn(i)='
  write(6,'(10F10.6)') (TempMagn(i),i=1,nTempMagn)
  write(6,*) 'm_paranoid =',m_paranoid
  write(6,*) 'm_accurate =',m_accurate
  write(6,*) '     smagn =',smagn
  do k=1,nneq
    write(6,'(A,i6)') 'site ',k
    write(6,'(A,i2,A)') 'ESO(',k,')'
    write(6,'(10F12.6)') (ESO(k,i),i=1,nss(k))
    write(6,'(A,i2,A)') 'S_SO(',k,')'
    do i=1,nss(k)
      do j=1,nss(k)
        write(6,'(3(2F15.10,3x))') (s_so(k,l,i,j),l=1,3)
      end do
    end do
    write(6,'(A,i2,A)') 'DIPSO(',k,')'
    do i=1,nss(k)
      do j=1,nss(k)
        write(6,'(3(2F15.10,3x))') (dipso(k,l,i,j),l=1,3)
      end do
    end do
  end do !k
end if

write(6,'(2X,A,i3,A)') 'Molar magnetization will be calculated in ',nH,' points, equally distributed in magnetic field range'
write(6,'(2X,F4.1,1x,a,1X,F4.1,a,5(F6.3,a))') HMIN,'--',HMAX,' T., at the following temperatures:'
do i=1,nTempMagn,10
  j = min(nTempMagn,i+9)
  write(6,'(10(F8.4,A))') (TempMagn(k),' K.;',k=i,j)
end do
write(6,'(2X,A,I4,A)') 'Powder molar magnetization will be averaged on ',nP,' directions of the applied magnetic field.'
write(6,'(2x,10A)') ('--------',i=1,10)
if (nsymm == 1) then
  write(6,'(23x,A)') 'Lebedev-Laikov grid on a hemisphere:'
  write(6,'(38x,A)') 'z >= 0;'
else if (nsymm == 2) then
  write(6,'(23x,A)') 'Lebedev-Laikov grid on a 4th-of-a-sphere:'
  write(6,'(34x,A)') 'x >= 0; z >= 0;'
else if (nsymm == 3) then
  write(6,'(23x,A)') 'Lebedev-Laikov grid on a 8th-of-a-sphere:'
  write(6,'(30x,A)') 'x >= 0; y >= 0; z >= 0;'
end if
write(6,'(2x,10A)') ('--------',i=1,10)
write(6,'(2x,A,12x,A,2(18x,A),16x,A)') 'Nr.','x','y','z','weight'
do i=1,nP
  write(6,'(i4,2x,4(F18.12,1x))') i,dHX(i+nDir+nDirZee),dHY(i+nDir+nDirZee),dHZ(i+nDir+nDirZee),dHW(i+nDir+nDirZee)
end do
write(6,'(2x,10A)') ('--------',i=1,10)
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
if (m_accurate) then
  write(6,'(2x,A)') 'The contribution of local excited states  is computed exactly.'
  if ((zJ /= 0.0_wp) .and. m_paranoid .and. (nTempMagn > 1)) then
    write(6,'(2x,A)') 'The average spin is computed exactly for each temperature point.'
  else if ((zJ /= 0.0_wp) .and. (.not. m_paranoid) .and. (nTempMagn > 1)) then
    write(6,'(2x,A,F9.3)') 'The average spin is computed exactly only for the temperature point T= ',maxval(TempMagn(:))
    write(6,'(2x,A     )') 'We consider this to be a good approximation.'
  end if
else
  if ((zJ /= 0.0_wp) .and. m_paranoid .and. (nTempMagn > 1)) then
    write(6,'(2x,A)') 'The average spin is computed exactly for each temperature point.'
  else if ((zJ /= 0.0_wp) .and. (.not. m_paranoid) .and. (nTempMagn > 1)) then
    write(6,'(2x,A,F9.3)') 'The average spin is computed exactly only for the temperature  point T= ',TempMagn(1)
    write(6,'(2x,A)') 'We consider this to be a good approximation.'
    write(6,'(2x,A)') 'Use MPAR keyword to compute average spin exactly, for each temperature point, in case higher accuracy '// &
                      'is needed.'
  end if
  write(6,'(2x,A)') 'The contribution of local excited states  is computed approximately (by using the susceptibility data). '
end if
if (compute_Mdir_vector) then
  write(6,'(2x,A,i2,a)') 'The magnetization vector for ',nDir,' directions of the applied magnetic field will be calculated.'
else
  write(6,'(2X,A)') 'The magnetization vector was not calculated.'
end if
if (zeeman_energy) then
  write(6,'(2x,A,i2,a)') 'The Zeeman splitting for ',nDirZee,' directions of the applied magnetic field will be calculated.'
  write(6,'(2x,a)') 'The Zeeman energies for each direction of the applied magnetic field are written in files "energy_XX.txt".'
  if (DBG) then
    do i=1,nDirZee
      write(6,'(A,I4,A,I4)') 'LuZee(',i,' )=',LUZee(i)
    end do
  end if
else
  write(6,'(2X,A)') 'Computation of the Zeeman splitting was not requested.'
end if
!if (maxval(TempMagn) < W(exch)) then
!  write(6,'(2x,A)') 'Contribution to molar magnetization coming from local excited states is taken into account.'
!else
!  write(6,'(2x,A)') 'Contribution to molar magnetization coming from local excited states is NOT taken into account.'
!  write(6,'(2x,A)') 'TMAG is requesting to compute magnetization at a larger temperature than the highest exchange state.'
!  write(6,'(2x,A)') 'Please include more states into the exchange coupling.'
!end if

! /// opening the loop over the field points
do iH=1,nH
  ! /// ----------------------------------------------------------------
  if (HINPUT) then
    H(iH) = HEXP(iH)
    if (H(iH) == 0.0_wp) then
      H(iH) = 0.0001_wp
    end if
  else
    DLTH = (HMAX-HMIN)/dble(nH-1)
    if (iH == 1) then
      H(iH) = HMIN+0.0001_wp
    else
      H(iH) = HMIN+DLTH*dble(iH-1)
    end if
    if (H(iH) == 0.0_wp) then
      H(iH) = 0.0001_wp
    end if
  end if

  if (DBG) write(6,'(A,i0,A,F10.5)') 'MAGNETIZATION::  H(',iH,') = ',H(iH)

  ! ///  opening the loop over different directions of the magnetic field
  do IM=1,NDIRTOT
    call dcopy_(nM,[0.0_wp],0,Wex,1)
    call dcopy_(nneq*nLoc,[0.0_wp],0,WL,1)
    call dcopy_(nneq*nLoc,[0.0_wp],0,WR,1)
    call dcopy_(nTempMagn,[0.0_wp],0,Zex,1)
    call dcopy_(nneq*nTempMagn,[0.0_wp],0,ZL,1)
    call dcopy_(nneq*nTempMagn,[0.0_wp],0,ZR,1)
    call dcopy_(3*nTempMagn,[0.0_wp],0,Sex,1)
    call dcopy_(3*nTempMagn,[0.0_wp],0,Mex,1)
    call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,SL,1)
    call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,ML,1)
    call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,SR,1)
    call dcopy_(3*nneq*nTempMagn,[0.0_wp],0,MR,1)
    if (DBG) write(6,'(A,F20.10)') 'zJ = ',zJ
    ! exchange magnetization:
    call MAGN(exch,NM,dHX(iM),dHY(iM),dHZ(iM),H(iH),W,zJ,THRS,DIPEXCH,S_EXCH,nTempMagn,TempMagn,smagn,Wex,Zex,Sex,Mex,m_paranoid, &
              DBG)
    if (DBG) write(6,'(A,3F11.7)') 'MEX:',(Mex(l,1),l=1,3)
    if (IM == 1) then
      call Add_Info('MEX_MAGN    ',[dnrm2_(3*nTempMagn,Mex,1)],1,8)
      call Add_Info('MR_MAGN  Wex',[dnrm2_(nM,Wex,1)],1,8)
    end if
    ! compute local magnetizations:
    if (m_accurate) then
      do i=1,nneq
        ! all states:
        if (NSS(i) > NEXCH(i)) then
          ! this check is to avoid the unnecessary computation, in cases when no local excited states are present
          call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),H(iH),ESO(i,1:NSS(i)),zJ,THRS,DIPSO(i,1:3,1:NSS(i),1:NSS(i)), &
                    S_SO(i,1:3,1:NSS(i),1:NSS(i)),nTempMagn,TempMagn(1:nTempMagn),smagn,WL(i,1:NEXCH(i)),ZL(i,1:nTempMagn), &
                    SL(i,1:3,1:nTempMagn),ML(i,1:3,1:nTempMagn),m_paranoid,DBG)
          if (IM == 2) then
            call Add_Info('MR_MAGN  WL',[dnrm2_(nexch(i),WL,1)],1,8)
          end if

          if (DBG) write(6,'(A,I2,A,3F11.7)') 'ML: site',i,' : ',(ML(i,l,1),l=1,3)
          ! only local "exchange states":
          call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),H(iH),ESO(i,1:NEXCH(i)),zJ,THRS,DIPSO(i,1:3,1:NEXCH(i),1:NEXCH(i)), &
                    S_SO(i,1:3,1:NEXCH(i),1:NEXCH(i)),nTempMagn,TempMagn,smagn,WR(i,1:Nexch(i)),ZR(i,1:nTempMagn), &
                    SR(i,1:3,1:nTempMagn),MR(i,1:3,1:nTempMagn),m_paranoid,DBG)
          call Add_Info('MR_MAGN  WR',[dnrm2_(nexch(i),WR,1)],1,8)
          if (DBG) write(6,'(A,I2,A,3F11.7)') 'MR: site',i,' : ',(MR(i,l,1),l=1,3)
        end if
      end do

      if (IM == 3) then
        call Add_Info('ML_MAGN',[dnrm2_(3*nTempMagn*nneq,ML,1)],1,8)
        call Add_Info('MR_MAGN',[dnrm2_(3*nTempMagn*nneq,MR,1)],1,8)
      end if

      ! expand the basis and rotate local vectors to the general
      ! coordinate system:
      call dcopy_(nCenter*nTempMagn,[0.0_wp],0,ZRT,1)
      call dcopy_(nCenter*nTempMagn,[0.0_wp],0,ZLT,1)
      call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,MRT,1)
      call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,MLT,1)
      call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,SRT,1)
      call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,SLT,1)
      isite = 0
      do i=1,NNEQ
        do j=1,NEQ(i)
          isite = isite+1
          ! statistical distributions
          do iT=1,nTempMagn
            ZLT(isite,iT) = ZL(i,iT)
            ZRT(isite,iT) = ZR(i,iT)
          end do
         ! magnetizations:
         !  use R_rot matrices, which have determinant +1.
         !  note that  R_lg matrices have arbitrary determinant.
          do iT=1,nTempMagn
            do l=1,3
              do n=1,3
                MLT(isite,l,iT) = MLT(isite,l,iT)+r_rot(i,j,l,n)*ML(i,n,iT)
                SLT(isite,l,iT) = SLT(isite,l,iT)+r_rot(i,j,l,n)*SL(i,n,iT)
                MRT(isite,l,iT) = MRT(isite,l,iT)+r_rot(i,j,l,n)*MR(i,n,iT)
                SRT(isite,l,iT) = SRT(isite,l,iT)+r_rot(i,j,l,n)*SR(i,n,iT)
              end do
            end do
          end do
        end do ! j, neq(i)
      end do ! i, nneq
    end if ! m_accurate

    ! compute the total magnetizations according to the derived formulas:
    if (m_accurate) then
      do iT=1,nTempMagn
        if (smagn) then
          call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLT(:,iT),SRT(:,:,iT),ZRT(:,iT),iopt,ST(:,iH,iT),ZT(iH,iT))
        end if
        call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLT(:,iT),MRT(:,:,iT),ZRT(:,iT),iopt,MT(:,iH,iT),ZT(iH,iT))
      end do
    else
      ! add the contribution from local excited states using the approximate
      ! X*H expression:
      call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,SRT,1)
      call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,SLT,1)
      do iT=1,nTempMagn
        do isite=1,nCenter
          do l=1,3
            MRT(isite,l,iT) = (XRM(isite,iT,l,1)*H(iH)*dHX(iM)+XRM(isite,iT,l,2)*H(iH)*dHY(iM)+XRM(isite,iT,l,3)*H(iH)*dHZ(iM))/ &
                              cm3tomB
            MLT(isite,l,iT) = (XLM(isite,iT,l,1)*H(iH)*dHX(iM)+XLM(isite,iT,l,2)*H(iH)*dHY(iM)+XLM(isite,iT,l,3)*H(iH)*dHZ(iM))/ &
                              cm3tomB
          end do
        end do
        if (smagn) then
          call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLM(:,iT),SRT(:,:,iT),ZRM(:,iT),iopt,ST(:,iH,iT),ZT(iH,iT))
        end if
        call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLM(:,iT),MRT(:,:,iT),ZRM(:,iT),iopt,MT(:,iH,iT),ZT(iH,iT))
      end do
    end if
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! print out hte Zeeman eigenstates
    if (zeeman_energy) then
      if ((iH == 1) .and. (iM == nDir+1)) write(6,'(A)') 'Energies of the Zeeman Hamiltonian for the following directions of '// &
                                                         'the applied field:'
      if ((iH == 1) .and. (iM > nDir) .and. (iM <= nDir+nDirZee)) then
        write(6,'(A,I3,A,3F10.6,3x,5A)') 'direction Nr.',iM-nDir,' : ',dHX(iM),dHY(iM),dHZ(iM),'written in file "zeeman_energy_', &
                                         char(48+mod(int((iM-nDir)/100),10)),char(48+mod(int((iM-nDir)/10),10)), &
                                         char(48+mod((iM-nDir),10)),'.txt".'
        write(LUZee(iM-nDir),'(A,3F24.15)') '# direction of the applied magnetic field:',dHX(iM),dHY(iM),dHZ(iM)
        write(LuZee(iM-nDir),'(A,6x,A,1000(I4,6x) )') '# H(T)',' State =>',(i,i=1,nm)
      end if

      if ((iM > nDir) .and. (iM <= nDir+nDirZee)) then
        write(LUZee(iM-nDir),'(F8.4,1000F10.3)') H(IH),(Wex(I),I=1,NM)
      end if
    end if !zeeman_energy
    ! ------------------------------------------------------------------
    ! computing the AVERAGE MOMENTS calculated at different temperatures
    ! (TempMagn(i))
    do iT=1,nTempMagn
      do l=1,3
        MVEC(iM,iH,iT,l) = MT(l,iH,iT)
        SVEC(iM,iH,iT,l) = ST(l,iH,iT)
      end do !l

      !MAV(iH,iT) = MAV(iH,iT)+MT(1,iH,iT)*dHX(iM)*dHW(iM)+MT(2,iH,iT)*dHY(iM)*dHW(iM)+MT(3,iH,iT)*dHZ(iM)*dHW(iM)
      !SAV(iH,iT) = SAV(iH,iT)+ST(1,iH,iT)*dHX(iM)*dHW(iM)+ST(2,iH,iT)*dHY(iM)*dHW(iM)+ST(3,iH,iT)*dHZ(iM)*dHW(iM)
      ! accumulate contributions:
      call daxpy_(1,dHX(iM)*dHW(iM),MT(1,iH,iT),1,MAV(iH,iT),1)
      call daxpy_(1,dHY(iM)*dHW(iM),MT(2,iH,iT),1,MAV(iH,iT),1)
      call daxpy_(1,dHZ(iM)*dHW(iM),MT(3,iH,iT),1,MAV(iH,iT),1)
      call daxpy_(1,dHX(iM)*dHW(iM),ST(1,iH,iT),1,SAV(iH,iT),1)
      call daxpy_(1,dHY(iM)*dHW(iM),ST(2,iH,iT),1,SAV(iH,iT),1)
      call daxpy_(1,dHZ(iM)*dHW(iM),ST(3,iH,iT),1,SAV(iH,iT),1)
    end do !iT
    ! ///  closing the loops over field strengths and directions
  end do ! iM
end do ! iH
! Close Zeeman files, if opened
if (Zeeman_Energy) then
  do i=1,nDirZee
    close(LUZee(i))
  end do
end if

! -------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
! -------------------------------------------------------------------
if (smagn) then
  do iT=1,nTempMagn
    !write(6,*)
    do iDir=1,nDir+nDirZee
      !write(6,*)
      write(6,'(A,A,1x,A)') '--------|','------------------------------------------------------------|', &
                            '|------------------------------------------------------------|'
      write(6,'(A,i3,26x,A,1x,A,60x,A)') 'Direction of the applied magnetic field:',iDir,'|','|','|'
      write(6,'(A,F18.14,44x,A,1x,A,60x,A)') 'proj X=',dHX(iDIR),'|','|','|'
      write(6,'(A,F18.14,44x,A,1x,A,60x,A)') 'proj Y=',dHY(iDir),'|','|','|'
      write(6,'(A,F18.14,44x,A,1x,A,60x,A)') 'proj Z=',dHZ(iDir),'|','|','|'
      write(6,'(A,F7.4,A,41x,A,1x,A,60x,A)') 'Temperature = ',TempMagn(iT),' kelvin','|','|','|'
      write(6,'(A,A,1x,A)') '--------|','------------------------------------------------------------|', &
                            '|------------------------------------------------------------|'
      write(6,'(2x,A,12x,2A,1x,A,10x,2A)') 'Field |','Magnetization Vector            |','   Total Magn. |','|', &
                                           'Spin Magnetization Vector         |','   Total Magn. |'
      write(6,'(5A,1x,5A)') '--------|','--- proj X ---|','--- proj Y ---|','--- proj Z ---|','- in this dir.-|', &
                            '|--- proj X ---|','--- proj Y ---|','--- proj Z ---|','- in this dir.-|'
      do iH=1,nH
        write(6,'(F7.3,1x,A,3(ES13.6,1x,A),ES14.7,1x,A,1x,A,3(ES13.6,1x,A),ES14.7,1x,A)') &
          H(iH),'|',MVEC(iDir,iH,iT,1),' ',MVEC(iDir,iH,iT,2),' ',MVEC(iDir,iH,iT,3),'|', &
          (MVEC(iDir,iH,iT,1)*dHX(iDir)+MVEC(iDir,iH,iT,2)*dHY(iDir)+MVEC(iDir,iH,iT,3)*dHZ(iDir)),'|','|',SVEC(iDir,iH,iT,1),' ', &
          SVEC(iDir,iH,iT,2),' ',SVEC(iDir,iH,iT,3),'|', &
          (SVEC(iDir,iH,iT,1)*dHX(iDir)+SVEC(iDir,iH,iT,2)*dHY(iDir)+SVEC(iDir,iH,iT,3)*dHZ(iDir)),'|'
      end do
      write(6,'(A,A,1x,A)') '--------|','------------------------------------------------------------|', &
                            '|------------------------------------------------------------|'
    end do !iDir
  end do !iT

else

  do iT=1,nTempMagn
    do iDir=1,nDir+nDirZee
      write(6,'(2A)') '--------|','------------------------------------------------------------|'
      write(6,'(A,i3,26x,A)') 'Direction of the applied magnetic field:',iDir,'|'
      write(6,'(A,F18.14,44x,A)') 'proj X=',dHX(iDIR),'|'
      write(6,'(A,F18.14,44x,A)') 'proj Y=',dHY(iDir),'|'
      write(6,'(A,F18.14,44x,A)') 'proj Z=',dHZ(iDir),'|'
      write(6,'(A,F7.4,A,41x,A)') 'Temperature = ',TempMagn(iT),' kelvin','|'
      write(6,'(2A)') '--------|','------------------------------------------------------------|'
      write(6,'(2x,A,12x,2A)') 'Field |','Magnetization Vector            |','   Total Magn. |'
      write(6,'(5A)') '--------|','--- proj X ---|','--- proj Y ---|','--- proj Z ---|','- in this dir.-|'
      do iH=1,nH
        write(6,'(F7.3,1x,A,3(ES13.6,1x,A),ES14.7,1x,A)') &
          H(iH),'|',MVEC(iDir,iH,iT,1),' ',MVEC(iDir,iH,iT,2),' ',MVEC(iDir,iH,iT,3),'|', &
          (MVEC(iDir,iH,iT,1)*dHX(iDir)+MVEC(iDir,iH,iT,2)*dHY(iDir)+MVEC(iDir,iH,iT,3)*dHZ(iDir)),'|'
      end do
      write(6,'(2A)') '--------|','------------------------------------------------------------|'
    end do !iDir

  end do !iT
end if !(smagn)
! ---------------------------------------------------------------------
! COMPUTING THE STANDARD DEVIATION OF THE MAGNETIZATION
!if (HINPUT) then
!  do iT=1,nTempMagn
!    STDEV(iT) = 0.0_wp
!    STDEV(iT) = dev(nH,MAV(:,iT),Mexp(:,iT))
!  end do
!end if

write(6,*)
write(6,'(15X,A)') 'HIGH-FIELD POWDER MAGNETIZATION'
write(6,'(20X,A)') '(Units: Bohr magneton)'
write(6,*)
do iT=1,nTempMagn,5
  iTEnd = min(nTempMagn,iT+4)
  write(6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  write(6,'(A,11(F10.3,A))') '    H(T)   |STATISTICAL SUM|',(TempMagn(i),' K.  |',i=iT,iTEnd)
  write(6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  do iH=1,nH
    write(6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))') H(iH),'|',ZT(iH,1),'|',(MAV(iH,i),'|',i=iT,iTEnd)
  end do
  write(6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  if (HINPUT) then
    write(6,'(A,15x,A, 11(f14.10,1x,A) )') 'ST.DEV.M   |','|',(dev(nH,MAV(:,i),Mexp(:,i)),'|',i=iT,iTEnd)
    write(6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  end if
end do
if (smagn) then
  write(6,*)
  write(6,'(15X,A)') 'HIGH-FIELD POWDER SPIN MAGNETIZATION'
  write(6,'(20X,A)') '(Units: Bohr magneton)'
  write(6,*)
  do iT=1,nTempMagn,5
    iTEnd = min(nTempMagn,iT+4)
    write(6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
    write(6,'(A,11(F10.3,A))') '   H(T)   |STATISTICAL SUM|',(TempMagn(i),' K.  |',i=iT,iTEnd)
    write(6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
    do iH=1,nH
      write(6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))') H(iH),'|',ZT(iH,1),'|',(SAV(iH,i),'|',i=iT,iTEnd)
    end do
    write(6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  end do
end if !smagn
! add some verification:
call Add_Info('MR_MAGN',[dnrm2_(3*nTempMagn*nneq,MR,1)],1,8)
call Add_Info('MR_MAGN',[dnrm2_(3*nTempMagn*nneq,MR,1)],1,8)
call Add_Info('MR_MAGN',[dnrm2_(3*nTempMagn*nneq,MR,1)],1,8)

call Add_Info('H_MAGN       ',[dnrm2_(nH,H,1)],1,6)
call Add_Info('MAGN_AVERAGED',[dnrm2_(nH*nTempMagn,MAV,1)],1,6)
call Add_Info('ZTL_MAGN',[dnrm2_(nH*nTempMagn,ZT,1)],1,8)
do iH=1,nH
  write(lbl_X,'(A,i3)') 'MAGN_VECT X ',iH
  write(lbl_Y,'(A,i3)') 'MAGN_VECT Y ',iH
  write(lbl_Z,'(A,i3)') 'MAGN_VECT Z ',iH
  ibuf = nDirTot*nTempMagn
  call Add_Info(lbl_X,[dnrm2_(ibuf,MVEC(:,iH,:,1),1)],1,8)
  call Add_Info(lbl_Y,[dnrm2_(ibuf,MVEC(:,iH,:,2),1)],1,8)
  call Add_Info(lbl_Z,[dnrm2_(ibuf,MVEC(:,iH,:,3),1)],1,8)
end do

!-------------------------  PLOTs -------------------------------------!
if (DoPlot) then
  if (hinput) then
    call plot_MH_with_Exp(nH,H,nTempMagn,TempMagn,MAV,Mexp)
  else
    call plot_MH_no_Exp(nH,H,nTempMagn,TempMagn,MAV)
  end if

  !if (zeeman_energy) call plot_zeeman(nH,nM,nDirZee,H,LuZee )
end if
!------------------------- END PLOTs ----------------------------------!

!-----------------------------------------------------------------------
call mma_deallocate(Wex)
call mma_deallocate(Zex)
call mma_deallocate(Sex)
call mma_deallocate(Mex)
call mma_deallocate(ZL)
call mma_deallocate(SL)
call mma_deallocate(ML)
call mma_deallocate(ZR)
call mma_deallocate(SR)
call mma_deallocate(MR)
call mma_deallocate(WL)
call mma_deallocate(WR)
call mma_deallocate(ZRT)
call mma_deallocate(ZLT)
call mma_deallocate(MRT)
call mma_deallocate(MLT)
call mma_deallocate(SRT)
call mma_deallocate(SLT)
call mma_deallocate(ZT)
call mma_deallocate(ST)
call mma_deallocate(MT)
call mma_deallocate(SAV)
call mma_deallocate(MAV)
call mma_deallocate(SVEC)
call mma_deallocate(MVEC)
call mma_deallocate(dHX)
call mma_deallocate(dHY)
call mma_deallocate(dHZ)
call mma_deallocate(dHW)

call mma_deallocate(H)

return

end subroutine magnetization_pa
