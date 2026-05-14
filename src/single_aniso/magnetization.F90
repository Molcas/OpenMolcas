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

subroutine magnetization(nss,nM,nTempMagn,nDirTot,nDir,nDirZee,nH,iPrint,LUZee,mem,nsymm,ngrid,compute_Mdir_vector,zeeman_energy, &
                         hinput,m_paranoid,smagn,doplot,TempMagn,eso,dirX,dirY,dirZ,dir_weight,hexp,magn_exp,zJ,hmin,hmax,EM,thrs, &
                         dipm,sm,dbg)
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
!               scalar real                                            *
!     EM      : cut-off energy (energy of the last s-o state which is  *
!               included in the Zeeman matrix                          *
!               scalar real                                            *
!     IFINAL  : integer                                                *
!               termination flag                                       *
!                                                                      *
!     local variables:                                                 *
!     MVEC, SVEC : magnetization and spin vector                       *
!     MAV, SAV   : average powder M and S                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     Liviu Ungur                                                      *
!     University of Leuven, Belgium, 2008-2017                         *
!                                                                      *
!----------------------------------------------------------------------*
!     History:                                                         *
!     Liviu Ungur, 2008-2017 various modifications                     *
!***********************************************************************

use Lebedev_quadrature, only: order_table
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: nss, nM, nTempMagn, nDirTot, nDir, nDirZee, nH, iprint, LUZee(nDirZee), mem, nsymm, ngrid
logical(kind=iwp), intent(in) :: compute_Mdir_vector, zeeman_energy, hinput, m_paranoid, smagn, DoPlot, dbg
real(kind=wp), intent(in) :: TempMagn(nTempMagn), eso(nss), dirX(nDir), dirY(nDir), dirZ(nDir), dir_weight(nDirZee,3), hexp(nH), &
                             magn_exp(nH,nTempMagn), zj, hmin, hmax, EM, thrs
complex(kind=wp), intent(in) :: dipm(3,nss,nss), sm(3,nss,nss)
integer(kind=iwp) :: I, IDIR, IH, IM, iT, iTemp, iTEnd, J, L, mem_local, nP
real(kind=wp) :: DLTH, mv, sv
character(len=99) :: STLNE2
real(kind=wp), allocatable :: dHW(:), dHX(:), dHY(:), dHZ(:), H(:), MAV(:,:), MT(:,:,:), MT_TMP(:,:), MVEC(:,:,:,:), SAV(:,:), &
                              ST(:,:,:), ST_TMP(:,:), STDEV(:), SVEC(:,:,:,:), WM(:), ZT(:,:), ZT_TMP(:)
real(kind=wp), external :: dev

!-----------------------------------------------------------------------
! Allocate necessary memory
mem_local = 0

! Zeeman exchange energy spectrum
call mma_allocate(WM,nM,'W')
mem_local = mem_local+size(WM)*RtoB

call mma_allocate(MT,3,nH,nTempMagn,'MT')
call mma_allocate(MT_TMP,3,nTempMagn,'MT_TMP')
mem_local = mem_local+size(MT)*RtoB
mem_local = mem_local+size(MT_TMP)*RtoB

call mma_allocate(ST,3,nH,nTempMagn,'ST')
call mma_allocate(ST_TMP,3,nTempMagn,'ST_TMP')
mem_local = mem_local+size(ST)*RtoB
mem_local = mem_local+size(ST_TMP)*RtoB

call mma_allocate(MAV,nH,nTempMagn,'MAV')
MAV(:,:) = Zero
mem_local = mem_local+size(MAV)*RtoB

call mma_allocate(SAV,nH,nTempMagn,'SAV')
SAV(:,:) = Zero
mem_local = mem_local+size(SAV)*RtoB

call mma_allocate(ZT,nH,nTempMagn,'ZT')
call mma_allocate(ZT_TMP,nTempMagn,'ZT_TMP')
mem_local = mem_local+size(ZT)*RtoB
mem_local = mem_local+size(ZT_TMP)*RtoB

call mma_allocate(MVEC,nDirTot,nH,nTempMagn,3,'MVEC')
call mma_allocate(SVEC,nDirTot,nH,nTempMagn,3,'SVEC')
mem_local = mem_local+size(MVEC)*RtoB
mem_local = mem_local+size(SVEC)*RtoB

call mma_allocate(H,nH,'H')
mem_local = mem_local+size(H)*RtoB

call mma_allocate(STDEV,nTempMagn,'STDEV')
mem_local = mem_local+size(STDEV)*RtoB

call mma_allocate(dHX,nDirTot,'dHX')
call mma_allocate(dHY,nDirTot,'dHY')
call mma_allocate(dHZ,nDirTot,'dHZ')
call mma_allocate(dHW,nDirTot,'dHW')
mem_local = mem_local+size(dHX)*RtoB
mem_local = mem_local+size(dHY)*RtoB
mem_local = mem_local+size(dHZ)*RtoB
mem_local = mem_local+size(dHW)*RtoB
if (dbg) then
  write(u6,*) 'MAGNETIZATION:  memory allocated (local):'
  write(u6,*) 'mem_local=',mem_local
  write(u6,*) 'MAGNETIZATION:  memory allocated (total):'
  write(u6,*) 'mem_total=',mem+mem_local
end if
!-----------------------------------------------------------------------
write(u6,*)
write(u6,'(A)') repeat('%',96)
write(u6,'(40X,A)') 'CALCULATION OF THE MOLAR MAGNETIZATION'
write(u6,'(A)') repeat('%',96)
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

nP = order_table(nsymm,ngrid)

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

  !/// opening the loop over different directions of the magnetic field
  do iM=1,NDIRTOT
    ! Entry into monitor: Status line
    write(STLNE2,'(A,I4,A,I4,A,I4,A,I4)') 'Field: ',IH,' from ',nH,' at direction ',IM,' from ',NDIRTOT
    call StatusLine('SINGLE_ANISO: powder magnetization: ',STLNE2)
    ! actual calculation of the MT and ST, ZT
    call MAGN(NSS,NM,dHX(iM),dHY(iM),dHZ(iM),H(iH),ESO,zJ,THRS,DIPM,SM,nTempMagn,TempMagn,smagn,WM,ZT_TMP,ST_TMP,MT_TMP, &
              m_paranoid,DBG)
    ZT(iH,:) = ZT_TMP(:)
    ST(:,iH,:) = ST_TMP(:,:)
    MT(:,iH,:) = MT_TMP(:,:)
    if (DBG .and. (iH == nH) .and. (iM == 23)) then
      write(u6,'(A,3ES16.8)') 'iM:',dHX(iM),dHY(iM),dHZ(iM)
      write(u6,'(2(A,3ES16.8,1x),A,ES16.8)') 'MT:',(MT(l,iH,1),l=1,3),'ST:',(ST(l,iH,1),l=1,3),'ZSTAT:',ZT(iH,1)
      write(u6,'(A,3ES16.8)') 'WM:',(WM(l),l=1,nM)
    end if
    !-------------------------------------------------------------------
    if (zeeman_energy) then
      if ((iH == 1) .and. (iM == nDir+1)) &
        write(u6,'(A)') 'Energies of the Zeeman Hamiltonian for the following directions of the applied field:'
      if ((iH == 1) .and. (iM > nDir) .and. (iM <= nDir+nDirZee)) then
        write(u6,'(A,I3,A,3F10.6,3x,5A)') 'direction Nr.',iM-nDir,' : ',dHX(iM),dHY(iM),dHZ(iM),'written in file "zeeman_energy_', &
                                          char(48+mod((iM-nDir)/100,10)),char(48+mod((iM-nDir)/10,10)), &
                                          char(48+mod(iM-nDir,10)),'.txt".'

        write(LUZee(iM-nDir),'(A,3F24.15)') '# direction of the applied magnetic field:',dHX(iM),dHY(iM),dHZ(iM)
        write(LUZee(iM-nDir),'(A,6x,A,1000(I4,6x) )') '# H(T)',' State =>',(i,i=1,nm)
      end if

      if ((iM > nDir) .and. (iM <= nDir+nDirZee)) write(LUZee(iM-nDir),'(F8.4,1000F10.3)') H(IH),(WM(I),I=1,NM)
    end if !zeeman_energy
    !-------------------------------------------------------------------
    ! computing the AVERAGE MOMENTS calculated at different temperatures
    ! (TempMagn(i))
    do iT=1,nTempMagn
      MVEC(iM,iH,iT,:) = MT(:,iH,iT)
      SVEC(iM,iH,iT,:) = ST(:,iH,iT)

      if (iM > nDir+nDirZee) then
        ! accumulate contributions:
        MAV(iH,iT) = MAV(iH,iT)+dHW(iM)*(dHX(iM)*MT(1,iH,iT)+dHY(iM)*MT(2,iH,iT)+dHZ(iM)*MT(3,iH,iT))
        SAV(iH,iT) = SAV(iH,iT)+dHW(iM)*(dHX(iM)*ST(1,iH,iT)+dHY(iM)*ST(2,iH,iT)+dHZ(iM)*ST(3,iH,iT))
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

call Add_Info('MAGN_AVERAGED',MAV,nH*nTempMagn,5)
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
call mma_deallocate(MT_TMP)
call mma_deallocate(ST)
call mma_deallocate(ST_TMP)
call mma_deallocate(MAV)
call mma_deallocate(SAV)
call mma_deallocate(ZT)
call mma_deallocate(ZT_TMP)
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
