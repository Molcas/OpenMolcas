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

subroutine magnetization_pa(exch,nLoc,nM,nH,nneq,neq,neqv,nCenter,nTempMagn,nDir,nDirZee,nDirTot,nss,nexch,iopt,LUZee,nsymm,ngrid, &
                            TempMagn,hexp,mexp,hmin,hmax,em,zJ,thrs,dirX,dirY,dirZ,dir_weight,w,dipexch,s_exch,dipso,s_so,eso, &
                            hinput,r_rot,XLM,ZLM,XRM,ZRM,zeeman_energy,compute_Mdir_vector,m_paranoid,m_accurate,smagn,mem,doplot)
! constants defining the sizes
!  NM : number of states included in the exchange Zeeman matrix, ( Nex <= exch)
!  exch, nLoc, nH, nCenter, nTempMagn, nDir, nDirZee, nDirTot, nneq, neqv, neq, nss, nexch, iopt, mem, LUZee, hinput,
!  zeeman_energy, compute_Mdir_vector, m_paranoid, m_accurate, smagn, doplot, R_ROT, W
! exchange energies printed out in the previous part
!  ESO
! spin-orbit energies from ANISO files
!  Hexp, Mexp, thrs, XLM, ZLM, XRM, ZRM, dirX, dirY, dirZ, dir_weight, TempMagn, zJ, hmin, hmax, em, DIPEXCH, S_EXCH, dipso, s_so
!
! exchange data:
!  Wex : Zeeman exchange energies
!  Zex : exchange statistical sum, Boltzmann distribution
!  Sex : spin magnetisation, from the exchange block;
!  Mex : magnetisation, from the exchange block
! data for individual sites (all states):
!  ZL : local statistical sum, Boltzmann distribution
!  WL : Zeeman local energies
!  SL : spin magnetisation, from the local sites, using ALL states ;
!  ML : magnetisation, from local sites, using ALL states;
! data for individual sites (only states that enter exchange):
!  ZR : local statistical sum, Boltzmann distribution, using only Nexch states
!  WR : Zeeman local reduced energies, using only Nexch states;
!  SR : spin magnetisation, from the local sites, using only Nexch states ;
!  MR : magnetisation, from local sites, using only Nexch states;
!  ZRT, ZLT, MRT, MLT, SRT, SLT
! data for total system:
!  ZT : total statistical sum, Boltzmann distribution
!  ST : total spin magnetisation,
!  MT : total magnetisation
! magnetic field strength and orientation data:
!  dltH, H, dHX, dHY, dHZ, dHW
! total average M and average S data:
!  MAV, SAV, MVEC, SVEC

use Lebedev_quadrature, only: order_table
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Ten, mBohr, rNAVO
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: exch, nLoc, nM, nH, nneq, neq(nneq), neqv, nCenter, nTempMagn, nDir, nDirZee, nDirTot, nss(nneq), &
                                 nexch(nneq), iopt, LUZee(nDirZee), nsymm, ngrid, mem
real(kind=wp), intent(in) :: TempMagn(nTempMagn), Hexp(nH), Mexp(nH,nTempMagn), hmin, hmax, em, zJ, thrs, dirX(nDir), dirY(nDir), &
                             dirZ(nDir), dir_weight(nDirZee,3), W(exch), ESO(nneq,nLoc), R_ROT(nneq,neqv,3,3), &
                             XLM(nCenter,nTempMagn,3,3), ZLM(nCenter,nTempMagn), XRM(nCenter,nTempMagn,3,3), ZRM(nCenter,nTempMagn)
complex(kind=wp), intent(in) :: DIPEXCH(3,exch,exch), S_EXCH(3,exch,exch), dipso(nneq,3,nLoc,nLoc), s_so(nneq,3,nLoc,nLoc)
logical(kind=iwp), intent(in) :: hinput, zeeman_energy, compute_Mdir_vector, m_paranoid, m_accurate, smagn, doplot
integer(kind=iwp) :: I, ibuf, iDir, iH, IM, isite, it, itEnd, J, k, mem_local, n, nP
real(kind=wp) :: dltH
character(len=15) :: lbl_X, lbl_Y, lbl_Z
real(kind=wp), allocatable :: dHW(:), dHX(:), dHY(:), dHZ(:), ESO_TMP(:), H(:), MAV(:,:), Mex(:,:), ML(:,:,:), MLT(:,:,:), &
                              MR(:,:,:), MRT(:,:,:), MT(:,:,:), MVEC(:,:,:,:), SAV(:,:), Sex(:,:), SL(:,:,:), SLT(:,:,:), &
                              SR(:,:,:), SRT(:,:,:), ST(:,:,:), SVEC(:,:,:,:), Wex(:), WL(:,:), WR(:,:), Zex(:), ZL(:,:), &
                              ZLT(:,:), ZR(:,:), ZRT(:,:), ZT(:,:)
complex(kind=wp), allocatable :: dipso_tmp(:,:,:), s_so_tmp(:,:,:)
real(kind=wp), parameter :: cm3tomB = rNAVO*mBohr/Ten ! in cm3 * mol-1 * T
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
integer(kind=iwp) :: l
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: DBG = _DBG_
real(kind=wp), external :: dev, dnrm2_

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(mem)
#endif

write(u6,*)
write(u6,'(A)') repeat('%',96)
write(u6,'(40X,A)') 'CALCULATION OF THE MOLAR MAGNETIZATION'
write(u6,'(A)') repeat('%',96)
write(u6,*)
!----------------------------------------------------------------------
mem_local = 0
#ifdef _DEBUGPRINT_
write(u6,*) 'MAGN:        nM=',nM
write(u6,*) 'MAGN:      exch=',exch
write(u6,*) 'MAGN:      nLoc=',nLoc
write(u6,*) 'MAGN: nTempMagn=',nTempMagn
#endif

! Zeeman exchange energy spectrum
call mma_allocate(Wex,nM,'Wex')
mem_local = mem_local+size(Wex)*RtoB

! exchange statistical sum, Boltzmann distribution
call mma_allocate(Zex,nTempMagn,'Zex')
mem_local = mem_local+size(Zex)*RtoB
! spin magnetisation, from the exchange block
call mma_allocate(Sex,3,nTempMagn,'Sex')
mem_local = mem_local+size(Sex)*RtoB
! magnetisation, from the exchange block
call mma_allocate(Mex,3,nTempMagn,'Mex')
mem_local = mem_local+size(Mex)*RtoB

! local statistical sum, Boltzmann distribution
call mma_allocate(ZL,nTempMagn,nneq,'ZL')
mem_local = mem_local+size(ZL)*RtoB
! spin magnetisation, from the local sites, using ALL states
call mma_allocate(SL,3,nTempMagn,nneq,'SL')
mem_local = mem_local+size(SL)*RtoB
! magnetisation, from local sites, using ALL states
call mma_allocate(ML,3,nTempMagn,nneq,'ML')
mem_local = mem_local+size(ML)*RtoB

! local statistical sum, Boltzmann distribution, using only Nexch states
call mma_allocate(ZR,nTempMagn,nneq,'ZR')
mem_local = mem_local+size(ZR)*RtoB
! spin magnetisation, from the local sites, using only Nexch states
call mma_allocate(SR,3,nTempMagn,nneq,'SR')
mem_local = mem_local+size(SR)*RtoB
! magnetisation, from local sites, using only Nexch states
call mma_allocate(MR,3,nTempMagn,nneq,'MR')
mem_local = mem_local+size(MR)*RtoB

! Zeeman local energies
call mma_allocate(WL,nLoc,nneq,'WL')
mem_local = mem_local+size(WL)*RtoB
! Zeeman local reduced energies, using only Nexch states
call mma_allocate(WR,nLoc,nneq,'WR')
mem_local = mem_local+size(WR)*RtoB

! ZRT(nCenter,nTempMagn)
call mma_allocate(ZRT,nCenter,nTempMagn,'ZRT')
mem_local = mem_local+size(ZRT)*RtoB
! ZLT(nCenter,nTempMagn)
call mma_allocate(ZLT,nCenter,nTempMagn,'ZLT')
mem_local = mem_local+size(ZLT)*RtoB
! MRT(nCenter,3,nTempMagn)
call mma_allocate(MRT,nCenter,3,nTempMagn,'MRT')
mem_local = mem_local+size(MRT)*RtoB
! MLT(nCenter,3,nTempMagn)
call mma_allocate(MLT,nCenter,3,nTempMagn,'MLT')
mem_local = mem_local+size(MLT)*RtoB
! SRT(nCenter,3,nTempMagn)
call mma_allocate(SRT,nCenter,3,nTempMagn,'SRT')
mem_local = mem_local+size(SRT)*RtoB
! SLT(nCenter,3,nTempMagn)
call mma_allocate(SLT,nCenter,3,nTempMagn,'SLT')
mem_local = mem_local+size(SLT)*RtoB

! total statistical sum, Boltzmann distribution
call mma_allocate(ZT,nH,nTempMagn,'ZT')
mem_local = mem_local+size(ZT)*RtoB
! total spin magnetisation
call mma_allocate(ST,3,nH,nTempMagn,'ST')
ST(:,:,:) = Zero
mem_local = mem_local+size(ST)*RtoB
! total magnetisation
call mma_allocate(MT,3,nH,nTempMagn,'MT')
mem_local = mem_local+size(MT)*RtoB
! total spin magnetisation
call mma_allocate(SAV,nH,nTempMagn,'SAV')
SAV(:,:) = Zero
mem_local = mem_local+size(SAV)*RtoB
! total magnetisation
call mma_allocate(MAV,nH,nTempMagn,'MAV')
MAV(:,:) = Zero
mem_local = mem_local+size(MAV)*RtoB
! total spin magnetisation vector
call mma_allocate(SVEC,nDirTot,nTempMagn,nH,3,'SVEC')
mem_local = mem_local+size(SVEC)*RtoB
! total magnetisation vector
call mma_allocate(MVEC,nDirTot,nTempMagn,nH,3,'MVEC')
mem_local = mem_local+size(MVEC)*RtoB

! orientation of the field
call mma_allocate(dHX,nDirTot,'dHX')
call mma_allocate(dHY,nDirTot,'dHY')
call mma_allocate(dHZ,nDirTot,'dHZ')
call mma_allocate(dHW,nDirTot,'dHW')
mem_local = mem_local+size(dHX)*RtoB
mem_local = mem_local+size(dHY)*RtoB
mem_local = mem_local+size(dHZ)*RtoB
mem_local = mem_local+size(dHW)*RtoB

call mma_allocate(H,nH,'H  field')
mem_local = mem_local+size(H)*RtoB
#ifdef _DEBUGPRINT_
write(u6,*) 'MAGN:  memory allocated (local):'
write(u6,*) 'mem_local=',mem_local
write(u6,*) 'MAGN:  memory allocated (total):'
write(u6,*) 'mem_total=',mem+mem_local
#endif
!-----------------------------------------------------------------------
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
nP = order_table(nsymm,ngrid)

call hdir(nDir,nDirZee,dirX,dirY,dirZ,dir_weight,nP,nsymm,ngrid,nDirTot,dHX,dHY,dHZ,dHW)
call Add_Info('MR_MAGN  dHX',dHX,min(5,nDirTot),10)
call Add_Info('MR_MAGN  dHY',dHY,min(5,nDirTot),10)
call Add_Info('MR_MAGN  dHZ',dHZ,min(5,nDirTot),10)
call Add_Info('MR_MAGN  dWX',dHW,min(5,nDirTot),10)
#ifdef _DEBUGPRINT_
write(u6,'(A,F10.6)') '        zJ = ',zJ
write(u6,'(A,F10.6)') '      thrs = ',thrs
write(u6,'(A,   I6)') '      iopt = ',iopt
write(u6,'(A, 10I6)') '     nss(i)=',(nss(i),i=1,nneq)
write(u6,'(A, 10I6)') '   nexch(i)=',(nexch(i),i=1,nneq)
write(u6,'(A,   I6)') '  nTempMagn= ',nTempMagn
write(u6,'(A      )') 'TempMagn(i)='
write(u6,'(10F10.6)') (TempMagn(i),i=1,nTempMagn)
write(u6,*) 'm_paranoid =',m_paranoid
write(u6,*) 'm_accurate =',m_accurate
write(u6,*) '     smagn =',smagn
do k=1,nneq
  write(u6,'(A,i6)') 'site ',k
  write(u6,'(A,i2,A)') 'ESO(',k,')'
  write(u6,'(10F12.6)') (ESO(k,i),i=1,nss(k))
  write(u6,'(A,i2,A)') 'S_SO(',k,')'
  do i=1,nss(k)
    do j=1,nss(k)
      write(u6,'(3(2F15.10,3x))') (s_so(k,l,i,j),l=1,3)
    end do
  end do
  write(u6,'(A,i2,A)') 'DIPSO(',k,')'
  do i=1,nss(k)
    do j=1,nss(k)
      write(u6,'(3(2F15.10,3x))') (dipso(k,l,i,j),l=1,3)
    end do
  end do
end do !k
#endif

write(u6,'(2X,A,i3,A)') 'Molar magnetization will be calculated in ',nH,' points, equally distributed in magnetic field range'
write(u6,'(2X,F4.1,1x,a,1X,F4.1,a,5(F6.3,a))') HMIN,'--',HMAX,' T., at the following temperatures:'
do i=1,nTempMagn,10
  j = min(nTempMagn,i+9)
  write(u6,'(10(F8.4,A))') (TempMagn(k),' K.;',k=i,j)
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
if (m_accurate) then
  write(u6,'(2x,A)') 'The contribution of local excited states  is computed exactly.'
  if ((zJ /= Zero) .and. m_paranoid .and. (nTempMagn > 1)) then
    write(u6,'(2x,A)') 'The average spin is computed exactly for each temperature point.'
  else if ((zJ /= Zero) .and. (.not. m_paranoid) .and. (nTempMagn > 1)) then
    write(u6,'(2x,A,F9.3)') 'The average spin is computed exactly only for the temperature point T= ',maxval(TempMagn(:))
    write(u6,'(2x,A     )') 'We consider this to be a good approximation.'
  end if
else
  if ((zJ /= Zero) .and. m_paranoid .and. (nTempMagn > 1)) then
    write(u6,'(2x,A)') 'The average spin is computed exactly for each temperature point.'
  else if ((zJ /= Zero) .and. (.not. m_paranoid) .and. (nTempMagn > 1)) then
    write(u6,'(2x,A,F9.3)') 'The average spin is computed exactly only for the temperature  point T= ',TempMagn(1)
    write(u6,'(2x,A)') 'We consider this to be a good approximation.'
    write(u6,'(2x,A)') 'Use MPAR keyword to compute average spin exactly, for each temperature point, in case higher accuracy '// &
                       'is needed.'
  end if
  write(u6,'(2x,A)') 'The contribution of local excited states  is computed approximately (by using the susceptibility data). '
end if
if (compute_Mdir_vector) then
  write(u6,'(2x,A,i2,a)') 'The magnetization vector for ',nDir,' directions of the applied magnetic field will be calculated.'
else
  write(u6,'(2X,A)') 'The magnetization vector was not calculated.'
end if
if (zeeman_energy) then
  write(u6,'(2x,A,i2,a)') 'The Zeeman splitting for ',nDirZee,' directions of the applied magnetic field will be calculated.'
  write(u6,'(2x,a)') 'The Zeeman energies for each direction of the applied magnetic field are written in files "energy_XX.txt".'
# ifdef _DEBUGPRINT_
  do i=1,nDirZee
    write(u6,'(A,I4,A,I4)') 'LuZee(',i,' )=',LUZee(i)
  end do
# endif
else
  write(u6,'(2X,A)') 'Computation of the Zeeman splitting was not requested.'
end if
!if (maxval(TempMagn) < W(exch)) then
!  write(u6,'(2x,A)') 'Contribution to molar magnetization coming from local excited states is taken into account.'
!else
!  write(u6,'(2x,A)') 'Contribution to molar magnetization coming from local excited states is NOT taken into account.'
!  write(u6,'(2x,A)') 'TMAG is requesting to compute magnetization at a larger temperature than the highest exchange state.'
!  write(u6,'(2x,A)') 'Please include more states into the exchange coupling.'
!end if

! /// opening the loop over the field points
do iH=1,nH
  ! /// ----------------------------------------------------------------
  if (HINPUT) then
    H(iH) = HEXP(iH)
    if (H(iH) == Zero) H(iH) = 0.0001_wp
  else
    DLTH = (HMAX-HMIN)/real(nH-1,kind=wp)
    if (iH == 1) then
      H(iH) = HMIN+0.0001_wp
    else
      H(iH) = HMIN+DLTH*real(iH-1,kind=wp)
    end if
    if (H(iH) == Zero) H(iH) = 0.0001_wp
  end if

# ifdef _DEBUGPRINT_
  write(u6,'(A,i0,A,F10.5)') 'MAGNETIZATION::  H(',iH,') = ',H(iH)
# endif

  ! ///  opening the loop over different directions of the magnetic field
  do IM=1,NDIRTOT
    ZL(:,:) = Zero
    ZR(:,:) = Zero
    SL(:,:,:) = Zero
    ML(:,:,:) = Zero
    SR(:,:,:) = Zero
    MR(:,:,:) = Zero
#   ifdef _DEBUGPRINT_
    write(u6,'(A,F20.10)') 'zJ = ',zJ
#   endif
    ! exchange magnetization:
    call MAGN(exch,NM,dHX(iM),dHY(iM),dHZ(iM),H(iH),W,zJ,THRS,DIPEXCH,S_EXCH,nTempMagn,TempMagn,smagn,Wex,Zex,Sex,Mex,m_paranoid, &
              DBG)
#   ifdef _DEBUGPRINT_
    write(u6,'(A,3F11.7)') 'MEX:',(Mex(l,1),l=1,3)
#   endif
    if (IM == 1) then
      call Add_Info('MEX_MAGN    ',[dnrm2_(3*nTempMagn,Mex,1)],1,8)
      call Add_Info('MR_MAGN  Wex',[dnrm2_(nM,Wex,1)],1,8)
    end if
    ! compute local magnetizations:
    if (m_accurate) then
      do i=1,nneq
        ! all states:
        if (NSS(i) > NEXCH(i)) then
          call mma_allocate(ESO_TMP,nss(i),label='ESO_TMP')
          call mma_allocate(dipso_tmp,3,nss(i),nss(i),label='dipso_tmp')
          call mma_allocate(s_so_tmp,3,nss(i),nss(i),label='s_so_tmp')
          ESO_TMP(:) = ESO(i,1:NSS(i))
          dipso_tmp(:,:,:) = DIPSO(i,:,1:NSS(i),1:NSS(i))
          s_so_tmp(:,:,:) = S_SO(i,:,1:NSS(i),1:NSS(i))
          ! this check is to avoid the unnecessary computation, in cases when no local excited states are present
          call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),H(iH),ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempMagn,TempMagn,smagn, &
                    WL(1:NEXCH(i),i),ZL(:,i),SL(:,:,i),ML(:,:,i),m_paranoid,DBG)
          call mma_deallocate(ESO_TMP)
          call mma_deallocate(dipso_tmp)
          call mma_deallocate(s_so_tmp)
          if (IM == 2) call Add_Info('MR_MAGN  WL',[dnrm2_(nexch(i),WL(:,i),1)],1,8)
#         ifdef _DEBUGPRINT_
          write(u6,'(A,I2,A,3F11.7)') 'ML: site',i,' : ',(ML(l,1,i),l=1,3)
#         endif
          call mma_allocate(ESO_TMP,nexch(i),label='ESO_TMP')
          call mma_allocate(dipso_tmp,3,nexch(i),nexch(i),label='dipso_tmp')
          call mma_allocate(s_so_tmp,3,nexch(i),nexch(i),label='s_so_tmp')
          ESO_TMP(:) = ESO(i,1:NEXCH(i))
          dipso_tmp(:,:,:) = DIPSO(i,:,1:NEXCH(i),1:NEXCH(i))
          s_so_tmp(:,:,:) = S_SO(i,:,1:NEXCH(i),1:NEXCH(i))
          ! only local "exchange states":
          call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),H(iH),ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempMagn,TempMagn,smagn, &
                    WR(1:NEXCH(i),i),ZR(:,i),SR(:,:,i),MR(:,:,i),m_paranoid,DBG)
          call mma_deallocate(ESO_TMP)
          call mma_deallocate(dipso_tmp)
          call mma_deallocate(s_so_tmp)
          if (IM == 2) call Add_Info('MR_MAGN  WR',[dnrm2_(nexch(i),WR(:,i),1)],1,8)
#         ifdef _DEBUGPRINT_
          write(u6,'(A,I2,A,3F11.7)') 'MR: site',i,' : ',(MR(l,1,i),l=1,3)
#         endif
        end if
      end do

      if (IM == 3) then
        call Add_Info('ML_MAGN',[dnrm2_(3*nTempMagn*nneq,ML,1)],1,8)
        call Add_Info('MR_MAGN',[dnrm2_(3*nTempMagn*nneq,MR,1)],1,8)
      end if

      ! expand the basis and rotate local vectors to the general
      ! coordinate system:
      MRT(:,:,:) = Zero
      MLT(:,:,:) = Zero
      SRT(:,:,:) = Zero
      SLT(:,:,:) = Zero
      isite = 0
      do i=1,NNEQ
        do j=1,NEQ(i)
          isite = isite+1
          ! statistical distributions
          ZLT(isite,:) = ZL(:,i)
          ZRT(isite,:) = ZR(:,i)
          ! magnetizations:
          !  use R_rot matrices, which have determinant +1.
          !  note that  R_lg matrices have arbitrary determinant.
          do iT=1,nTempMagn
            do n=1,3
              MLT(isite,:,iT) = MLT(isite,:,iT)+r_rot(i,j,:,n)*ML(n,iT,i)
              SLT(isite,:,iT) = SLT(isite,:,iT)+r_rot(i,j,:,n)*SL(n,iT,i)
              MRT(isite,:,iT) = MRT(isite,:,iT)+r_rot(i,j,:,n)*MR(n,iT,i)
              SRT(isite,:,iT) = SRT(isite,:,iT)+r_rot(i,j,:,n)*SR(n,iT,i)
            end do
          end do
        end do ! j, neq(i)
      end do ! i, nneq
      ZRT(isite:,:) = Zero
      ZLT(isite:,:) = Zero
    end if ! m_accurate

    ! compute the total magnetizations according to the derived formulas:
    if (m_accurate) then
      do iT=1,nTempMagn
        if (smagn) call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLT(:,iT),SRT(:,:,iT),ZRT(:,iT),iopt,ST(:,iH,iT),ZT(iH,iT))
        call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLT(:,iT),MRT(:,:,iT),ZRT(:,iT),iopt,MT(:,iH,iT),ZT(iH,iT))
      end do
    else
      ! add the contribution from local excited states using the approximate X*H expression:
      do iT=1,nTempMagn
        do isite=1,nCenter
          MRT(isite,:,iT) = H(iH)*(XRM(isite,iT,:,1)*dHX(iM)+XRM(isite,iT,:,2)*dHY(iM)+XRM(isite,iT,:,3)*dHZ(iM))/cm3tomB
          MLT(isite,:,iT) = H(iH)*(XLM(isite,iT,:,1)*dHX(iM)+XLM(isite,iT,:,2)*dHY(iM)+XLM(isite,iT,:,3)*dHZ(iM))/cm3tomB
        end do
        if (smagn) then
          SRT(:,:,:) = Zero
          SLT(:,:,:) = Zero
          call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLM(:,iT),SRT(:,:,iT),ZRM(:,iT),iopt,ST(:,iH,iT),ZT(iH,iT))
        end if
        call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLM(:,iT),MRT(:,:,iT),ZRM(:,iT),iopt,MT(:,iH,iT),ZT(iH,iT))
      end do
    end if
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! print out hte Zeeman eigenstates
    if (zeeman_energy) then
      if ((iH == 1) .and. (iM == nDir+1)) &
        write(u6,'(A)') 'Energies of the Zeeman Hamiltonian for the following directions of the applied field:'
      if ((iH == 1) .and. (iM > nDir) .and. (iM <= nDir+nDirZee)) then
        write(u6,'(A,I3,A,3F10.6,3x,5A)') 'direction Nr.',iM-nDir,' : ',dHX(iM),dHY(iM),dHZ(iM),'written in file "zeeman_energy_', &
                                          char(48+mod((iM-nDir)/100,10)),char(48+mod((iM-nDir)/10,10)),char(48+mod(iM-nDir,10)), &
                                          '.txt".'
        write(LUZee(iM-nDir),'(A,3F24.15)') '# direction of the applied magnetic field:',dHX(iM),dHY(iM),dHZ(iM)
        write(LuZee(iM-nDir),'(A,6x,A,1000(I4,6x) )') '# H(T)',' State =>',(i,i=1,nm)
      end if

      if ((iM > nDir) .and. (iM <= nDir+nDirZee)) write(LUZee(iM-nDir),'(F8.4,1000F10.3)') H(IH),(Wex(I),I=1,NM)
    end if !zeeman_energy
    ! ------------------------------------------------------------------
    ! computing the AVERAGE MOMENTS calculated at different temperatures
    ! (TempMagn(i))
    do iT=1,nTempMagn
      MVEC(iM,iT,iH,:) = MT(:,iH,iT)
      SVEC(iM,iT,iH,:) = ST(:,iH,iT)
    end do !iT

    ! accumulate contributions:
    MAV(iH,:) = MAV(iH,:)+dHW(iM)*(dHX(iM)*MT(1,iH,:)+dHY(iM)*MT(2,iH,:)+dHZ(iM)*MT(3,iH,:))
    SAV(iH,:) = SAV(iH,:)+dHW(iM)*(dHX(iM)*ST(1,iH,:)+dHY(iM)*ST(2,iH,:)+dHZ(iM)*ST(3,iH,:))
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
    !write(u6,*)
    do iDir=1,nDir+nDirZee
      !write(u6,*)
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
        write(u6,'(F7.3,1x,A,3(ES13.6,1x,A),ES14.7,1x,A,1x,A,3(ES13.6,1x,A),ES14.7,1x,A)') &
          H(iH),'|',MVEC(iDir,iT,iH,1),' ',MVEC(iDir,iT,iH,2),' ',MVEC(iDir,iT,iH,3),'|', &
          (MVEC(iDir,iT,iH,1)*dHX(iDir)+MVEC(iDir,iT,iH,2)*dHY(iDir)+MVEC(iDir,iT,iH,3)*dHZ(iDir)),'|','|',SVEC(iDir,iT,iH,1),' ', &
          SVEC(iDir,iT,iH,2),' ',SVEC(iDir,iT,iH,3),'|', &
          (SVEC(iDir,iT,iH,1)*dHX(iDir)+SVEC(iDir,iT,iH,2)*dHY(iDir)+SVEC(iDir,iT,iH,3)*dHZ(iDir)),'|'
      end do
      write(u6,'(A,A,1x,A)') '--------|','------------------------------------------------------------|', &
                             '|------------------------------------------------------------|'
    end do !iDir
  end do !iT

else

  do iT=1,nTempMagn
    do iDir=1,nDir+nDirZee
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
        write(u6,'(F7.3,1x,A,3(ES13.6,1x,A),ES14.7,1x,A)') &
          H(iH),'|',MVEC(iDir,iT,iH,1),' ',MVEC(iDir,iT,iH,2),' ',MVEC(iDir,iT,iH,3),'|', &
          (MVEC(iDir,iT,iH,1)*dHX(iDir)+MVEC(iDir,iT,iH,2)*dHY(iDir)+MVEC(iDir,iT,iH,3)*dHZ(iDir)),'|'
      end do
      write(u6,'(2A)') '--------|','------------------------------------------------------------|'
    end do !iDir

  end do !iT
end if !(smagn)
! ---------------------------------------------------------------------
! COMPUTING THE STANDARD DEVIATION OF THE MAGNETIZATION
!if (HINPUT) then
!  do iT=1,nTempMagn
!    STDEV(iT) = dev(nH,MAV(:,iT),Mexp(:,iT))
!  end do
!end if

write(u6,*)
write(u6,'(15X,A)') 'HIGH-FIELD POWDER MAGNETIZATION'
write(u6,'(20X,A)') '(Units: Bohr magneton)'
write(u6,*)
do iT=1,nTempMagn,5
  iTEnd = min(nTempMagn,iT+4)
  write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  write(u6,'(A,11(F10.3,A))') '    H(T)   |STATISTICAL SUM|',(TempMagn(i),' K.  |',i=iT,iTEnd)
  write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  do iH=1,nH
    write(u6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))') H(iH),'|',ZT(iH,1),'|',(MAV(iH,i),'|',i=iT,iTEnd)
  end do
  write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
  if (HINPUT) then
    write(u6,'(A,15x,A, 11(f14.10,1x,A) )') 'ST.DEV.M   |','|',(dev(nH,MAV(:,i),Mexp(:,i)),'|',i=iT,iTEnd)
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
    write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
    write(u6,'(A,11(F10.3,A))') '   H(T)   |STATISTICAL SUM|',(TempMagn(i),' K.  |',i=iT,iTEnd)
    write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
    do iH=1,nH
      write(u6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))') H(iH),'|',ZT(iH,1),'|',(SAV(iH,i),'|',i=iT,iTEnd)
    end do
    write(u6,'(A,11A)') '-----------|',('---------------|',i=iT,iTEnd+1)
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
  call Add_Info(lbl_X,[dnrm2_(ibuf,MVEC(:,:,iH,1),1)],1,8)
  call Add_Info(lbl_Y,[dnrm2_(ibuf,MVEC(:,:,iH,2),1)],1,8)
  call Add_Info(lbl_Z,[dnrm2_(ibuf,MVEC(:,:,iH,3),1)],1,8)
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
