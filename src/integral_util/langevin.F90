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

!#define _DEBUGPRINT_
subroutine Langevin(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: DBSC, ncnttp
use Center_Info, only: DC
use Langevin_arrays, only: CAVxyz, DAVxyz, DField, DIP, DipEF, Field, Grid, PolEF, RAVxyz
use External_Centers, only: iXPOLType, nOrd_XF, nXF, nXMolNr, XEle, XF, XMolNr
use Symmetry_Info, only: nIrrep
use rctfld_module, only: CordSI, DipCutOff, DIPSI, DistSparse, Done_Lattice, lDamping, lDipRestart, lFirstIter, lGridAverage, &
                         lLangevin, nCavxyz, nGrid, nGrid_Eff, nGridAverage, nGridSeed, nSparse, PolSI, RadLat, RotAlpha, RotBeta, &
                         RotGamma, Scala, TK
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, auTokJ, kBoltzmann
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nh1
real(kind=wp), intent(inout) :: h1(nh1), TwoHam(nh1), RepNuc
real(kind=wp), intent(in) :: D(nh1)
logical(kind=iwp), intent(in) :: First, Dff, NonEq
integer(kind=iwp) :: mdc, ndc, iCnttp, jCnttp, MaxAto, iCnt, jCnt, nC, mCnt, i, nAV, iSeed, iAV, Lu, j, Inc, iOrdOp, iXF, nCnt, &
                     iEle, k
integer(kind=iwp), save :: nAnisopol, nPolComp
logical(kind=iwp) :: Exists
real(kind=wp) :: ATOD, ATRad, SUMREPNUC, SumRepNuc2, XA, YA, Z, ZA
real(kind=wp), allocatable :: Atom_R(:), Chrg(:), Cord(:,:), D1ao(:), pField(:,:), tmpField(:,:)
real(kind=wp), parameter :: auToK = auTokJ/kBoltzmann*1.0e3_wp
real(kind=wp), external :: CovRadT, Random_molcas

#include "macros.fh"
unused_var(NonEq)

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of all atoms
!
! Cord: list of all atoms
! Atom_R: associated effective atomic radius

mdc = 0
MaxAto = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) cycle
  nCnt = dbsc(iCnttp)%nCntr
  do iCnt=1,nCnt
    mdc = mdc+1
    MaxAto = MaxAto+nIrrep/dc(mdc)%nStab
  end do
end do

call mma_allocate(Cord,3,MaxAto,Label='Cord')
call mma_allocate(Chrg,MaxAto,Label='Chrg')
call mma_allocate(Atom_R,MaxAto,Label='Atom_R')

ndc = 0
nc = 0
do jCnttp=1,nCnttp
  if (dbsc(jCnttp)%Aux .or. dbsc(jCnttp)%Frag) cycle
  Z = dbsc(jCnttp)%Charge
  mCnt = dbsc(jCnttp)%nCntr
  if (dbsc(jCnttp)%AtmNr >= 1) then
    !Atod = CovRad(dbsc(jCnttp)%AtmNr)
    Atod = CovRadT(dbsc(jCnttp)%AtmNr)
  else
    Atod = Zero
  end if
  do jCnt=1,mCnt
    ndc = ndc+1
    do i=0,nIrrep/dc(ndc)%nStab-1
      nc = nc+1
      call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),Cord(:,nc))
      Atom_R(nc) = Atod
      Chrg(nc) = Z
      !write(u6,*) 'Z=',Z
    end do
  end do
end do

if (LGridAverage) then
  nAv = nGridAverage
  if (nGridSeed == 0) then
    iSeed = 0
  else if (nGridSeed == -1) then
    iSeed = 0
    ! Use system_clock only for some systems
    !call System_clock(iSeed,j,k)
  else
    iSeed = nGridSeed
  end if
else
  nAv = 1
end if
sumRepNuc = Zero
sumRepNuc2 = Zero

do iAv=1,nAv
  if (LGridAverage) then
    cordsi(1,1) = Random_Molcas(iSeed)
    cordsi(2,1) = Random_Molcas(iSeed)
    cordsi(3,1) = Random_Molcas(iSeed)
    rotAlpha = Random_Molcas(iSeed)*180.0_wp
    rotBeta = Random_Molcas(iSeed)*180.0_wp
    rotGamma = Random_Molcas(iSeed)*180.0_wp
    write(u6,'(a,i4,a,6f10.4)') 'GRID NR',iAv,': ',cordsi(1,1),cordsi(2,1),cordsi(3,1),rotAlpha,rotBeta,rotGamma
    Done_Lattice = .false.
    RepNuc = Zero
  end if

  if (.not. Done_Lattice) then
    lFirstIter = .true.
    Done_Lattice = .true.
    if (iXPolType == 2) then
      nAnisopol = nXF
      nPolComp = 6
    else
      nAnisopol = 0
      nPolComp = 1
    end if

    ! Compute Effective Polarizabilities on the Langevin grid,
    ! flag also if points on the grid are excluded.
    ! This information is computed once!

    ! Grid : The Langevin grid
    ! PolEf: Effective Polarizability on grid
    ! DipEf: Effective Dipole moment on grid

    nGrid_Eff = 0
    ! Both these subroutines can increase nGrid_Eff
    if (iXPolType > 0) call lattXPol(Grid,nGrid,nGrid_Eff,PolEf,DipEf,XF,nXF,nOrd_XF,nPolComp)
    ! Note: Gen_Grid is now a part of the lattcr subroutine
    if (lLangevin) call lattcr(Grid,nGrid,nGrid_Eff,PolEf,DipEf,Cord,maxato,Atom_R,nPolComp,XF,nXF,nOrd_XF,XEle,iXPolType)
    !write(u6,*) 'nGrid,  nGrid_Eff',nGrid,nGrid_Eff

  end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute Electric Field from the Quantum Chemical System on the
  ! Langevin grid.

  ! cavxyz: MM expansion of the total charge of the QM system
  ! ravxyz: scratch
  ! dField: EF on Langevin grid due to QM system

  ! Get the total 1st order AO density matrix

  call mma_allocate(D1ao,nh1,Label='D1ao')
  call Get_dArray_chk('D1ao',D1ao,nh1)

  ! Save field from permanent multipoles for use in ener
  call mma_allocate(pField,4,nGrid_Eff,Label='pField')
  call mma_allocate(tmpField,4,nGrid_Eff,Label='tmpField')

  pField(:,:) = Zero

  call eperm(D1ao,nh1,Ravxyz,Cavxyz,nCavxyz,dField,Grid,nGrid_Eff,Cord,MaxAto,Chrg,pField)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Save system info, to be used by visualisation program

  Lu = 21
  call OpnFl('LANGINFO',Lu,Exists)
  if (.not. Exists) then
    write(Lu,*) nc
    do i=1,nc
      write(Lu,11) int(Chrg(i)),Atom_R(i),(Cord(j,i),j=1,3)
    end do
    write(Lu,*) nXF
    Inc = 3
    do iOrdOp=0,nOrd_XF
      Inc = Inc+nTri_Elem1(iOrdOp)
    end do
    if (iXPolType > 0) Inc = Inc+6
    do iXF=1,nXF
      xa = XF(1,iXF)
      ya = XF(2,iXF)
      za = XF(3,iXF)
      if (XEle(iXF) <= 0) then
        atrad = -real(XEle(iXF),kind=wp)/1000.0_wp
        iele = 0
      else
        iele = XEle(iXF)
        atrad = CovRadT(iele)
      end if
      write(Lu,11) iele,atrad,xa,ya,za
    end do
    write(Lu,*) nGrid_eff,nAnisopol
    do i=1,nGrid_eff
      write(Lu,12) (Grid(j,i),j=1,3),PolEf(:,i),DipEf(i),(dField(j,i),j=1,3),(pField(j,i),j=1,3)
    end do
    write(Lu,*) polsi,dipsi,scala,auToK/tK
    write(Lu,*) (cordsi(k,1),k=1,3)
    write(Lu,*) rotAlpha,rotBeta,rotGamma
    write(Lu,*) radlat,nSparse,distSparse
    write(Lu,*) lDamping,dipCutoff
  end if
  close(Lu)

  tmpField(:,:) = dField(:,1:nGrid_eff)

  if (lDiprestart .or. lFirstIter) then
    Field(:,:) = Zero
    Dip(:,:) = Zero
    Davxyz(:) = Zero
  end if

  ! Subtract the static MM from the previous iteration
  ! from the static MM of this iteration, and save the
  ! untouched static MM of this iteration into Davxyz
  ! for use in the next iteration. Ravxyz is
  ! just a temporary array

  Ravxyz(:) = Cavxyz(:)
  Cavxyz(:) = Cavxyz(:)-Davxyz(:)
  Davxyz(:) = Ravxyz(:)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Equation solver: compute the Langevin dipole moments and the
  !                  counter charge on the boundary of the cavity.

  ! Field : total EF of the Langevin grid
  ! Dip   : dipole momement on the Langevin grid

  call edip(Field,Dip,dField,PolEf,DipEf,Grid,nGrid_Eff,nPolComp,nAnisopol,nXF,iXPolType,nXMolnr,XMolnr)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute contributions to RepNuc, h1, and TwoHam

  call Ener(h1,TwoHam,D,RepNuc,nh1,First,Dff,D1ao,Grid,nGrid_Eff,Dip,Field,DipEf,PolEf,Cord,MaxAto,Chrg,nPolComp,nAnisopol,pField, &
            tmpField)

  ! Subtract the static field from the self-consistent field
  ! This gives the field from the induced dipoles (saved
  ! in Field, to be used in the next iteration if
  ! not DRES has been requested

  Field(:,:) = Field(:,:)-tmpField(:,:)
  lFirstIter = .false.

  call mma_deallocate(pField)
  call mma_deallocate(tmpField)
  if (LGridAverage) then
    write(u6,'(a,i4,a,f18.10)') 'Solvation energy (Grid nr. ',iAv,'):',RepNuc
    sumRepNuc = sumRepNuc+RepNuc
    sumRepNuc2 = sumRepNuc2+RepNuc**2
  end if
end do
if (LGridAverage) write(u6,'(a,f18.10,f18.10)') 'Average solvation energy and stdev: ',sumRepNuc/real(nAv,kind=wp), &
                                                sqrt(sumRepNuc2/real(nAv,kind=wp)-(sumRepNuc/real(nAv,kind=wp))**2)
call mma_deallocate(Atom_R)
call mma_deallocate(D1ao)
call mma_deallocate(Chrg)
call mma_deallocate(Cord)
!                                                                      *
!***********************************************************************
!                                                                      *
return

11 format(i3,f10.4,3f16.8)
12 format(11f20.10)

end subroutine Langevin
