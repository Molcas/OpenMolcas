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
subroutine DrvPCM(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq)

use Basis_Info, only: DBSC, nBas, nCnttp
use Center_Info, only: DC
use PCM_arrays, only: C_Tessera, DiagScale, nTiles, PCMDM, PCMTess, Q_Tessera
use Symmetry_Info, only: iChBas, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use rctfld_module, only: Eps, EpsInf, nTs
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nh1
real(kind=wp), intent(inout) :: h1(nh1), TwoHam(nh1), RepNuc
real(kind=wp), intent(in) :: D(nh1)
logical(kind=iwp), intent(in) :: First, Dff, NonEq
integer(kind=iwp) :: i, jCnt, jCnttp, MaxAto, mCnt, nC, nDC
real(kind=wp) :: Z
real(kind=wp), allocatable :: Chrg(:), Cord(:,:), PCM_charge(:,:), Q_Slow(:), V_Save(:,:), V_Slow(:), V_Tile(:,:)

#include "macros.fh"
unused_var(Dff)

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate list of all atoms

! Cord: list of all atoms

call Get_nAtoms_All(MaxAto)

call mma_allocate(Cord,3,MaxAto,Label='Cord')
call mma_allocate(Chrg,MaxAto,Label='Chrg')

ndc = 0
nc = 1
do jCnttp=1,nCnttp
  Z = dbsc(jCnttp)%Charge
  mCnt = dbsc(jCnttp)%nCntr
  if (dbsc(jCnttp)%Aux) mCnt = 0
  do jCnt=1,mCnt
    ndc = ndc+1
    do i=0,nIrrep/dc(ndc)%nStab-1
      call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),Cord(1:3,nc))
      Chrg(nc) = Z
      nc = nc+1
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(PCM_Charge,2,nTS,Label='PCM_Charge')
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(V_Tile,2,nTs,Label='V_Tile')
call mma_allocate(V_Save,2,nTs,Label='V_Save')
call mma_allocate(Q_Slow,nTs,Label='Q_Slow')
call mma_allocate(V_Slow,nTs,Label='V_Slow')

call DrvPCM_Internal(Chrg,Cord,MaxAto,PCMTess,PCMDM,V_Tile,V_Save,PCM_Charge,Q_Slow,V_Slow,nTs,Eps,EpsInf)

call mma_deallocate(V_Slow)
call mma_deallocate(Q_Slow)
call mma_deallocate(V_Save)
call mma_deallocate(V_Tile)
!                                                                      *
!***********************************************************************
!                                                                      *
! Put the current set of PCM charges on the run file.

call Put_dArray('PCM Charges',PCM_Charge,2*nTs)
call mma_deallocate(PCM_Charge)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Chrg)
call mma_deallocate(Cord)
!                                                                      *
!***********************************************************************
!                                                                      *
return

contains

subroutine DrvPCM_Internal(Z_Nuc,Cord,MaxAto,Tessera,DMat,VTessera,VSave,QTessera,QTessera_Slow,VSlow,nTs,Eps,EpsInf)

  use Gateway_global, only: PrPrt
  use Integral_Interfaces, only: int_kernel, int_mem, OneEl_Integrals
  use Constants, only: Zero, One, Two, Half, Pi

  implicit none
  integer(kind=iwp), intent(in) :: MaxAto, nTs
  real(kind=wp), intent(in) :: Z_Nuc(MaxAto), Cord(3,MaxAto), Tessera(4,nTs), Eps, EpsInf
  real(kind=wp), intent(inout) :: DMat(nTs,nTs)
  real(kind=wp), intent(out) :: VTessera(2,nTs), VSave(2,nTs), QTessera(2,nTs), QTessera_Slow(nTs), VSlow(nTs)
  integer(kind=iwp) :: ip(1), iTile, jTile, kOper, lOper, n_Int, nComp, nOrdOp
  real(kind=wp) :: Alpha, Dij, EEE, EEN, ENE, ENN, Fact, Origin(3), QInf, rHrmt, Rij, W_0_or_El, W_0_or_Inf, W_or_El, W_or_Inf, &
                   W_or_InfNuc, W_or_Nuc, Xi, Xj, Yi, Yj, Zi, Zj
  logical(kind=iwp) :: Save_tmp
  character(len=8) :: Label
  integer(kind=iwp), allocatable :: lOper2(:)
  real(kind=wp), allocatable :: FactOp(:), Integrals(:)
  procedure(int_kernel) :: PCMInt
  procedure(int_mem) :: NaMem
  integer(kind=iwp), external :: n2Tri

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Reaction field a la PCM

  ! Set up some parameters.

  rHrmt = One
  nOrdOp = 0
  nComp = 1
  lOper = 255
  kOper = iChBas(1)
  Origin(:) = Zero
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Pick up slow charges and W's. These are present if we are doing
  ! the final state.

  if (NonEq) then

    ! Read slow components originating from the initial state

    !write(u6,*) 'Rd:',QTessera_Slow(1)
    call Get_dArray('RCTFLD',QTessera_Slow,nTs)
    call Get_dScalar('W_or_el',W_0_or_el)
    call Get_dScalar('W_or_Inf',W_0_or_Inf)

    ! Compute the electrostatic potential due to slow charges

    VSlow(:) = Zero
    do iTile=1,nTs
      XI = Tessera(1,iTile)
      YI = Tessera(2,iTile)
      ZI = Tessera(3,iTile)
      do jTile=1,nTs
        if (jTile == iTile) then
          Dij = DiagScale*Two*sqrt(Pi/Tessera(4,iTile))
        else
          XJ = Tessera(1,jTile)
          YJ = Tessera(2,jTile)
          ZJ = Tessera(3,jTile)
          RIJ = sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
          Dij = One/RIJ
        end if
        VSlow(iTile) = VSlow(iTile)+Dij*QTessera_Slow(jTile)
      end do
    end do
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Evaluate the potential field on the tiles.

  Save_tmp = PrPrt
  PrPrt = .true.

  ! Do the nuclear contribution

  do iTile=1,nTs
    call EFNuc(Tessera(1,iTile),Z_Nuc,Cord,MaxAto,VTessera(1,iTile),nOrdOp)
    VTessera(2,iTile) = Zero
  end do

  ! Do the electronic contribution

  call mma_allocate(FactOp,nTs,Label='FactOp')
  call mma_allocate(lOper2,nTs,Label='lOper2')
  FactOp(:) = One
  lOper2(:) = 255

  call Drv1_PCM(FactOp,nTs,D,nh1,Tessera,lOper2,VTessera,nOrdOp)

  call mma_deallocate(lOper2)
  call mma_deallocate(FactOp)

  ! Save the electrostatic potential

  call dcopy_(2*nTs,VTessera,1,VSave,1)

  ! Add the slow charge contribution to the nuclear potential

  if (NonEq) call DaXpY_(nTs,One,VSlow,1,VTessera(1,1),2)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Evaluate the charges on the cavity boundary, nuclear and
  ! electronic.

  call PCM_Driver(DMat,VTessera,QTessera,nTs)

  ! Make the slow charges (also called orientational charges or
  ! frozen charges). This is always done regardless if they ever
  ! will be used.  The charges can be used in non-equilibrium
  ! calculations.

  if ((EpsInf > Zero) .and. (.not. NonEq)) then
    Fact = (Eps-EpsInf)/(Eps-One)
    call dcopy_(nTs,QTessera(1,1),2,QTessera_Slow,1)
    call DaXpY_(nTs,One,QTessera(2,1),2,QTessera_Slow,1)
    call DScal_(nTs,Fact,QTessera_Slow,1)

    ! Compute the electrostatic potential due to slow charges

    call dcopy_(nTs,[Zero],0,VSlow,1)
    do iTile=1,nTs
      XI = Tessera(1,iTile)
      YI = Tessera(2,iTile)
      ZI = Tessera(3,iTile)
      do jTile=1,nTs
        if (jTile == iTile) then
          Dij = DiagScale*Two*sqrt(Pi/Tessera(4,iTile))
        else
          XJ = Tessera(1,jTile)
          YJ = Tessera(2,jTile)
          ZJ = Tessera(3,jTile)
          RIJ = sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
          Dij = One/RIJ
        end if
        VSlow(iTile) = VSlow(iTile)+Dij*QTessera_Slow(jTile)
      end do
    end do

  end if

  ! Recover the nuclear potential (discarding the possible variations
  ! due to slow charges)

  call dcopy_(nTs,VSave(1,1),2,VTessera(1,1),2)

  W_or_el = Zero
  W_or_nuc = Zero
  W_or_Inf = Zero
  W_or_InfNuc = zero
  do iTile=1,nTs
    W_or_el = W_or_el+QTessera_Slow(iTile)*VTessera(2,iTile)
    W_or_nuc = W_or_nuc+QTessera_Slow(iTile)*VTessera(1,iTile)
    if (NonEq) then
      QInf = QTessera(1,iTile)+QTessera(2,iTile)
    else
      Fact = (Eps-Epsinf)/(Eps-One)
      QInf = (QTessera(1,iTile)+QTessera(2,iTile))*(One-Fact)
    end if
    W_or_Inf = W_or_Inf+QInf*VSlow(iTile)
    W_or_InfNuc = W_or_InfNuc+QTessera(1,iTile)*VSlow(iTile)
  end do

  ! Write out the slow charges and constants if this is an
  ! equilibrium calculations.

  if ((EpsInf > Zero) .and. (.not. NonEq)) then
    !write(u6,*) 'Wr:',QTessera_Slow(1)
    call Put_dArray('RCTFLD',QTessera_Slow,nTs)
    call Put_dScalar('W_or_el',W_or_el)
    call Put_dScalar('W_or_Inf',W_or_Inf)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now add terms to RepNuc

  ! Interaction terms:
  ! nuclei-nuclear solvation charge       (ENN)
  ! nuclei-electronic solvation charge    (ENE)
  ! electrons-nuclear solvation charge    (EEN)
  ! electrons-electronic solvation charge (EEE)

  ENN = Zero
  ENE = Zero
  EEN = Zero
  EEE = Zero
  do iTile=1,nTs
    ENN = ENN+QTessera(1,iTile)*VTessera(1,iTile)
    ENE = ENE+QTessera(1,iTile)*VTessera(2,iTile)
    EEN = EEN+QTessera(2,iTile)*VTessera(1,iTile)
    EEE = EEE+QTessera(2,iTile)*VTessera(2,iTile)
  end do
  if (First) then
    RepNuc = RepNuc+Half*ENN
    if (NonEq) RepNuc = RepNuc+Half*W_or_nuc+Half*W_or_InfNuc-Half*W_0_or_el-Half*W_0_or_Inf
    Label = 'PotNuc00'
    call Put_Temp(Label,[RepNuc],1)
  end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now add terms to h1 and TwoHam!

  call mma_allocate(C_Tessera,3,nTs,Label='C_Tessera')
  nTiles = nTs
  C_Tessera(:,:) = Tessera(1:3,:)
  call mma_allocate(Q_Tessera,nTs,Label='Q_Tessera')
  Label = '<Q|V>'

  if (First) then

    ! PCM-integrals weighted by Q(1)
    ! h1 + correction

    Q_Tessera(:) = QTessera(1,:)
    if (NonEq) call DaXpY_(nTs,One,QTessera_Slow,1,Q_Tessera,1)
    call OneEl_Integrals(PCMInt,NaMem,Label,ip,[lOper],nComp,Origin,nOrdOp,rHrmt,[kOper],Integrals)
    n_Int = n2Tri(lOper)
    call CmpInt(Integrals(ip(1)),n_Int,nBas,nIrrep,lOper)
    Alpha = One
    call DaXpY_(n_Int,Alpha,Integrals(ip(1)),1,h1,1)
    call mma_deallocate(Integrals)

    ! Save the modified h1 matrix

    Label = 'h1_raw  '
    call Put_Temp(Label,h1,nh1)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! PCM-integrals weighted by Q(2)
  ! TwoHam + correction

  Q_Tessera(:) = QTessera(2,:)
  call OneEl_Integrals(PCMInt,NaMem,Label,ip,[lOper],nComp,Origin,nOrdOp,rHrmt,[kOper],Integrals)
  n_Int = n2Tri(lOper)
  call CmpInt(Integrals(ip(1)),n_Int,nBas,nIrrep,lOper)
  Alpha = One
  call DaXpY_(n_Int,Alpha,Integrals(ip(1)),1,TwoHam,1)
  call mma_deallocate(Integrals)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_deallocate(Q_Tessera)
  call mma_deallocate(C_Tessera)
  PrPrt = Save_tmp
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  return

end subroutine DrvPCM_Internal

end subroutine DrvPCM
