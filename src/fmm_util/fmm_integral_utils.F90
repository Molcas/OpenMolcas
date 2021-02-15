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

module fmm_integral_utils

use fmm_global_paras, only: INTK, REALK, fmm_basis, fmm_prim_batch, One, Two, Half
implicit none
private
! Public procedures
public :: fmm_get_prim_batch, &
          fmm_build_Ecoef1, &
          fmm_build_Ecoef2

contains

!-------------------------------------------------------------------------------

subroutine fmm_get_prim_batch(basis,Ish,Jsh,batch,NPrim)

  implicit none
  type(fmm_basis), intent(in)       :: basis
  integer(INTK), intent(in)         :: Ish, Jsh
  type(fmm_prim_batch), intent(out) :: batch(:)
  integer(INTK), intent(out)        :: NPrim

!fixme
  real(REALK), parameter :: ThrFac = -2.30258_REALK*13
  integer(INTK) :: I, J
  integer(INTK) :: IPrim1, JPrim1
  integer(INTK) :: IPrim2, JPrim2, JPTemp2
  real(REALK) :: Acentr(3), Bcentr(3), P(3), RAB(3)
  real(REALK) :: ExpA, ExpB, ExpP, ExpPI, ExpAR2, R2AB, ExpKAB

  Acentr(:) = basis%Centr(:,basis%KAtom(Ish))
  IPrim1 = basis%KStart(Ish)
  IPrim2 = IPrim1+basis%KontG(Ish)-1

  Bcentr(:) = basis%Centr(:,basis%KAtom(Jsh))
  JPrim1 = basis%KStart(Jsh)
  JPrim2 = JPrim1+basis%KontG(Jsh)-1

  RAB(:) = Acentr(:)-Bcentr(:)
  R2AB = dot_product(RAB,RAB)

  NPrim = 0
  do I=IPrim1,IPrim2
    ExpA = basis%Expnt(I)
    ExpAR2 = ExpA*R2AB
    JPTemp2 = JPrim2
    if (Ish == Jsh) JPTemp2 = I
    do J=JPrim1,JPTemp2
      ExpB = basis%Expnt(J)
      ExpP = ExpA+ExpB
      ExpPI = One/ExpP
      ExpKAB = -ExpAR2*ExpB*ExpPI
      if (ExpKAB >= ThrFac) then
        NPrim = NPrim+1
        batch(NPrim)%ExpntP = ExpP
        batch(NPrim)%ExpPHalf = half*ExpPI
        batch(NPrim)%PreFactAB = exp(ExpKAB)
        batch(NPrim)%CCoefAB = basis%CCoef(I)*basis%CCoef(J)
        if ((Ish == Jsh) .and. (I /= J)) then
          batch(NPrim)%CCoefAB = two*batch(NPrim)%CCoefAB
        end if
        P(:) = (ExpA*Acentr(:)+ExpB*Bcentr(:))*ExpPI
        batch(NPrim)%P(:) = P(:)
        batch(NPrim)%PA(:) = P(:)-Acentr(:)
        batch(NPrim)%PB(:) = P(:)-Bcentr(:)
        !print *, 'ExpP:', ExpP
        !print *, 'PreFactAB:', batch(NPrim)%PreFactAB
        !print *, 'CCoefAB:', batch(NPrim)%CCoefAB
        !print *, 'CntrP:', batch(NPrim)%P
      end if
    end do
  end do

end subroutine fmm_get_prim_batch

!-------------------------------------------------------------------------------

subroutine fmm_build_Ecoef1(batch,NPrim,IAngl,JAngl,ECoefX,ECoefY,ECoefZ)

  implicit none
  type(fmm_prim_batch), intent(in) :: batch(:)
  integer(INTK), intent(in)        :: NPrim, IAngl, JAngl
  real(REALK), intent(out)         :: ECoefX(0:,0:,0:,:)
  real(REALK), intent(out)         :: ECoefY(0:,0:,0:,:)
  real(REALK), intent(out)         :: ECoefZ(0:,0:,0:,:)

  integer(INTK) :: I, J, It, IJ, JTemp
  real(REALK) :: PXAX(NPrim), PXBX(NPrim)
  real(REALK) :: PYAY(NPrim), PYBY(NPrim)
  real(REALK) :: PZAZ(NPrim), PZBZ(NPrim)
  real(REALK) :: ExpHalf(NPrim)
  real(REALK) :: PreFact(NPrim)

  ! o Generate E-coefficients (Hermite expansion coefficients) with the three-term recurrence relation

  PXAX(:) = batch(1:Nprim)%PA(1)
  PYAY(:) = batch(1:Nprim)%PA(2)
  PZAZ(:) = batch(1:Nprim)%PA(3)
  PXBX(:) = batch(1:Nprim)%PB(1)
  PYBY(:) = batch(1:Nprim)%PB(2)
  PZBZ(:) = batch(1:Nprim)%PB(3)
  ExpHalf(:) = batch(1:Nprim)%ExpPHalf
  PreFact(:) = batch(1:Nprim)%PreFactAB

  ! o E(0,0,0,NPrim), The prefactor is incorporated.

  ECoefX(0,0,0,1:NPrim) = One
  ECoefY(0,0,0,1:NPrim) = One
  ECoefZ(0,0,0,1:NPrim) = PreFact(1:NPrim)

  ! o Recurrence relation
  !   E(I,J,t,NPrim), 0 <= t <= I + J, IAngl >= JAngl

  if (IAngl == 0) return
  do I=1,IAngl
    JTemp = I
    if (I == IAngl) JTemp = JAngl
    do J=0,JTemp

      ! * Case of It = 0

      It = 0
      if (I+J <= 1) then
        do IJ=1,NPrim
          ECoefX(I,J,It,IJ) = PXAX(IJ)*ECoefX(I-1,J,It,IJ)
          ECoefY(I,J,It,IJ) = PYAY(IJ)*ECoefY(I-1,J,It,IJ)
          ECoefZ(I,J,It,IJ) = PZAZ(IJ)*ECoefZ(I-1,J,It,IJ)
        end do
      else
        do IJ=1,NPrim
          ECoefX(I,J,It,IJ) = PXAX(IJ)*ECoefX(I-1,J,It,IJ)+ECoefX(I-1,J,It+1,IJ)
          ECoefY(I,J,It,IJ) = PYAY(IJ)*ECoefY(I-1,J,It,IJ)+ECoefY(I-1,J,It+1,IJ)
          ECoefZ(I,J,It,IJ) = PZAZ(IJ)*ECoefZ(I-1,J,It,IJ)+ECoefZ(I-1,J,It+1,IJ)
        end do
      end if

      ! * Case of It >= 1

      do It=1,I+J
        if (It == I+J) then
          do IJ=1,NPrim
            ECoefX(I,J,It,IJ) = ExpHalf(IJ)*ECoefX(I-1,J,It-1,IJ)
            ECoefY(I,J,It,IJ) = ExpHalf(IJ)*ECoefY(I-1,J,It-1,IJ)
            ECoefZ(I,J,It,IJ) = ExpHalf(IJ)*ECoefZ(I-1,J,It-1,IJ)
          end do
        else if (It == I+J-1) then
          do IJ=1,NPrim
            ECoefX(I,J,It,IJ) = ExpHalf(IJ)*ECoefX(I-1,J,It-1,IJ)+PXAX(IJ)*ECoefX(I-1,J,It,IJ)
            ECoefY(I,J,It,IJ) = ExpHalf(IJ)*ECoefY(I-1,J,It-1,IJ)+PYAY(IJ)*ECoefY(I-1,J,It,IJ)
            ECoefZ(I,J,It,IJ) = ExpHalf(IJ)*ECoefZ(I-1,J,It-1,IJ)+PZAZ(IJ)*ECoefZ(I-1,J,It,IJ)
          end do
        else
          do IJ=1,NPrim
            ECoefX(I,J,It,IJ) = ExpHalf(IJ)*ECoefX(I-1,J,It-1,IJ)+PXAX(IJ)*ECoefX(I-1,J,It,IJ)+(It+1)*ECoefX(I-1,J,It+1,IJ)
            ECoefY(I,J,It,IJ) = ExpHalf(IJ)*ECoefY(I-1,J,It-1,IJ)+PYAY(IJ)*ECoefY(I-1,J,It,IJ)+(It+1)*ECoefY(I-1,J,It+1,IJ)
            ECoefZ(I,J,It,IJ) = ExpHalf(IJ)*ECoefZ(I-1,J,It-1,IJ)+PZAZ(IJ)*ECoefZ(I-1,J,It,IJ)+(It+1)*ECoefZ(I-1,J,It+1,IJ)
          end do
        end if
      end do

      if (I /= j) then

        ! * Case of It = 0

        It = 0
        if (I+J <= 1) then
          do IJ=1,NPrim
            ECoefX(J,I,It,IJ) = PXBX(IJ)*ECoefX(J,I-1,It,IJ)
            ECoefY(J,I,It,IJ) = PYBY(IJ)*ECoefY(J,I-1,It,IJ)
            ECoefZ(J,I,It,IJ) = PZBZ(IJ)*ECoefZ(J,I-1,It,IJ)
          end do
        else
          do IJ=1,NPrim
            ECoefX(J,I,It,IJ) = PXBX(IJ)*ECoefX(J,I-1,It,IJ)+ECoefX(J,I-1,It+1,IJ)
            ECoefY(J,I,It,IJ) = PYBY(IJ)*ECoefY(J,I-1,It,IJ)+ECoefY(J,I-1,It+1,IJ)
            ECoefZ(J,I,It,IJ) = PZBZ(IJ)*ECoefZ(J,I-1,It,IJ)+ECoefZ(J,I-1,It+1,IJ)
          end do
        end if

        ! * Case of It >= 1

        do It=1,I+J
          if (It == I+J) then
            do IJ=1,NPrim
              ECoefX(J,I,It,IJ) = ExpHalf(IJ)*ECoefX(J,I-1,It-1,IJ)
              ECoefY(J,I,It,IJ) = ExpHalf(IJ)*ECoefY(J,I-1,It-1,IJ)
              ECoefZ(J,I,It,IJ) = ExpHalf(IJ)*ECoefZ(J,I-1,It-1,IJ)
            end do
          else if (It == I+J-1) then
            do IJ=1,NPrim
              ECoefX(J,I,It,IJ) = ExpHalf(IJ)*ECoefX(J,I-1,It-1,IJ)+PXBX(IJ)*ECoefX(J,I-1,It,IJ)
              ECoefY(J,I,It,IJ) = ExpHalf(IJ)*ECoefY(J,I-1,It-1,IJ)+PYBY(IJ)*ECoefY(J,I-1,It,IJ)
              ECoefZ(J,I,It,IJ) = ExpHalf(IJ)*ECoefZ(J,I-1,It-1,IJ)+PZBZ(IJ)*ECoefZ(J,I-1,It,IJ)
            end do
          else
            do IJ=1,NPrim
              ECoefX(J,I,It,IJ) = ExpHalf(IJ)*ECoefX(J,I-1,It-1,IJ)+PXBX(IJ)*ECoefX(J,I-1,It,IJ)+(It+1)*ECoefX(J,I-1,It+1,IJ)
              ECoefY(J,I,It,IJ) = ExpHalf(IJ)*ECoefY(J,I-1,It-1,IJ)+PYBY(IJ)*ECoefY(J,I-1,It,IJ)+(It+1)*ECoefY(J,I-1,It+1,IJ)
              ECoefZ(J,I,It,IJ) = ExpHalf(IJ)*ECoefZ(J,I-1,It-1,IJ)+PZBZ(IJ)*ECoefZ(J,I-1,It,IJ)+(It+1)*ECoefZ(J,I-1,It+1,IJ)
            end do
          end if
        end do

      end if

    end do
  end do

end subroutine fmm_build_Ecoef1

!-------------------------------------------------------------------------------

subroutine fmm_build_Ecoef2(batch,NPrim,IAngl,JAngl,ECoefX,ECoefY,ECoefZ)

  implicit none
  type(fmm_prim_batch), intent(in) :: batch(:)
  integer(INTK), intent(in)        :: NPrim, IAngl, JAngl
  real(REALK), intent(out)         :: ECoefX(0:,0:,0:,:)
  real(REALK), intent(out)         :: ECoefY(0:,0:,0:,:)
  real(REALK), intent(out)         :: ECoefZ(0:,0:,0:,:)

  integer(INTK) :: I, J, It, IJ, JTemp
  real(REALK) :: PXAX(NPrim), PXBX(NPrim)
  real(REALK) :: PYAY(NPrim), PYBY(NPrim)
  real(REALK) :: PZAZ(NPrim), PZBZ(NPrim)
  real(REALK) :: ExpHalf(NPrim)
  real(REALK) :: PreFact(NPrim)

  ! o Generate E-coefficients (Hermite expansion coefficients) with the three-term recurrence relation

  PXAX(:) = batch(1:Nprim)%PA(1)
  PYAY(:) = batch(1:Nprim)%PA(2)
  PZAZ(:) = batch(1:Nprim)%PA(3)
  PXBX(:) = batch(1:Nprim)%PB(1)
  PYBY(:) = batch(1:Nprim)%PB(2)
  PZBZ(:) = batch(1:Nprim)%PB(3)
  ExpHalf(:) = batch(1:Nprim)%ExpPHalf
  PreFact(:) = batch(1:Nprim)%PreFactAB

  ! o E(0,0,0,NPrim), The prefactor is incorporated.

  ECoefX(0,0,0,1:NPrim) = One
  ECoefY(0,0,0,1:NPrim) = One
  ECoefZ(0,0,0,1:NPrim) = PreFact(1:NPrim)

  ! o Recurrence relation
  !   E(I,J,t,NPrim), 0 <= t <= I + J, IAngl >= JAngl

  if (JAngl == 0) return
  do J=1,JAngl
    JTemp = J
    if (J == JAngl) JTemp = IAngl
    do I=0,JTemp

      ! * Case of It = 0

      It = 0
      if (I+J <= 1) then
        do IJ=1,NPrim
          ECoefX(I,J,It,IJ) = PXBX(IJ)*ECoefX(I,j-1,It,IJ)
          ECoefY(I,J,It,IJ) = PYBY(IJ)*ECoefY(I,j-1,It,IJ)
          ECoefZ(I,J,It,IJ) = PZBZ(IJ)*ECoefZ(I,j-1,It,IJ)
        end do
      else
        do IJ=1,NPrim
          ECoefX(I,J,It,IJ) = PXBX(IJ)*ECoefX(I,j-1,It,IJ)+ECoefX(I,j-1,It+1,IJ)
          ECoefY(I,J,It,IJ) = PYBY(IJ)*ECoefY(I,j-1,It,IJ)+ECoefY(I,j-1,It+1,IJ)
          ECoefZ(I,J,It,IJ) = PZBZ(IJ)*ECoefZ(I,j-1,It,IJ)+ECoefZ(I,j-1,It+1,IJ)
        end do
      end if

      ! * Case of It >= 1

      do It=1,I+J
        if (It == I+J) then
          do IJ=1,NPrim
            ECoefX(I,J,It,IJ) = ExpHalf(IJ)*ECoefX(I,j-1,It-1,IJ)
            ECoefY(I,J,It,IJ) = ExpHalf(IJ)*ECoefY(I,j-1,It-1,IJ)
            ECoefZ(I,J,It,IJ) = ExpHalf(IJ)*ECoefZ(I,j-1,It-1,IJ)
          end do
        else if (It == I+J-1) then
          do IJ=1,NPrim
            ECoefX(I,J,It,IJ) = ExpHalf(IJ)*ECoefX(I,j-1,It-1,IJ)+PXBX(IJ)*ECoefX(I,j-1,It,IJ)
            ECoefY(I,J,It,IJ) = ExpHalf(IJ)*ECoefY(I,j-1,It-1,IJ)+PYBY(IJ)*ECoefY(I,j-1,It,IJ)
            ECoefZ(I,J,It,IJ) = ExpHalf(IJ)*ECoefZ(I,j-1,It-1,IJ)+PZBZ(IJ)*ECoefZ(I,j-1,It,IJ)
          end do
        else
          do IJ=1,NPrim
            ECoefX(I,J,It,IJ) = ExpHalf(IJ)*ECoefX(I,j-1,It-1,IJ)+PXBX(IJ)*ECoefX(I,j-1,It,IJ)+(It+1)*ECoefX(I,j-1,It+1,IJ)
            ECoefY(I,J,It,IJ) = ExpHalf(IJ)*ECoefY(I,j-1,It-1,IJ)+PYBY(IJ)*ECoefY(I,j-1,It,IJ)+(It+1)*ECoefY(I,j-1,It+1,IJ)
            ECoefZ(I,J,It,IJ) = ExpHalf(IJ)*ECoefZ(I,j-1,It-1,IJ)+PZBZ(IJ)*ECoefZ(I,j-1,It,IJ)+(It+1)*ECoefZ(I,j-1,It+1,IJ)
          end do
        end if
      end do

      if (I /= j) then

        ! * Case of It = 0

        It = 0
        if (I+J <= 1) then
          do IJ=1,NPrim
            ECoefX(J,I,It,IJ) = PXAX(IJ)*ECoefX(J-1,I,It,IJ)
            ECoefY(J,I,It,IJ) = PYAY(IJ)*ECoefY(J-1,I,It,IJ)
            ECoefZ(J,I,It,IJ) = PZAZ(IJ)*ECoefZ(J-1,I,It,IJ)
          end do
        else
          do IJ=1,NPrim
            ECoefX(J,I,It,IJ) = PXAX(IJ)*ECoefX(J-1,I,It,IJ)+ECoefX(J-1,I,It+1,IJ)
            ECoefY(J,I,It,IJ) = PYAY(IJ)*ECoefY(J-1,I,It,IJ)+ECoefY(J-1,I,It+1,IJ)
            ECoefZ(J,I,It,IJ) = PZAZ(IJ)*ECoefZ(J-1,I,It,IJ)+ECoefZ(J-1,I,It+1,IJ)
          end do
        end if

        ! * Case of It >= 1

        do It=1,I+J
          if (It == I+J) then
            do IJ=1,NPrim
              ECoefX(J,I,It,IJ) = ExpHalf(IJ)*ECoefX(J-1,I,It-1,IJ)
              ECoefY(J,I,It,IJ) = ExpHalf(IJ)*ECoefY(J-1,I,It-1,IJ)
              ECoefZ(J,I,It,IJ) = ExpHalf(IJ)*ECoefZ(J-1,I,It-1,IJ)
            end do
          else if (It == I+J-1) then
            do IJ=1,NPrim
              ECoefX(J,I,It,IJ) = ExpHalf(IJ)*ECoefX(J-1,I,It-1,IJ)+PXAX(IJ)*ECoefX(J-1,I,It,IJ)
              ECoefY(J,I,It,IJ) = ExpHalf(IJ)*ECoefY(J-1,I,It-1,IJ)+PYAY(IJ)*ECoefY(J-1,I,It,IJ)
              ECoefZ(J,I,It,IJ) = ExpHalf(IJ)*ECoefZ(J-1,I,It-1,IJ)+PZAZ(IJ)*ECoefZ(J-1,I,It,IJ)
            end do
          else
            do IJ=1,NPrim
              ECoefX(J,I,It,IJ) = ExpHalf(IJ)*ECoefX(J-1,I,It-1,IJ)+PXAX(IJ)*ECoefX(J-1,I,It,IJ)+(It+1)*ECoefX(J-1,I,It+1,IJ)
              ECoefY(J,I,It,IJ) = ExpHalf(IJ)*ECoefY(J-1,I,It-1,IJ)+PYAY(IJ)*ECoefY(J-1,I,It,IJ)+(It+1)*ECoefY(J-1,I,It+1,IJ)
              ECoefZ(J,I,It,IJ) = ExpHalf(IJ)*ECoefZ(J-1,I,It-1,IJ)+PZAZ(IJ)*ECoefZ(J-1,I,It,IJ)+(It+1)*ECoefZ(J-1,I,It+1,IJ)
            end do
          end if
        end do

      end if

    end do
  end do

end subroutine fmm_build_Ecoef2

!-------------------------------------------------------------------------------

end module fmm_integral_utils
