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

subroutine DerChg(nAt,nAt3,nTs,nS,Eps,IAtm,AtmC,AtmChg,DM,DerMat,Tessera,Q,Qtot,QDer,DerTes,DerPunt,DerCentr,DerRad,Der1,Der2, &
                  VDer,Sphere,ISphe)

use Constants, only: One, Pi
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nAt, nAt3, nTs, nS, IAtm(nAt), ISphe(*)
real(kind=wp), intent(in) :: Eps, AtmC(3,nAt), AtmChg(nAt), DM(nTs,*), Tessera(4,*), Q(2,*), DerTes(nTs,nAt,3), &
                             DerPunt(nTs,nAt,3,3), DerCentr(nS,nAt,3,3), DerRad(nS,nAt,3), VDer(nTs,*), Sphere(4,*)
real(kind=wp), intent(_OUT_) :: DerMat(nTs,*), Qtot(*), QDer(3,nAt,*), Der1(*), Der2(*)
integer(kind=iwp) :: iAt, iCoord, idx, iTs
real(kind=wp) :: Diag, Sc_Cond

Diag = -1.0694_wp*sqrt(PI)
Sc_Cond = (Eps-One)/Eps

! Total charges

do iTs=1,nTs
  Qtot(iTs) = Q(1,iTs)+Q(2,iTs)
end do

! Compute the derivative of PCM matrix

! Loop on atoms and coordinates
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    ! Conductor-like case
    call DMat_CPCM(iAt,iCoord,Eps,nTs,nS,nAt,Diag,Tessera,DerMat,DerTes,DerPunt,DerCentr,iSphe)

    ! First product: DerDM * Qtot

    call PrMatVec(.false.,.false.,DerMat,-One,nTs,nTs,Qtot,Der1)
    ! pcm_solvent
    !if ((iat == 5) .and.( icoord == 1)) then
    !  write(u6,'(a)') 'In DerChg first contribution for 5, 1'
    !  do its=1,nts
    !    write(u6,'(i4,f20.12)') its,der1(its)
    !  end do
    !end if
    ! pcm_solvent end

    ! Total deriv. of the potential summed up the quantity already
    ! computed (DerDM*Qtot)

    do iTs=1,nTs
      Der1(iTs) = Der1(iTs)+Sc_Cond*VDer(iTs,idx)
    end do

    ! Last product: - DM^-1 * (V^x + DM^x*q)

    ! pcm_solvent
    !if ((iat == 5) .and. (icoord == 1)) then
    !  write(u6,'(a)') 'In DerChg second contribution for 5, 1'
    !  do its=1,nts
    !    write(u6,'(i4,f20.12)') its,der1(its)
    !  end do
    !end if
    ! pcm_solvent end
    call PrMatVec(.false.,.false.,DM,-One,nTs,nTs,Der1,Der2)
    call FillQDer(nAt,nTs,iAt,iCoord,Der2,QDer)
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nAt3)
  call Unused_integer_array(IAtm)
  call Unused_real_array(AtmC)
  call Unused_real_array(AtmChg)
  call Unused_real_array(DerRad)
  call Unused_real_array(Sphere)
end if

end subroutine DerChg
!====
subroutine DMat_CPCM(iAt,iC,Eps,nTs,nS,nAt,fact,Tessera,DerMat,DerTes,DerPunt,DerCentr,iSphe)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAt, iC, nTs, nS, nAt, iSphe(*)
real(kind=wp), intent(in) :: Eps, fact, Tessera(4,*), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), DerCentr(nS,nAt,3,3)
real(kind=wp), intent(_OUT_) :: DerMat(nTs,*)
integer(kind=iwp) :: ITs, JTs, L, LJ
real(kind=wp) :: DIJ, DXIJ, DYIJ, DZIJ, PROD, XIJ, YIJ, ZIJ

! Compute the derivative of the CPCM matrix wrt atom iat, coord. ic

! Loop on tesserae
do ITs=1,NTs
  L = ISPHE(ITs)
  do JTs=1,NTs
    LJ = ISPHE(JTs)
    ! Diagonal elements
    if (ITs == JTs) then
      DerMat(ITs,ITs) = fact*DERTES(ITs,iAt,IC)/(Tessera(4,ITs)*sqrt(Tessera(4,ITs)))
    else
      ! Off diagonal elements
      XIJ = Tessera(1,ITs)-Tessera(1,JTs)
      YIJ = Tessera(2,ITs)-Tessera(2,JTs)
      ZIJ = Tessera(3,ITs)-Tessera(3,JTs)
      DIJ = sqrt(XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ)
      DXIJ = DERPUNT(ITs,iAt,IC,1)+DERCENTR(L,iAt,IC,1)-DERPUNT(JTs,iAt,IC,1)-DERCENTR(LJ,iAt,IC,1)
      DYIJ = DERPUNT(ITs,iAt,IC,2)+DERCENTR(L,iAt,IC,2)-DERPUNT(JTs,iAt,IC,2)-DERCENTR(LJ,iAt,IC,2)
      DZIJ = DERPUNT(ITs,iAt,IC,3)+DERCENTR(L,iAt,IC,3)-DERPUNT(JTs,iAt,IC,3)-DERCENTR(LJ,iAt,IC,3)
      PROD = (XIJ*DXIJ+YIJ*DYIJ+ZIJ*DZIJ)/DIJ**3
      DerMat(ITs,JTs) = -PROD
    end if
  end do
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_real(Eps)

end subroutine DMat_CPCM
!====
subroutine PrMatVec(Dag,DoSym,Mat,f,n,m,Vec,Res)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
logical(kind=iwp), intent(in) :: Dag, DoSym
integer(kind=iwp), intent(in) :: n, m
real(kind=wp), intent(in) :: Mat(n,*), f, Vec(*)
real(kind=wp), intent(_OUT_) :: Res(*)
integer(kind=iwp) :: i, j
real(kind=wp) :: Elm

! Do the matrix vector product: Res(n,1)=f*Mat(n,m)*Vec(m,1)
! possibly transposed (if Dag): Res(1,n)=f*Vec(1,m)*Mat(m,n)
! If DoSym symmetrize the matrix elements

do i=1,n
  Res(i) = Zero
  do j=1,m
    if (DoSym) then
      Elm = (Mat(i,j)+Mat(j,i))/Two
    else
      if (Dag) Elm = Mat(j,i)
      if (.not. Dag) Elm = Mat(i,j)
    end if
    Res(i) = Res(i)+f*Elm*Vec(j)
  end do
end do

return

end subroutine PrMatVec
!====
subroutine FillQDer(nAt,nTs,iAt,iC,Der,QDer)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs, iAt, iC
real(kind=wp), intent(in) :: Der(*)
real(kind=wp), intent(_OUT_) :: QDer(3,nAt,*)
integer(kind=iwp) :: iTs

do iTs=1,nTs
  QDer(iC,iAt,iTs) = Der(iTs)
end do

return

end subroutine FillQDer
!====
subroutine testq(nAt,nTs,VDer,Q,QTot)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs
real(kind=wp), intent(in) :: Q(2,*)
real(kind=wp), intent(_OUT_) :: VDer(nTs,*), QTot(*)
integer(kind=iwp) :: iAt, iCoord, idx, iTs, Lu
real(kind=wp) :: rsum

Lu = 1
call Molcas_open(Lu,'DerPt.dat')
!open(1,file='DerPot.dat',status='old',form='formatted')
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    do iTs=1,nTs
      read(1,*) VDer(iTs,idx)
    end do
  end do
end do
close(1)
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    rsum = Zero
    do iTs=1,nts
      QTot(iTs) = Q(1,iTs)+Q(2,its)
      rsum = rsum+QTot(iTs)*VDer(iTs,idx)
    end do
    write(u6,'("Charges times VDer",i4,f20.12)') idx,rsum
  end do
end do

return

end subroutine testq
