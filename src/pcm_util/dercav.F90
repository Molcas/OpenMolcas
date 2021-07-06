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

subroutine DerCav(ToAng,nTs,nAt,nS,nAt3,Eps,Tessera,Q,Qtot,Der1,DerTes,DerPunt,DerCentr,DerRad,QDer,Sphere,iSphe,nOrd)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nTs, nAt, nS, nAt3, iSphe(*), nOrd(*)
real(kind=wp), intent(in) :: ToAng, Eps, Tessera(4,*), Q(2,*), QTot(*), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), &
                             DerCentr(nS,nAt,3,3), DerRad(nS,nAt,3), QDer(3,nAt,*), Sphere(4,*)
real(kind=wp), intent(_OUT_) :: Der1(*)
integer(kind=iwp) :: iAt1, iAt2, iAt2_S, iCoord1, iCoord2, Index1, Index2, iS, iTs, L
real(kind=wp) :: dCent, DerQ, dN, Sum1, Sum2, XN, YN, ZN

! Derivative of the cavity factor U_x(q)=2 Pi Eps/(Eps-1) sum_i [Qtot**2 * n_x]

dN = Zero ! Dummy initialization.
! Double loop on atoms and coordinates
do Index1=1,nAt3
  iAt1 = int((Index1-1)/3)+1
  iCoord1 = Index1-3*(iAt1-1)
  do Index2=1,nAt3
    iAt2 = int((Index2-1)/3)+1
    iCoord2 = Index2-3*(iAt2-1)
    ! Derivative of the normal factor n_x
    call Der_Norm(ToAng,iAt1,iCoord1,iAt2,iCOord2,nTs,nAt,nS,Tessera,Der1,DerRad,DerTes,DerPunt,Sphere,iSphe,nOrd)
    ! Find out if atom iAt2 has a sphere around
    iAt2_S = 0
    do iS=1,nS
      if (iAt2 == nOrd(iS)) iAt2_S = iS
    end do
    Sum1 = Zero
    Sum2 = Zero
    ! Loop on tesserae
    do iTs=1,nTs
      L = iSphe(iTs)
      XN = -(Sphere(1,L)-Tessera(1,iTs))/Sphere(4,L)
      YN = -(Sphere(2,L)-Tessera(2,iTs))/Sphere(4,L)
      ZN = -(Sphere(3,L)-Tessera(3,iTs))/Sphere(4,L)
      if (L == iAt2_S) then
        if (iCoord2 == 1) dN = XN
        if (iCoord2 == 2) dN = YN
        if (iCoord2 == 3) dN = ZN
      else
        dCent = XN*DerCentr(L,iAt2,iCoord2,1)+YN*DerCentr(L,iAt2,iCoord2,2)+ZN*DerCentr(L,iAt2,iCoord2,3)
        dN = DerRad(L,iAt2,iCoord2)+dCent
      end if
      DerQ = QDer(iCoord1,iAt1,iTs)
      Sum1 = Sum1+Two*Qtot(iTs)*DerQ*dN/Tessera(4,iTs)
      Sum2 = Sum2+Qtot(iTs)*Qtot(iTs)*Der1(iTs)
    end do
  end do
end do
! pcm_solvent
!write(u6,'(a)') 'In DerCav tesserae coord., area, qnuc and qel'
!do Its=1,nTs
!  write(u6,'(6f12.6)') (Tessera(j,iTs),j=1,4),q(1,iTs),q(2,iTs)
!end do
! pcm_solvent end

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_(Eps)
  call Unused_real_array(Q)
end if

end subroutine DerCav
!====
subroutine Der_Norm(ToAng,iAt1,iCoord1,iAt2,iCoord2,nTs,nAt,nS,Tessera,Der1,DerRad,DerTes,DerPunt,Sphere,iSphe,nOrd)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iAt1, iCoord1, iAt2, iCoord2, nTs, nAt, nS, iSphe(*), nOrd(*)
real(kind=wp), intent(in) :: ToAng, Tessera(4,*), DerRad(nS,nAt,3), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), Sphere(4,*)
real(kind=wp), intent(_OUT_) :: Der1(*)
integer(kind=iwp) :: iAt1_S, iS, iTs, L
real(kind=wp) :: AnToAU, Der_1, dN, dN_Its, SqArea

AnToAU = One/ToAng
! Find out if atom iAt1 has a sphere around
iAt1_S = 0
do iS=1,nS
  if (iAt1 == nOrd(iS)) iAt1_S = iS
end do

! Compute the derivative of the normal vector over the tessera area

! Loop on tiles
dN = Zero ! Dummy initialization.
do iTs=1,nTs
  L = iSphe(iTs)
  Der1(iTs) = Zero
  if (L == iAt1_S) then
    if (iCoord1 == 1) dN = (Sphere(1,L)-Tessera(1,iTs))/Sphere(4,L)
    if (iCoord1 == 2) dN = (Sphere(2,L)-Tessera(2,iTs))/Sphere(4,L)
    if (iCoord1 == 3) dN = (Sphere(3,L)-Tessera(3,iTs))/Sphere(4,L)
    dN_Its = DerPunt(iTs,iAt2,iCoord2,iCoord1)+dN*DerRad(L,iAt2,iCoord2)
    dN_Its = -dN_Its/Sphere(4,L)
  else
    dN = Zero
    dN_Its = Zero
  end if
  SqArea = Tessera(4,iTs)*Tessera(4,iTs)
  Der_1 = -dN*DerTes(iTs,iAt2,iCoord2)*AnToAU/SqArea
  Der1(its) = -dN_Its/Tessera(4,iTs)-Der_1
end do

return

end subroutine Der_Norm
