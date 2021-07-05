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

subroutine Cav_Hss(nAt,nAt3,nTs,nS,Eps,Sphere,iSphe,nOrd,Tessera,Q,DM,Der1,DerDM,Temp,DerTes,DerPunt,DerRad,DerCentr,Hess,nHess)

implicit real*8(a-h,o-z)
dimension Sphere(4,*), iSphe(*), nOrd(*)
dimension Tessera(4,*), Q(2,*), Der1(*)
dimension DM(nTs,*), DerDM(nTs,*), Temp(nTs,*)
dimension DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3)
dimension DerRad(nS,nAt,3), DerCentr(nS,nAt,3,3)
dimension Hess(nAt3,*)
#include "angstr.fh"
#include "real.fh"

! Derivative of the cavity factor U_x(q)=2 Pi Eps/(Eps-1) sum_i [Qtot**2 * n_x]

Fact = Two*PI*Eps/(Eps-One)
FPI = Four*PI
Diag = -1.0694d0*sqrt(FPI)/Two
dN = Zero
! Double loop on atoms and coordinates
do Index1=1,nAt3
  iAt1 = int((Index1-1)/3)+1
  iCoord1 = Index1-3*(iAt1-1)
  ! Derivative of the PCM matrix
  call DMat_CPCM(iAt1,iCoord1,Eps,nTs,nS,nAt,Diag,Tessera,DerDM,DerTes,DerPunt,DerCentr,iSphe)
  ! Matrix product: derivative of PCM matrix times the inverted matrix
  call DGEMM_('N','N',nTs,nTs,nTs,One,DerDM,nTs,DM,nTs,Zero,Temp,nTs)

  do Index2=1,nAt3
    iAt2 = int((Index2-1)/3)+1
    iCoord2 = Index2-3*(iAt2-1)
    ! Derivative of the normal factor n_x
    call Der_Norm(Angstr,iAt1,iCoord1,iAt2,iCOord2,nTs,nAt,nS,Tessera,Der1,DerRad,DerTes,DerPunt,Sphere,iSphe,nOrd)
    ! Find out if atom iAt2 has a sphere around
    iAt2_S = 0
    do iS=1,nS
      if (iAt2 == nOrd(iS)) iAt2_S = iS
    end do
    Sum1 = Zero
    Sum2 = Zero
    ! Loop on tesserae
    do iTs=1,nTs
      QtotI = Q(1,iTs)+Q(2,iTs)
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
      dNI = dN/Tessera(4,iTs)
      Sum1 = Sum1+QtotI*QtotI*Der1(iTs)
      do jTs=1,nTs
        QtotJ = Q(1,jTs)+Q(2,jTs)
        Sum2 = Sum2+Two*QtotI*dNI*Temp(iTs,jTs)*QtotJ
      end do
    end do
    Hess(Index1,Index2) = Fact*(Sum1+Sum2)
  end do
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nHess)

end subroutine Cav_Hss
