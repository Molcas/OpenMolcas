************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Cav_Hss(nAt,nAt3,nTs,nS,Eps,Sphere,iSphe,nOrd,
     &     Tessera,Q,DM,Der1,DerDM,Temp,DerTes,DerPunt,DerRad,DerCentr,
     &     Hess,nHess)
      Implicit Real*8(a-h,o-z)
      Dimension Sphere(4,*),iSphe(*),nOrd(*)
      Dimension Tessera(4,*),Q(2,*),Der1(*)
      Dimension DM(nTs,*),DerDM(nTs,*),Temp(nTs,*)
      Dimension DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3)
      Dimension DerRad(nS,nAt,3),DerCentr(nS,nAt,3,3)
      Dimension Hess(nAt3,*)
#include "angstr.fh"
#include "real.fh"
*
*---- Derivative of the cavity factor U_x(q)=2 Pi Eps/(Eps-1) sum_i [Qtot**2 * n_x]
*
      Fact = Two * PI * Eps / (Eps-One)
      FPI = Four * PI
      Diag = - 1.0694d0 * Sqrt(FPI) / Two
      dN = Zero
c     Double loop on atoms and coordinates
      Do 100 Index1 = 1, nAt3
        iAt1 = Int( (Index1-1)/3 ) + 1
        iCoord1 = Index1 - 3 * (iAt1-1)
c       Derivative of the PCM matrix
        Call DMat_CPCM(iAt1,iCoord1,Eps,nTs,nS,nAt,Diag,Tessera,
     &                 DerDM,DerTes,DerPunt,DerCentr,iSphe)
c       Matrix product: derivative of PCM matrix times the inverted matrix
        Call DGEMM_('N','N',nTs,nTs,nTs,One,DerDM,nTs,
     &             DM,nTs,Zero,Temp,nTs)

        Do 101 Index2 = 1, nAt3
          iAt2 = Int( (Index2-1)/3 ) + 1
          iCoord2 = Index2 - 3 * (iAt2-1)
c         Derivative of the normal factor n_x
          Call Der_Norm(Angstr,iAt1,iCoord1,iAt2,iCOord2,nTs,nAt,nS,
     &    Tessera,Der1,DerRad,DerTes,DerPunt,Sphere,iSphe,nOrd)
c         Find out if atom iAt2 has a sphere around
          iAt2_S = 0
          Do iS = 1, nS
           If(iAt2.eq.nOrd(iS)) iAt2_S = iS
          EndDo
          Sum1 = Zero
          Sum2 = Zero
c         Loop on tesserae
          Do 200 iTs = 1, nTs
            QtotI = Q(1,iTs) + Q(2,iTs)
            L = iSphe(iTs)
            XN = - (Sphere(1,L) - Tessera(1,iTs)) / Sphere(4,L)
            YN = - (Sphere(2,L) - Tessera(2,iTs)) / Sphere(4,L)
            ZN = - (Sphere(3,L) - Tessera(3,iTs)) / Sphere(4,L)
            If(L.eq.iAt2_S) then
             if(iCoord2.eq.1)dN = XN
             if(iCoord2.eq.2)dN = YN
             if(iCoord2.eq.3)dN = ZN
            Else
              dCent = XN * DerCentr(L,iAt2,iCoord2,1)
     &              + YN * DerCentr(L,iAt2,iCoord2,2)
     &              + ZN * DerCentr(L,iAt2,iCoord2,3)
              dN = DerRad(L,iAt2,iCoord2) + dCent
            EndIf
            dNI = dN / Tessera(4,iTs)
            Sum1 = Sum1 + QtotI * QtotI * Der1(iTs)
            Do  201 jTs = 1, nTs
              QtotJ = Q(1,jTs) + Q(2,jTs)
              Sum2 = Sum2 + Two * QtotI * dNI * Temp(iTs,jTs) * QtotJ
  201       Continue
  200     Continue
          Hess(Index1,Index2) = Fact * (Sum1 + Sum2)
  101   Continue
  100 Continue
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nHess)
      End
