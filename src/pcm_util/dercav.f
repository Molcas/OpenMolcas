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
      Subroutine DerCav(ToAng,nTs,nAt,nS,nAt3,Eps,Tessera,Q,Qtot,Der1,
     $  Dertes,DerPunt,DerCentr,DerRad,QDer,Sphere,iSphe,
     $  nOrd)
      Implicit Real*8(a-h,o-z)
      Dimension Sphere(4,*),iSphe(*),nOrd(*)
      Dimension Tessera(4,*),Q(2,*),Qtot(*)
      Dimension Der1(*)
      Dimension DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3)
      Dimension DerRad(nS,nAt,3),DerCentr(nS,nAt,3,3)
      Dimension QDer(3,nAt,*)
      Save Zero,Two
      Data Zero/0.0d0/,Two/2.0d0/
*
*---- Derivative of the cavity factor U_x(q)=2 Pi Eps/(Eps-1) sum_i [Qtot**2 * n_x]
*
      dN=Zero ! Dummy initialization.
c     Double loop on atoms and coordinates
      Do 100 Index1 = 1, nAt3
        iAt1 = Int( (Index1-1)/3 ) + 1
        iCoord1 = Index1 - 3 * (iAt1-1)
        Do 101 Index2 = 1, nAt3
          iAt2 = Int( (Index2-1)/3 ) + 1
          iCoord2 = Index2 - 3 * (iAt2-1)
c         Derivative of the normal factor n_x
          Call Der_Norm(ToAng,iAt1,iCoord1,iAt2,iCOord2,nTs,nAt,nS,
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
            DerQ = QDer(iCoord1,iAt1,iTs)
            Sum1 = Sum1 + Two * Qtot(iTs) * DerQ * dN / Tessera(4,iTs)
            Sum2 = Sum2 + Qtot(iTs) * Qtot(iTs) * Der1(iTs)
  200     Continue
  101   Continue
  100 Continue
cpcm_solvent
c      write(6,
c    & '(''In DerCav tesserae coord., area, qnuc and qel'')')
c     Do Its = 1, nTs
c       write(6,'(6f12.6)')(Tessera(j,iTs),j=1,4),q(1,iTs),q(2,iTs)
c     EndDo
cpcm_solvent end
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_(Eps)
         Call Unused_real_array(Q)
      End If
      End

      Subroutine Der_Norm(ToAng,iAt1,iCoord1,iAt2,iCOord2,nTs,nAt,nS,
     &  Tessera,Der1,DerRad,DerTes,DerPunt,Sphere,iSphe,nOrd)
      Implicit Real*8(A-H,O-Z)
      Dimension Der1(*),Sphere(4,*),iSphe(*),nOrd(*),Tessera(4,*)
      Dimension DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3)
      Dimension DerRad(nS,nAt,3)
      Data Zero, One /0.0d0, 1.0d0/
c
      AnToAU = One / ToAng
c     Find out if atom iAt1 has a sphere around
      iAt1_S = 0
      Do iS = 1, nS
        If(iAt1.eq.nOrd(iS)) iAt1_S = iS
      EndDo
*
*---- Compute the derivative of the normal vector over the tessera area
*
c     Loop on tiles
      dN = Zero ! Dummy initialization.
      Do 100 iTs = 1, nTs
        L = iSphe(iTs)
        Der1(iTs) = Zero
        If(L.eq.iAt1_S) then
          If(iCoord1.eq.1) dN = (Sphere(1,L)-Tessera(1,iTs))/Sphere(4,L)
          If(iCoord1.eq.2) dN = (Sphere(2,L)-Tessera(2,iTs))/Sphere(4,L)
          If(iCoord1.eq.3) dN = (Sphere(3,L)-Tessera(3,iTs))/Sphere(4,L)
          dN_Its = DerPunt(iTs,iAt2,iCoord2,iCoord1)
     &           + dN * DerRad(L,iAt2,iCoord2)
          dN_Its = - dN_Its / Sphere(4,L)
        else
          dN = Zero
          dN_Its = Zero
        EndIf
        SqArea = Tessera(4,iTs) * Tessera(4,iTs)
        Der_1 = - dN * DerTes(iTs,iAt2,iCoord2) * AnToAU / SqArea
        Der1(its) = - dN_Its / Tessera(4,iTs) - Der_1
  100  Continue
      Return
      End
