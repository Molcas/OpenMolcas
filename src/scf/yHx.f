************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2022, Roland Lindh                                     *
************************************************************************
      SubRoutine yHx(X,Y,nXY)
************************************************************************
*                                                                      *
*     purpose: multiply an approximation of the orbital Hessian times  *
*              a trial vector, X.                                      *
*                                                                      *
*     Y = H(approximate) x  X                                          *
*                                                                      *
************************************************************************
      use Orb_Type, only: OrbType
      use InfSCF, only: nSym, nFro, nOrb, nOcc
      use SCF_Arrays, only: FockMO
      use Constants, only: Zero, Four
      Implicit None
      Integer nXY
      Real*8, Target:: X(nXY), Y(nXY)
*
*     declaration local variables
      Integer nD, iD
      Integer iSym,iOcc,iVir,nOccmF,nOrbmF, iOff_F
      Integer jOcc, jVir, iOff_XY
      Real*8 Tmp
      Real*8, Parameter:: Hii_Min=0.05D0
      Real*8, Parameter:: Hii_Max=1.00D0
      Real*8, Pointer:: Fock(:,:), XP(:,:), YP(:,:)
*
*----------------------------------------------------------------------*
*

*
      nD = Size(FockMO,2)
      iOff_XY = 0
      Do iD = 1, nD
*
         iOff_F=0
         Do iSym=1,nSym
*
*            loop over all occ orbitals in sym block
*
             ! number of Occupied, excluding frozen
             nOccmF=nOcc(iSym,iD)-nFro(iSym)
             ! number of Orbitals, excluding frozen
             nOrbmF=nOrb(iSym)-nFro(iSym)

             Fock(1:nOrb(iSym),1:nOrb(iSym)) =>
     &            FockMO(iOff_F+1:iOff_F+nOrb(iSym)**2,iD)
             XP(1:nOccmF,nOccmF:nOrbmF) =>
     &            X(iOff_XY+1:iOff_XY+nOccmF*(nOrbmF-nOccmF))
             YP(1:nOccmF,nOccmF:nOrbmF) =>
     &            Y(iOff_XY+1:iOff_XY+nOccmF*(nOrbmF-nOccmF))
*
             Do iOcc= 1, nOccmF
                Do iVir=nOccmF+1,nOrbmF
*
                   Tmp = Zero
                   Do jOcc = 1, nOccmF
                      Do jVir=nOccmF+1,nOrbmF

                         If (OrbType(iVir,iD).eq.OrbType(iOcc,iD) .and.
     &                       OrbType(jVir,iD).eq.OrbType(jOcc,iD) .and.
     &                       OrbType(iVir,iD).eq.OrbType(jOcc,iD)) Then

                         If (iVir==jVir .and. iOcc==jOcc) Then

                            Tmp = Tmp + (Four
     &                          * (Fock(jVir,jVir)-Fock(jOcc,jOcc))
     &                          / DBLE(nD)) * XP(jOcc,jVir)

                            If (Tmp<Zero) Then
                               Write (6,*) 'Hii<0.0, Hii=',Tmp
                               Tmp=Max(Hii_Max,Abs(Tmp))
                            Else If (Abs(Tmp).lt.Hii_Min) Then
                               Tmp=Hii_Min
                               Write (6,*) 'Abs(Hii)<0.05'
                            End If

                         End If

                         End If

                      End Do  ! jVir
                   End Do     ! jOcc
                   YP(iOcc,iVir) = Tmp


                End Do  ! iVir
             End Do     ! iOcc
*
             Nullify(Fock,XP,YP)
             iOff_XY= iOff_XY + nOccmF*(nOrbmF-nOccmF)
             iOff_F = iOff_F + nOrb(iSym)**2
*
          End Do ! iSym
      End Do ! iD

      End Subroutine yHx
