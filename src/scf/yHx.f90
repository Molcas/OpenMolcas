!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine yHx(X,Y,nXY)
!***********************************************************************
!                                                                      *
!     purpose: multiply an approximation of the orbital Hessian times  *
!              a trial vector, X.                                      *
!                                                                      *
!     Y = H(approximate) x  X                                          *
!                                                                      *
!***********************************************************************
      use Orb_Type, only: OrbType
      use InfSCF, only: nSym, nFro, nOrb, nOcc
      use SCF_Arrays, only: FockMO
      use Constants, only: Zero, Four
      Implicit None
      Integer nXY
      Real*8, Target:: X(nXY), Y(nXY)
!
!     declaration local variables
      Integer nD, iD
      Integer iSym,iOcc,iVir,nOccmF,nOrbmF, iOff_F
      Integer jOcc, jVir, iOff_XY
      Real*8 Tmp, Hij
      Real*8, Parameter:: Hii_Min=0.05D0
      Real*8, Parameter:: Hii_Max=1.00D0
      Real*8, Pointer:: Fock(:,:), XP(:,:), YP(:,:)
!
!----------------------------------------------------------------------*
!
!     Write (6,*)
!     Call NrmClc(FockMO,SIZE(FockMO),'yHx','FockMO(:,:)')
!     Call NrmClc(X,SIZE(X),'yHx','X(:)')
!     Call RecPrt('yHx: FockMO',' ',FockMO,Size(FockMO,1),Size(FockMO,2))
!     Call RecPrt('yHx: X',' ',X,1,Size(X))

!
      nD = Size(FockMO,2)
      iOff_XY = 0
      Do iD = 1, nD
!
         iOff_F=0
         Do iSym=1,nSym
!
!            loop over all occ orbitals in sym block
!
             ! number of Occupied, excluding frozen
             nOccmF=nOcc(iSym,iD)-nFro(iSym)
             ! number of Orbitals, excluding frozen
             nOrbmF=nOrb(iSym)-nFro(iSym)

             Fock(1:nOrb(iSym),1:nOrb(iSym)) => FockMO(iOff_F+1:iOff_F+nOrb(iSym)**2,iD)
             XP(nOccmF+1:nOrbmF,1:nOccmF) => X(iOff_XY+1:iOff_XY+nOccmF*(nOrbmF-nOccmF))
             YP(nOccmF+1:nOrbmF,1:nOccmF) => Y(iOff_XY+1:iOff_XY+nOccmF*(nOrbmF-nOccmF))
!
             Do iOcc= 1, nOccmF
                Do iVir=nOccmF+1,nOrbmF
!
                   Tmp = Zero
                   Do jOcc = 1, nOccmF
                      Do jVir=nOccmF+1,nOrbmF

                         Hij = Zero
                         If (OrbType(iVir,iD).eq.OrbType(iOcc,iD) .and.   &
                             OrbType(jVir,iD).eq.OrbType(jOcc,iD) .and.   &
                             OrbType(iVir,iD).eq.OrbType(jOcc,iD)) Then

                         If (iVir==jVir .and. iOcc==jOcc) Then
                            Hij = ( Four * (Fock(iVir,jVir)-Fock(iOcc,jOcc)) / DBLE(nD) )

                            If (Hij<Zero) Then
!                              Write (6,*) 'Hii<0.0, Hii=',Hij
                               Hij=Max(Hii_Max,Abs(Hij))
                            Else If (Abs(Hij)<Hii_Min) Then
!                              Write (6,*) 'Abs(Hii)<0.05, Hii=',Hij
!                              Write (6,*) 'jVir,jOcc=',jVir,jOcc
!                              Write (6,*) 'Fock(jOcc,jOcc)=', Fock(jOcc,jOcc)
!                              Write (6,*) 'Fock(jVir,jVir)=', Fock(jVir,jVir)
                               Hij=Hii_Min
                            End If

                         Else If (iVir==jVir .and. iOcc/=jOcc) Then
                            Hij = ( Four * (               -Fock(iOcc,jOcc)) / DBLE(nD) )
                         Else If (iOcc==jOcc .and. iVir/=jVir) Then
                            Hij = ( Four * (Fock(iVir,jVir)                ) / DBLE(nD) )
                         End If
                         Tmp = Tmp + Hij * XP(jVir,jOcc)

                         End If

                      End Do  ! jVir
                   End Do     ! jOcc
                   YP(iVir,iOcc) = Tmp


                End Do  ! iVir
             End Do     ! iOcc
!
             Nullify(Fock,XP,YP)
             iOff_XY= iOff_XY + nOccmF*(nOrbmF-nOccmF)
             iOff_F = iOff_F + nOrb(iSym)**2
!
          End Do ! iSym
      End Do ! iD
#ifdef _DEBUGPRINT_
      Call NrmClc(Y,SIZE(Y),'yHx','Y(:)')
      Call RecPrt('yHx: Y',' ',Y,1,Size(Y))
#endif
      End Subroutine yHx
