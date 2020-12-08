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
      Subroutine ConInt(xyz,nCent,dE,Bf,lWrite_,Label,dBf,ldB,lIter)
      use Slapaf_Info, only: Gx, Gx0, Energy, Energy0
      Implicit Real*8  (a-h,o-z)
#include "info_slapaf.fh"
#include "real.fh"
#include "nadc.fh"
#include "constants.fh"
      Real*8   Bf(3,nCent), xyz(3,nCent), dBf(3*nCent,3*nCent)
      Logical lWrite_, ldB
      Character(LEN=8) Label
*
*
      E1 = Energy (lIter)
      E0 = Energy0(lIter)
*
c     iOpt=1 -> Linear
c     iOpt=2 -> Quadratic
c     iOpt=3 -> Absolute value
      If (NADC) Then
         If (ApproxNADC) Then
            iOpt=2
         Else
            iOpt=3
         End If
      Else
         iOpt=1
      End If
*
*     Observe that the program is storing the forces rather than the
*     gradients!
*
*     For a true conical intersection the storage is done a bit
*     differently (see init2.f). Here the average energy is stored
*     in E1 and the energy difference in E0. Ditto for the gradients.
*
      dE=Zero
      If (iOpt.eq.1) Then
C------- Linear ------------------
         If (NADC) Then
            dE=E0
         Else
            dE=E1-E0
         End If
      Else If (iOpt.eq.2) Then
C------- Quadratic ---------------
         If (NADC) Then
            dE=E0**2
         Else
            dE=(E1-E0)**2
         End If
      Else If (iOpt.eq.3) Then
C------- Absolute value ----------
         If (NADC) Then
            dE=Abs(E0)
         Else
            dE=Abs(E1-E0)
         End IF
      End If
      If (lWrite_) Then
         If (NADC) Then
            Write (6,'(2A,F18.8,A,F18.8,A)')
     &                    Label,' : Energy difference = ',
     &                    E0, ' hartree, ',
     &                    E0*CONV_AU_TO_KJ_PER_MOLE_,
     &                    ' kJ/mol'
            Write (6,'( A,F18.8,A)') '           Average energy    = ',
     &                               E1   , ' hartree'
         Else
            Write (6,'(2A,F18.8,A,F18.8,A)')
     &                    Label,' : Energy difference = ',
     &                    E1-E0, ' hartree, ',
     &                    (E1-E0)*CONV_AU_TO_KJ_PER_MOLE_,
     &                    ' kJ/mol'
            Write (6,'( A,F18.8,A)') '           E(i)              = ',
     &                               E1   , ' hartree'
            Write (6,'( A,F18.8,A)') '           E(j)              = ',
     &                               E0   , ' hartree'
         End If
      End If
*
*---- Compute the WDC B-matrix
*
C     Call RecPrt('Grad1',' ',Gx(1,1,lIter),3,nCent)
C     Call RecPrt('Grad0',' ',Gx0(1,1,lIter),3,nCent)
      Do iCent = 1, nCent
         Fact=DBLE(iDeg(xyz(1,iCent)))
C        Write (6,*) 'Fact=',Fact
         Do iCar = 1, 3
            Bf(iCar,iCent)=Zero
            If (iOpt.eq.1) Then
C------------- Linear ------------------
               If (NADC) Then
                  Bf(iCar,iCent)=-Gx0(iCar,iCent,lIter)
               Else
                  Bf(iCar,iCent)=-(Gx(iCar,iCent,lIter)
     &                          -Gx0(iCar,iCent,lIter))
               End If
            Else If (iOpt.eq.2) Then
C------------- Quadratic ---------------
*
*              When the energy difference becomes small, the true derivative vanishes.
*              In such case use simply a scaled-down energy difference gradient
*
               If (NADC) Then
                  If (Abs(E0).gt.1.0D-5) Then
                     Bf(iCar,iCent)=-Two*E0
     &                             *Gx0(iCar,iCent,lIter)
                  Else
                     Bf(iCar,iCent)=-Two*1.0D-5
     &                             *Gx0(iCar,iCent,lIter)
                  End If
               Else
                  If (Abs(E1-E0).gt.1.0D-5) Then
                     Bf(iCar,iCent)=-Two*(E1-E0)
     &                             *(Gx(iCar,iCent,lIter)
     &                              -Gx0(iCar,iCent,lIter))
                  Else
                     Bf(iCar,iCent)=-Two*1.0D-5
     &                             *(Gx(iCar,iCent,lIter)
     &                              -Gx0(iCar,iCent,lIter))
                  End If
               End If
            Else If (iOpt.eq.3) Then
C------------- Absolute value ----------
               If (NADC) Then
                  Bf(iCar,iCent)=-Sign(One,E0)*Gx0(iCar,iCent,lIter)
               Else
                  Bf(iCar,iCent)=-Sign(One,E1-E0)
     &                             *(Gx(iCar,iCent,lIter)
     &                              -Gx0(iCar,iCent,lIter))
               End If
            End If
*
            Bf(iCar,iCent)=Fact*Bf(iCar,iCent)
         End Do
      End Do
*     Call RecPrt('Bf',' ',Bf,3,nCent)
      If (lWrite_.and.iOpt.eq.1) Then
         XX=Sqrt(DDot_(3*nCent,Bf,1,Bf,1))
         If (XX.le.1.0D-3) Then
            Write (6,*)
            Write (6,*)
     &            '    Warning: PESs might be parallel!'
            Write (6,*)
         End If
      End If
*
*---- Compute the cartesian derivative of the B-Matrix.
*
      If (ldB) Then
         nGrad=3*nCent
         If (iOpt.eq.1) Then
C---------- Linear ------------------
            Call FZero(dBf,(3*nCent)**2)
         Else If (iOpt.eq.2) Then
C---------- Quadratic ---------------
            Call FZero(dBf,(3*nCent)**2)

            ix = 0
            Do iCent = 1, nCent
            Do i   = 1, 3
               ix = ix + 1

               iy=0
               Do jCent = 1, nCent
               Do j   = 1, 3
                 iy = iy + 1
                  If (NADC) Then
                     dBf(ix,iy)=-Two*
     &                               Gx0(i,iCent,lIter)
     &                             * Gx0(j,jCent,lIter)
                  Else
                     dBf(ix,iy)=-Two*
     &                              (Gx(i,iCent,lIter)-
     &                               Gx0(i,iCent,lIter))
     &                             *(Gx(j,jCent,lIter)-
     &                               Gx0(j,jCent,lIter))
                  End If
               End Do
               End Do
            End Do
            End Do

         Else If (iOpt.eq.3) Then
C------------- Absolute value ----------
            Call FZero(dBf,(3*nCent)**2)
*
         End If
*        Call RecPrt('dBf','(9F9.1)',dBf,3*nCent,3*nCent)
*
      End If
*
      Return
      End
