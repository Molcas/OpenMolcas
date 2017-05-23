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
      Subroutine ConInt(xyz,nCent,dE,Bf,lWrite_,lWarn,Label,dBf,ldB,
     &                  lIter)
      Implicit Real*8  (a-h,o-z)
#include "info_slapaf.fh"
#include "real.fh"
#include "nadc.fh"
#include "WrkSpc.fh"
#include "constants.fh"
      Real*8   Bf(3,nCent), xyz(3,nCent), dBf(3*nCent,3*nCent)
      Logical lWrite_, ldB, lWarn
      Character*8 Label
*
*     Call QEnter('ConInt')
*
      E1 = Work(ipEner +lIter-1)
      E0 = Work(ipEner0+lIter-1)
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
     &                       Label,' : Energy difference = ',
     &                       E0, ' hartree, ',
     &                       E0*CONV_AU_TO_KJ_PER_MOLE_,
     &                       ' kJ/mol'
            Write (6,'( A,F18.8,A)') '           Average energy    = ',
     &                               E1   , ' hartree'
         Else
            Write (6,'(2A,F18.8,A,F18.8,A)')
     &                       Label,' : Energy difference = ',
     &                       E1-E0, ' hartree, ',
     &                       (E1-E0)*CONV_AU_TO_KJ_PER_MOLE_,
     &                       ' kJ/mol'
            Write (6,'( A,F18.8,A)') '           E(i)              = ',
     &                               E1   , ' hartree'
            Write (6,'( A,F18.8,A)') '           E(j)              = ',
     &                               E0   , ' hartree'
         End If
      End If
*
*---- Compute the WDC B-matrix
*
      iOff = 0
      ipGrad1=(lIter-1)*3*nsAtom + ipGx
      ipGrad0=(lIter-1)*3*nsAtom + ipGx0
C     Call RecPrt('Grad1',' ',Work(ipGrad1),3,nCent)
C     Call RecPrt('Grad0',' ',Work(ipGrad0),3,nCent)
      Do iCent = 1, nCent
         Fact=DBLE(iDeg(xyz(1,iCent),iOper,nSym))
C        Write (6,*) 'Fact=',Fact
         Do iCar = 1, 3
            Bf(iCar,iCent)=Zero
            If (iOpt.eq.1) Then
C------------- Linear ------------------
               If (NADC) Then
                  Bf(iCar,iCent)=-Work(ipGrad0+iOff)
               Else
                  Bf(iCar,iCent)=-(Work(ipGrad1+iOff)
     &                            -Work(ipGrad0+iOff))
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
     &                             *Work(ipGrad0+iOff)
                  Else
                     Bf(iCar,iCent)=-Two*1.0D-5
     &                             *Work(ipGrad0+iOff)
                  End If
               Else
                  If (Abs(E1-E0).gt.1.0D-5) Then
                     Bf(iCar,iCent)=-Two*(E1-E0)
     &                             *(Work(ipGrad1+iOff)
     &                              -Work(ipGrad0+iOff))
                  Else
                     Bf(iCar,iCent)=-Two*1.0D-5
     &                             *(Work(ipGrad1+iOff)
     &                              -Work(ipGrad0+iOff))
                  End If
               End If
            Else If (iOpt.eq.3) Then
C------------- Absolute value ----------
               If (NADC) Then
                  Bf(iCar,iCent)=-Sign(One,E0)*(Work(ipGrad0+iOff))
               Else
                  Bf(iCar,iCent)=-Sign(One,E1-E0)
     &                          *(Work(ipGrad1+iOff)
     &                           -Work(ipGrad0+iOff))
               End If
            End If
*
            Bf(iCar,iCent)=Fact*Bf(iCar,iCent)
            iOff = iOff+1
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
            Do ix = 0, nGrad-1
               Do iy = 0, nGrad-1
                  If (NADC) Then
                     dBf(1+ix,1+iy)=-Two*
     &                               Work(ipGrad0+ix)
     &                             * Work(ipGrad0+iy)
                  Else
                     dBf(1+ix,1+iy)=-Two*
     &                              (Work(ipGrad1+ix)-
     &                               Work(ipGrad0+ix))
     &                             *(Work(ipGrad1+iy)-
     &                               Work(ipGrad0+iy))
                  End If
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
*     Call QExit('ConInt')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(lWarn)
      End
