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
! Copyright (C) 1996, Anders Bernhardsson                              *
!               2006, Giovanni Ghigo                                   *
!***********************************************************************
      SubRoutine TrGrd_Alaska(CGrad,CNames,GradIn,nGrad,iCen)
!***********************************************************************
!                                                                      *
!      Transforms a symmetry adapted gradient to unsymmetric form.     *
!                                                                      *
!       Written by Anders Bernhardsson                                 *
!       Adapted by Giovanni Ghigo                                      *
!       University of Torino, July 2006                                *
!***********************************************************************
      Use Basis_Info, only: nCnttp, DBSC
      Use Center_Info, only: DC
      use Symmetry_Info, only: nIrrep
      use Constants, only: Zero
      use Disp, only: IndDsp
      Implicit None
#include "Molcas.fh"
#include "SysDef.fh"
      Integer nGrad
      Real*8 CGrad(3,MxAtom)
      Real*8 GradIn(nGrad)
      Character(LEN=LENIN5)  CNames(MxAtom)
      Integer iCen

      Real*8, Parameter :: tol=1d-8
      Logical, External :: TF
      Integer mdc, iIrrep, iCnttp, iCnt, iCo, kOp, nDisps, iCar,
     &        iComp
      Integer, External :: iPrmt, NrOpr
      Real*8 Xr
!
      mdc=0
      iIrrep=0
!
      call dcopy_(3*MxAtom,[Zero],0,CGrad,1)
      iCen=0
!     nCnttp_Valence=0
!     Do iCnttp = 1, nCnttp
!        If (dbsc(iCnttp)%Aux) Go To 999
!        nCnttp_Valence = nCnttp_Valence+1
!     End Do
!999  Continue
!
      Do iCnttp=1,nCnttp
         If (.Not.(dbsc(iCnttp)%pChrg.or.
     &             dbsc(iCnttp)%Frag.or.
     &             dbsc(iCnttp)%Aux)) Then
            Do iCnt=1,dbsc(iCnttp)%nCntr
               mdc=mdc+1
               Do iCo=0,nIrrep/dc(mdc)%nStab-1
                  kop=dc(mdc)%iCoSet(iCo,0)
                  nDispS = IndDsp(mdc,iIrrep)
                  iCen=iCen+1
                  Do iCar=0,2
                     iComp = 2**iCar
                     If ( TF(mdc,iIrrep,iComp)) Then
                        nDispS = nDispS + 1
                        XR=DBLE(iPrmt(NrOpr(kop),icomp))
                        CGrad(iCar+1,iCen)=XR*GradIn(nDispS)
                     End If
                  End Do
                  CNames(iCen)=dc(mdc)%LblCnt
               End Do
            End Do
         End If
      End Do
!
      End SubRoutine TrGrd_Alaska
