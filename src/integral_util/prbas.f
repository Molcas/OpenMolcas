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
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      SubRoutine PrBas
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Authors: Roland Lindh, Martin Schuetz                            *
*              Dept. of Theoretical Chemistry,                         *
*              University of Lund, SWEDEN                              *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "angtp.fh"
#include "info.fh"
#include "WrkSpc.fh"
      Character*80 Lines(10)
*
*---- Statement Function
*
      IndSOff(iCnttp,iCnt)=(iCnttp-1)*Max_Cnt+iCnt
*
*     Call QEnter('PrBas',0)
*
*
*-----Loop over the basis sets and print all information
*
      mdc = 0
      iShell = 0
      jExp = 0
      iAOttp=0
*
*-----Loop over basis sets
*
      Do 10 iCnttp = 1, nCnttp
         Write (Lines(1),'(A,I3)') 'Unique basis set index:',iCnttp
         Write (Lines(2),*)
         Write (Lines(3),'(A,I3)')
     &         ' Number of different unique angular functions:',
     &         lOffAO(iCnttp)
         Call Banner(Lines,3,60)
         nTest = nVal_Shells(iCnttp)
         Write (6,*)
         Write (6,*) ' IndS'
         Write (6,*)
         iStr=iShell+1
         Do i = 1, dbsc(iCnttp)%nCntr
            iEnd = iStr + nTest-1
            Write (6,*) (IndS(j),j=iStr,iEnd)
            iStr=iEnd+1
         End Do
         Write (6,*)
         iShell = iShell + dbsc(iCnttp)%nCntr*nTest
         lComp = 0
         lSh = 0
*
*--------Loop over shells (s,p,d,..)
*
         Do 11 iAng = 0, nTest-1
            jSh = ipVal(iCnttp)+iAng
            nExpj=Shells(jSh)%nExp
            If (nBasis(jSh).eq.0) Go To 11
            Write (Lines(1),'(A,I3)') ' Shell of angular type ',iAng
            Write (Lines(2),*)
            Write (Lines(3),'(A,A,I3)')
     &         ' Number of different angular functions ',
     &         'preceding this shell',kOffAO(iCnttp,lSh)
            Call Banner(Lines,3,30)
            kCmp = (iAng+1)*(iAng+2)/2
            If (Prjct(jSh)) kCmp=2*iAng+1
*
*-----------Print exponents and contraction coefficients
*
            If (MaxPrm(iAng).gt.0 .and. nExpj.gt.0 .and.
     &             nBasis(jSh).gt.0) Then
               Write (6,*)
               If (Prjct(jSh).and.Transf(jSh)) Then
                  Write (6,*) ' Gaussian type: Spherical Harmonics'
               Else If (Transf(jSh)) Then
                  Write (6,*)
     &            ' Gaussian type: Spherical Harmonics and Contaminants'
               Else
                  Write (6,*) ' Gaussian type: Cartesians'
               End If
               Write (6,*)
               Write (6,*) '                 Type         '
               Write (6,'(19X,A)') AngTp(iAng)
               Write (6,*) '          No.      Exponent   ',
     &                     ' Contraction Coefficients'
            End If
            If (nBasis(jSh).gt.0) Then
               Do 13 kExp = 1, nExpj
                  jExp  = jExp  + 1
                  Write (6,'(10X,I3,1X,D16.9,10(1X,F9.4),'//
     &                     '3(/,30X,10(1X,F9.4)))')
     &                  jExp , Shells(jSh)%Exp(kExp),
     &                       ( Shells(jSh)%pCff(kExp,ib),
     &                  ib=1,nBasis(jSh))
 13            Continue
            End If
*
*-----------Print centers
*
            Write (6,*)
            Write (6,*) 'Charge:',Charge(iCnttp)
            Do 14 iCnt = 1, dbsc(iCnttp)%nCntr
               mShell = Ind_Shell(IndSOff(iCnttp,iCnt))+iAng+1
               Write (6,*)
               Write (6,*) 'iChCnt:',iChCnt(mdc+iCnt)
               Write (6,*)
               Write (6,*) ' CoSets'
               Do i = 0, 7
                  Write (6,'(8I1)') (iCoset(i,j,mdc+iCnt),j=0,7)
               End Do
               Write (6,*)
               Write (6,*) 'Stabilizer'
               Write (6,'(8I1)') (jStab(i,mdc+iCnt),
     &                            i=0,nStab(mdc+iCnt)-1)
               Write (6,*)
               Write (6,'(A,I3)') ' Unique shell index:', mShell
               Write (6,*)
               Write (6,'(A,I3)')
     &               ' Number of preseding unique angular functions',
     &               IndS(mShell)
               Write (6,*) ' Label   Cartesian Coordinates / Bohr'
               Write (6,*)
               Write (6,'(1X,A,1X,3F20.10)') LblCnt(mdc+iCnt),
     &                              (dbsc(iCnttp)%Coor(i,iCnt),i=1,3)
               Write (6,*)
*
*--------------Loop over the angular components of this unique shell
*
               iAO = iAOttp + (iCnt-1)*lOffAO(iCnttp)
     &             + kOffAO(iCnttp,iAng)
               Do 15 iCmp = 1, kCmp
                  Write (6,*)
                  Write (6,'(A,I3)')
     &              ' Unique angular index ', iAO+iCmp
                  Write (6,'(A,/,A,/,8(I4,1X))')
     &              ' Starting SO index for the first'//
     &              ' contracted basis function of this angular type',
     &              ' in each irrep',(iAOtSO(iAO+iCmp,iIrrep),
     &              iIrrep=0,nIrrep-1)
                  Write (6,'(A,I3)')
     &                     ' Symmetry label:',IrrCmp(IndS(mShell)+iCmp)
 15            Continue
 14         Continue
*
            lSh = lSh + 1
 11      Continue
         mdc = mdc + dbsc(iCnttp)%nCntr
         iAOttp = iAOttp + lOffAO(iCnttp)*dbsc(iCnttp)%nCntr
 10   Continue
      Write (6,*)
      Write (6,'(A42,8I4)') 'Number of basis functions in each irrep:',
     &                      (iIrrep,iIrrep=0,nIrrep-1)
      Write (6,'(42X,8I4)') (nBas(iIrrep),iIrrep=0,nIrrep-1)
      Write (6,*)
*
*
*     Call QExit('PrBas',0)
      Return
      End
