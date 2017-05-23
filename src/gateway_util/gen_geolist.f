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
      Subroutine Gen_GeoList(DInf,nDInf)
      use GeoList
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 DInf(nDInf)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Centr,3,mCentr,label='Centr')
      Call mma_allocate(Mass,mCentr,label='Mass')
      Call mma_allocate(Chrg,mCentr,label='Chrg')
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the center list.
*
      kCentr=0
*
      nc = 1
      Do jCnttp = 1, nCnttp
         mCnt = nCntr(jCnttp)
*
*        Do not include Auxiliary basis sets, or fragment basis sets
*
         If (AuxCnttp(jCnttp).or.FragCnttp(jCnttp)) Go To 1212
*
*        Do not include ECP basis sets which doesn't have any valence
*        basis set.
*
         If (ECP(jCnttp).and.nVal_Shells(jCnttp).eq.0) Go To 1212
*
         jxyz = ipCntr(jCnttp)
         Do jCnt = 1, mCnt
            ndc = jCnt + mdciCnttp(jCnttp)
            x1 = DInf(jxyz)
            y1 = DInf(jxyz+1)
            z1 = DInf(jxyz+2)
            Do i = 0, nIrrep/nStab(ndc) - 1
               iFacx=iPhase(1,iCoset(i,0,ndc))
               iFacy=iPhase(2,iCoset(i,0,ndc))
               iFacz=iPhase(3,iCoset(i,0,ndc))
               Centr(1,nc) = x1*DBLE(iFacx)
               Centr(2,nc) = y1*DBLE(iFacy)
               Centr(3,nc) = z1*DBLE(iFacz)
               nchr=iAtmNr(jCnttp)
               If (nchr.ge.0) Then
                  Mass(nc) = rMass(nchr)
               Else
                  Mass(nc) = Zero
               End If
               nchr=iAtmNr(jCnttp)
               If (nchr.ge.0) Then
                  Chrg(nc) = DBLE(nchr)
               Else
                  Chrg(nc) = Zero
               End If
               if (nc.gt.8*mxdc) Then
                  Call WarningMessage(2,'lblxxx too small')
                  Call Abend()
               End If
               lblxxx(nc)=lblcnt(ndc)(1:LENIN)
               nc = nc + 1
            End Do
            jxyz = jxyz + 3
            kCentr = kCentr + nIrrep/nStab(ndc)
         End Do
 1212    Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Compute Total Charge and Center of Charge centroid
*
      Call CoW(Centr,CoC,Chrg,kCentr,qNuc)
      If (iChCar(1).ne.0) CoC(1)=Zero
      If (iChCar(2).ne.0) CoC(2)=Zero
      If (iChCar(3).ne.0) CoC(3)=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     Modify charges to effective charges.
*
      nc = 1
      Do jCnttp = 1, nCnttp
         Z = Charge(jCnttp)
         mCnt = nCntr(jCnttp)
         If (AuxCnttp(jCnttp).or.FragCnttp(jCnttp)) Go To 1213
         jxyz = ipCntr(jCnttp)
         Do jCnt = 1, mCnt
            ndc = jCnt + mdciCnttp(jCnttp)
            Do i = 0, nIrrep/nStab(ndc) - 1
               nchr=iAtmNr(jCnttp)
               Chrg(nc) = Z
               nc = nc + 1
            End Do
            jxyz = jxyz + 3
         End Do
 1213    Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Compute Total Mass and Center of Mass
*
      Call CoW(Centr,CoM,Mass,kCentr,TMass)
      If (iChCar(1).ne.0) CoM(1)=Zero
      If (iChCar(2).ne.0) CoM(2)=Zero
      If (iChCar(3).ne.0) CoM(3)=Zero
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
