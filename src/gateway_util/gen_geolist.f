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
      Subroutine Gen_GeoList()
      use GeoList
      use Basis_Info
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
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
         mCnt = dbsc(jCnttp)%nCntr
*
*        Do not include Auxiliary basis sets, or fragment basis sets
*
         If (dbsc(jCnttp)%Aux.or.dbsc(jCnttp)%Frag) Cycle
*
*        Do not include ECP basis sets which does not have any valence
*        basis set.
*
         If (dbsc(jCnttp)%ECP.and.dbsc(jCnttp)%nVal.eq.0) Cycle
*
         Do jCnt = 1, mCnt
            ndc = jCnt + dbsc(jCnttp)%mdci
            x1 = dbsc(jCnttp)%Coor(1,jCnt)
            y1 = dbsc(jCnttp)%Coor(2,jCnt)
            z1 = dbsc(jCnttp)%Coor(3,jCnt)
            Do i = 0, nIrrep/nStab(ndc) - 1
               iFacx=iPhase(1,iCoset(i,0,ndc))
               iFacy=iPhase(2,iCoset(i,0,ndc))
               iFacz=iPhase(3,iCoset(i,0,ndc))
               Centr(1,nc) = x1*DBLE(iFacx)
               Centr(2,nc) = y1*DBLE(iFacy)
               Centr(3,nc) = z1*DBLE(iFacz)
               nchr=dbsc(jCnttp)%AtmNr
               If (nchr.ge.0) Then
                  Mass(nc) = dbsc(jCnttp)%CntMass
               Else
                  Mass(nc) = Zero
               End If
               nchr=dbsc(jCnttp)%AtmNr
               If (nchr.ge.0) Then
                  Chrg(nc) = DBLE(nchr)
               Else
                  Chrg(nc) = Zero
               End If
               nc = nc + 1
            End Do
            kCentr = kCentr + nIrrep/nStab(ndc)
         End Do
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
         Z = dbsc(jCnttp)%Charge
         mCnt = dbsc(jCnttp)%nCntr
         If (dbsc(jCnttp)%Aux.or.dbsc(jCnttp)%Frag) Cycle
         Do jCnt = 1, mCnt
            ndc = jCnt + dbsc(jCnttp)%mdci
            Do i = 0, nIrrep/nStab(ndc) - 1
               nchr=dbsc(jCnttp)%AtmNr
               Chrg(nc) = Z
               nc = nc + 1
            End Do
         End Do
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
