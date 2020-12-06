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
      Subroutine SphInt(xyz,nCent,iOfRef,RR0,Bf,l_Write,Label,dBf,ldB)
      use Slapaf_Info, only: Weights, RefGeo
      Implicit Real*8  (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "weighting.fh"
#include "info_slapaf.fh"
      Real*8   Bf(3,nCent), xyz(3,nCent), dBf(3,nCent,3,nCent)
      Logical l_Write, ldB
      Character*8 Label
*
      xyz0(i,j)=Work(ipRef_+(j-1)*3+i-1)
*
*
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the radius of the hypersphere
*
      If (iOfRef.eq.ip_Dummy) Then
         ipRef_=ip_of_Work(RefGeo(1,1))
      Else
         ipRef_=iOfRef
      End If
C     Call RecPrt('SphInt: xyz',' ',xyz,3,nCent)
C     Call RecPrt('Ref: xyz0',' ',Work(ipRef_),3,nCent)
      RR0=Zero
      TWeight=Zero
      Do iCent = 1, nCent
         Fact=DBLE(iDeg(xyz(1,iCent)))
         xWeight=Fact*Weights(iCent)
         TWeight=TWeight+xWeight
C        Write (*,*) 'xWeight=',xWeight
         Do ixyz = 1, 3
C           Write (*,*)xyz(ixyz,iCent),xyz0(ixyz,iCent)
            temp=xyz(ixyz,iCent)-xyz0(ixyz,iCent)
            RR0=RR0 + xWeight*temp**2
         End Do
      End Do
      RR0_unscaled=Sqrt(RR0)
*
*     RR0_unscaled is the real (weighthed) distance,
*     RR0 is scaled by 1/Sqrt(TWeight) and so are the derivatives
*
      SqInvTWeight=One/Sqrt(TWeight)
      RR0=RR0_unscaled*SqInvTWeight
*
      If (l_Write) Then
         Write (6,'(2A,F18.8,A)') Label,' : Radius of h-sphere= ',
     &                            RR0,
     &                            ' au (weighted/sqrt(total weight))'
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the WDC B-matrix
*
*FIXME: revise the symmetry
      Do iCent = 1, nCent
         Fact=DBLE(iDeg(xyz(1,iCent)))
         xWeight=Fact*Weights(iCent)
         Do iCar = 1, 3
            temp=xyz(iCar,iCent)-xyz0(iCar,iCent)
            If (RR0_unscaled.ne.Zero) Then
               temp=xyz(iCar,iCent)-xyz0(iCar,iCent)
               Bf(iCar,iCent)=xWeight*temp/RR0_unscaled*SqInvTWeight
            Else
*
*              If we are standing on the reference point the gradient
*              is not well defined.
*
               Bf(iCar,iCent)=Zero
            End If
         End Do
      End Do
c     Call RecPrt('Bf',' ',Bf,3,nCent)
*                                                                      *
************************************************************************
*                                                                      *
*
*---- Compute the cartesian derivative of the B-Matrix.
*
*FIXME: revise the symmetry
      If (ldB) Then
         Call FZero(dBf,(3*nCent)**2)
         If (RR0.eq.Zero) Go To 99
         Do iCent = 1, nCent
            Fact=DBLE(iDeg(xyz(1,iCent)))
            xWeight=Fact*Weights(iCent)
            Do ixyz = 1, 3
               tempi=xyz(ixyz,iCent)-xyz0(ixyz,iCent)
               Do jCent = 1, nCent
                  Fact=DBLE(iDeg(xyz(1,jCent)))
                  yWeight=Fact*Weights(jCent)
                  Do jxyz = 1, 3
                     tempj=xyz(jxyz,jCent)-xyz0(jxyz,jCent)
                     temp=Zero
                     If (ixyz.eq.jxyz.and.iCent.eq.jCent)
     &                  temp=RR0_unscaled
                     temp=temp-yWeight*tempi*tempj/RR0_unscaled
                     temp=(xWeight*temp)/RR0_unscaled**2
                     dBf(ixyz,iCent,jxyz,jCent) = temp*SqInvTWeight
                  End Do
               End Do
            End Do
         End Do
 99      Continue
C        Call RecPrt('dBf',' ',dBf,3*nCent,3*nCent)
*
      End If
*
      Return
      End
