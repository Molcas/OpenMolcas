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
      Subroutine Transverse(xyz,nCent,HDist,Bf,l_write,Label,dBf,ldB)
      use Slapaf_Info, only: Weights, RefGeo, R12, GradRef
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "weighting.fh"
      Real*8 Bf(3,nCent), xyz(3,nCent), dBf(3,nCent,3,nCent)
      Logical l_Write, ldB, lTrans
      Character(LEN=8) Label
      Real*8, Allocatable, Target:: TV(:,:)
      Real*8, Pointer:: r12_p(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*     The reference direction (normal to the hyperplane) is taken from
*     the input (GradRef) if allocated, or from the stored
*     transverse direction if there is one, or from the R-P vector
*     (R12) otherwise
*
      If (Allocated(GradRef)) Then
         lTrans = .False.
         r12_p => GradRef
c        write(6,*), 'Using Reference Gradient'
      Else
         Call qpg_dArray('Transverse',lTrans,nTrans)
         If (lTrans) Then
            Call mma_allocate(TV,3,nCent,Label='TV')
            Call Get_dArray('Transverse',TV,3*nCent)
            r12_p => TV
c           write(6,*), 'Using stored Transverse'
         Else
            r12_p => R12
c           write(6,*), 'Using R-P Vector'
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
c     Call RecPrt('Ref',' ',RefGeo,3,nCent)
c     Call RecPrt('R12',' ',r12_p,3,nCent)
*
*     Length of the direction vector (weighted)
*
      RR_R12 = Zero
      TWeight = Zero
      Do iCent = 1, nCent
         Fact = Dble(iDeg(xyz(1,iCent)))
         xWeight = Fact*Weights(iCent)
         TWeight = TWeight+xWeight
         Do i = 1, 3
            RR_R12 = RR_R12 + xWeight*(r12_p(i,iCent))**2
         End Do
      End Do
      RR_R12=Sqrt(RR_R12)
*
*     The distance scaling will be 1/Sqrt(TWeight)
*
      SqInvTWeight=One/Sqrt(TWeight)
*
*     Dot product between the x-x0 vector and the direction (weighted)
*
      f=Zero
      Do iCent = 1, nCent
         Fact = Dble(iDeg(xyz(1,iCent)))
         xWeight = Fact*Weights(iCent)
         Do i = 1, 3
            f = f + xWeight*(xyz(i,iCent)-RefGeo(i,iCent))
     &        *  r12_p(i,iCent)
         End Do
      End Do
*
*     The distance to the plane is the dot product / direction length
*
      If (RR_R12.eq.Zero) Then
         HDist = Zero
      Else
         HDist = f / RR_R12 * SqInvTWeight
      End If
c     Write (6,*) 'f, RR_R12=',f,RR_R12
*
      If (l_Write) Then
         Write (6,'(2A,F18.8,A)') Label,' : Hyperplane distance =',
     &                            HDist,
     &                            ' au (weighted/sqrt(total weight)'
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the WDC B-matrix
*     If the direction is null, the derivative is not defined
*
      Call FZero(Bf,3*nCent)
      If (RR_R12.gt.Zero) Then
*
*        The derivative is simply the unit direction vector
*        (with weighting and scaling accounted for)
*
         Do iCent = 1, nCent
            Fact = Dble(iDeg(xyz(1,iCent)))
            xWeight = Fact*Weights(iCent)
            Do i = 1, 3
               Bf(i,iCent) = xWeight*r12_p(i,iCent)/RR_R12*SqInvTWeight
            End Do
         End Do
      End If
c     Call RecPrt('Bf',' ',Bf,3,nCent)
*                                                                      *
************************************************************************
*                                                                      *
*     The second derivative is null, as the derivative is constant
*
      If (ldB) Then
         Call FZero(dBf,(3*nCent)**2)
      End If
*
      If (lTrans) Call mma_deallocate(TV)
      r12_p => Null()
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
