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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "weighting.fh"
#include "info_slapaf.fh"
      Real*8 Bf(3,nCent), xyz(3,nCent), dBf(3,nCent,3,nCent)

      Logical l_Write, ldB, lTrans
      Character*8 Label
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      xyz0(i,j)=Work(ipRef+(j-1)*3+i-1)
      r12(i,j)=Work(ipR12_+(j-1)*3+i-1)
*                                                                      *
************************************************************************
*                                                                      *
*     The reference direction (normal to the hyperplane) is taken from
*     the input (ipGradRef) if Ref_Grad=.True., or from the stored
*     transverse direction if there is one, or from the R-P vector
*     (ipR12) otherwise
*
      Call qpg_dArray('Transverse',lTrans,nTrans)
      If (Ref_Grad) Then
         lTrans = .False.
         ipR12_ = ipGradRef
c        write(6,*), 'Using Reference Gradient'
      Else If (lTrans) Then
         Call Allocate_Work(ipR12_,3*nCent)
         Call Get_dArray('Transverse',Work(ipR12_),3*nCent)
c        write(6,*), 'Using stored Transverse'
      Else
         ipR12_ = ipR12
c        write(6,*), 'Using R-P Vector'
      End If
*                                                                      *
************************************************************************
*                                                                      *
c     Call RecPrt('Ref',' ',Work(ipRef),3,nCent)
c     Call RecPrt('R12',' ',Work(ipR12_),3,nCent)
*
*     Length of the direction vector (weighted)
*
      RR_R12 = Zero
      TWeight = Zero
      Do iCent = 1, nCent
         Fact = Dble(iDeg(xyz(1,iCent)))
         xWeight = Fact*Work(ipWeights+iCent-1)
         TWeight = TWeight+xWeight
         Do i = 1, 3
            RR_R12 = RR_R12 + xWeight*(r12(i,iCent))**2
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
         xWeight = Fact*Work(ipWeights+iCent-1)
         Do i = 1, 3
            f = f + xWeight*(xyz(i,iCent)-xyz0(i,iCent))*r12(i,iCent)
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
            xWeight = Fact*Work(ipWeights+iCent-1)
            Do i = 1, 3
               Bf(i,iCent) = xWeight*r12(i,iCent)/RR_R12*SqInvTWeight
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
      If (lTrans) Call Free_Work(ipR12_)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
