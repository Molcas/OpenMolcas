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
      Subroutine Min_Mult_Error(EC,A,B,Ci,Cj,rMP,xrMP,xxrMP,xnrMP,lMax,
     &                          nij,nElem,iAtom,jAtom,nAtoms,nPert,
     &                          C_o_C,Scratch_New,Scratch_Org,iPlot,
     &                          T_Values,iWarnings,Num_Warnings)
      Implicit Real*8 (A-H, O-Z)
#include "real.fh"
      Dimension EC(3,nij),C_o_C(3),A(3,nij),B(3,nij),Ci(3),Cj(3),R_ij(3)
      Dimension rMP(nij,0:nElem-1,0:nPert-1),xnrMP(nij,nElem)
      Dimension xrMP(nij,nElem),xxrMP(nij,nElem),T_Values(nij)
      Dimension Scratch_Org(nij*(2+lMax+1)),Scratch_New(nij*(2+lMax+1))
      Dimension iWarnings(nij)
      Parameter (Error_Threshold = 1.0D-12, Delta_Threshold = 1.0D-12)
      External Error_for_t,Golden
*
      iDim = nij*nElem
      ij = iAtom*(iAtom-1)/2+jAtom
      iPrint_Errors   = 0
      l = lMax - 1
      Do i = 1,3
         R_ij(i) = Cj(i) - Ci(i)
      End Do
*
      t_min       = Zero
      t_max       = Zero
      Do i = 1,3
         If (R_ij(i) .ne. 0.0D0) Then
            t_min       = (Ci(i)-EC(i,ij))/R_ij(i)
            t_max       = (Cj(i)-EC(i,ij))/R_ij(i)
         End If
      End Do
      Delta_Orig  = 0.1D0
      Delta       = Delta_Orig
*
* Check the range of possible t-values for minima
*
      num_min     = 0
      iSlope      = 0
      If (iPlot .eq. 1) Then
         Write(6,*)
         Write(6,*) 'iAtom, jAtom = ',iAtom,jAtom
      End If
      t_temp = t_min
      Error_Best = -1.0D0
      Error_Old  = Zero
      t_best = 0.0D0
  50  Continue
         Error = Error_for_t(t_temp,rMP,xrMP,xxrMP,xnrMP,EC,A,R_ij,
     &                       C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,
     &                       Scratch_New,Scratch_Org,iPrint_Errors)
         If (iPlot .eq. 1) Then
            Write(6,'(1X,A,F5.2,F16.12)') 't, Error = ',t_temp,Error
            Call xFlush(6)
         End If
         Delta_Error = Error - Error_Old
         Error_Old   = Error
         iSlope_Old  = iSlope
         If (Abs(Delta_Error) .lt. Error_Threshold) Then
            iSlope = 0
         Else If (Delta_Error .lt. Zero) Then
            iSlope = -1
         Else
            iSlope = 1
         End If
         If (iSlope_Old .lt. 0 .And. iSlope .ge. 0) Then
            num_min = num_min + 1
         End If
         If (Error .lt. Error_Best .OR. Error_Best .lt. Zero) Then
            Error_Best = Error
            t_best     = t_temp
         End If
         t_temp = t_temp + Delta/10.0D0
      If (t_temp .le. t_max+Delta/100.0D0) Goto 50
*
* Any warnings from scan?
*
      If (num_min .gt. 1) Then
         iWarnings(ij) = 1
         Num_Warnings  = Num_Warnings + 1
      End If
*
* Find minima with Golden Section Search
*
* (It's assumed that the error function is either constant or has at
* least one minima somewhere, although not necessarily within the
* allowed range).
*
      ax = t_best
      bx = t_best + Delta
* First make an initial bracketing of the minima
      Call mnBrak(ax,bx,cx,fa,fb,fc,Error_for_t,rMP,xrMP,xxrMP,
     &           xnrMP,EC,A,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,
     &           Scratch_New,Scratch_Org,iPrint_Errors)
*
      If (abs(fa-fc) .lt. Error_Threshold) Then
         iWarnings(ij) = 4
         Num_Warnings = Num_Warnings + 1
         t_final = Zero
      Else
         Error = Golden(ax,bx,cx,Error_for_t,Delta_Threshold,
     &           Error_Threshold,t_final,rMP,xrMP,xxrMP,
     &           xnrMP,EC,A,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,
     &           Scratch_New,Scratch_Org,iPrint_Errors)
      End If
*
* Check that the minima is within the allowed range
*
      If (t_final .gt. t_max) Then
         t_final = t_max
         iWarnings(ij) = 2
         Num_Warnings  = Num_Warnings + 1
      Else If (t_final .lt. t_min) Then
         t_final = t_min
         iWarnings(ij) = 2
         Num_Warnings  = Num_Warnings + 1
      End If
*
      T_values(ij) = t_final
      B(1,ij) = EC(1,ij) + t_final*R_ij(1)
      B(2,ij) = EC(2,ij) + t_final*R_ij(2)
      B(3,ij) = EC(3,ij) + t_final*R_ij(3)
*
      Return
      End
