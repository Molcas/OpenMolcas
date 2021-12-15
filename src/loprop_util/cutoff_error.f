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
      Subroutine CutOff_Error(l,lMax,rMP,xrMP,nij,EC,C_o_C,nElem,
     &                        Scratch_New,Scratch_Org,nAtoms,iPrint,
     &                        Cut_Off_Error)
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "status.fh"
      Real*8 rMP(nij,nElem),xrMP(nij,nElem),EC(3,nij),C_o_C(3)
      Real*8 Scratch_New(nij*(2+lMax+1)),Scratch_Org(nij*(2+lMax+1))
      Character*80 Banner_Line
*     Statement function
*
      mElem(i)=(i+1)*(i+2)*(i+3)/6
*
      iEnd  = mElem(lMax)
      iStrt = mElem(l) + 1
      ij = 0
      Do iAtom = 1, nAtoms
         Do jAtom = 1, iAtom
            ij = ij + 1
            Call ReExpand(xrMP,nij,nElem,C_o_C,EC(1,ij),ij,lMax)
*
*--- Set all elements corresponding to l+1 = 0
            Do iElem = iStrt, iEnd
               xrMP(ij,iElem) = Zero
            End Do
*
            Call ReExpand(xrMP,nij,nElem,EC(1,ij),C_o_C,ij,lMax)
         End Do
      End Do
*
      If (iPrint .ge. 1) Then
         Write (6,*)
         Write(Banner_Line,'(A,I2)') 'Errors introduced by zeroing '
     &      // 'multipole moments greater than l = ',l
         Call Banner (Banner_Line,1,80)
      End If
      Sum    = 0.0D0
      iElem  = mElem(l) + 1
      Do k = l+1, lMax
         If (iPrint .ge. 1) Then
            Write (6,*)
            Write (6,'(A,I1)') 'l=',k
            Write (6,*)
            Write (6,*)
     &   ' m     Original       New            Error            Percent'
            Write (6,*)
         End If
*
         kDim = (k+1)*(k+2)/2
         Call DGEMM_('N','N',
     &               nij,2*k+1,kDim,
     &               1.0d0,xrMP(1,iElem),nij,
     &               RSph(ipSph(k)),kDim,
     &               0.0d0,Scratch_New,nij)
         Call DGEMM_('N','N',
     &               nij,2*k+1,kDim,
     &               1.0d0,rMP(1,iElem),nij,
     &               RSph(ipSph(k)),kDim,
     &               0.0d0,Scratch_Org,nij)


         iOff   = 1
         rms    = 0.0D0
         iCount = 0
         Do m = -k, k
            Original  = DDot_(nij,[One],0,Scratch_Org(iOff),1)
            Estimated = DDot_(nij,[One],0,Scratch_New(iOff),1)
            Error     = Original - Estimated

            Sum       = Sum + Error*Error
            rms       = rms + Error*Error
            Percent   = 0.0D0
            If (abs(Error) .lt. 1.0D-8) Then
               Percent   = 0.0D0
            Else If (abs(Original) .gt. 1.0D-13) Then
               Percent   = abs(Error/Original) * 100.0d0
            Else
               Percent   = -999.99D0
            End If

            If (iPrint .ge. 1) Then
               If (Percent .ge. Zero) Then
                  Write (6,'(I3,3F16.8,4X,F6.2)')
     &                         m, Original, Estimated, Error, Percent
               Else
                  Write (6,'(I3,3F16.8,4X,A)')
     &                         m, Original, Estimated, Error, 'Infinite'
               Endif
            End If

            iOff = iOff + nij
         End Do
         If (iPrint .ge. 1) Then
            rms = sqrt(rms/DBLE(2*k+1))
            Write(6,*)
            Write(6,'(A,F16.8)') 'Root mean square = ',rms
         End If
*
         iElem = iElem + (k+1)*(k+2)/2
      End Do
*
      Cut_Off_Error = Sum
      Return
      End
