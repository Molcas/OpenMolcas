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
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine EFNuc(CoOP,Chrg,Coor,nAtm,ESIT,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute the electricstatic interaction tensor contribution*
!         from the nuclei. In the case that the test charge coincide   *
!         with a nucleau we will remove that center.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, April '95.                                      *
!***********************************************************************
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nAtm, nOrdOp
      Real*8 Chrg(nAtm), Coor(3,nAtm), ESIT((nOrdOp+1)*(nOrdOp+2)/2)

      Integer, Allocatable:: C_ESIT(:)
      Real*8 CoOp(3)
      Integer nTot, iPowr, iAtom, ix, iy, iz
      Real*8 Fact, x, y, z, r2, Thr, r, eix, eiy, eiz, temp
#ifdef _DEBUGPRINT_
      Integer n, nElem
!
!---- Statement function
!
      nElem(n)=(n+1)*(n+2)/2
#endif
!
!     Compute the nuclear contribution to the electrostatic interation
!     tensor, ESIT.
!
      ESIT(:)=Zero
!
      nTot=(nOrdOp+1)**6
      Call mma_allocate(C_ESIT,nTot,Label='ESIT')
      Call InitIA(C_ESIT,nOrdOp)
!
      iPowR=2*nOrdOp+1
      Fact=One
      If (nOrdOp.ge.1) Fact=-One
      Do iAtom = 1, nAtm
         x = CoOp(1) - Coor(1,iAtom)
         y = CoOp(2) - Coor(2,iAtom)
         z = CoOp(3) - Coor(3,iAtom)
         r2 = x**2 + y**2 + z**2
         Thr=1.0D-12
         If (r2.gt.Thr) Then
            r  = Chrg(iAtom)/Sqrt(r2)**iPowR
            Do ix = nOrdOp, 0, -1
               Do iy = nOrdOp-ix, 0, -1
                  iz = nOrdOp - ix - iy
                  If (ix.eq.0) Then
                     EIx=One
                  Else
                     EIx=x**ix
                  End If
                  If (iy.eq.0) Then
                     EIy=One
                  Else
                     EIy=y**iy
                  End If
                  If (iz.eq.0) Then
                     EIz=One
                  Else
                     EIz=z**iz
                  End If
                  temp=Fact*EIx*EIy*EIz*r
!
                  Call ContEI(C_ESIT,nOrdOp,ESIT,ix,iy,iz,temp)
!
               End Do
            End Do       ! End loop over cartesian combinations
         End If
      End Do             ! End loop over atoms
!
      Call mma_deallocate(C_ESIT)
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' The Electrostatic Interaction'                      &
     &                 //' Tensor',' ',ESIT,nElem(nOrdOp),1)
#endif
      Return
      End SubRoutine EFNuc
