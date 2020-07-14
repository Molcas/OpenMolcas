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
      SubRoutine Dstncs(Lbls,xyz,mCentr,Angstr,Max_Center,iCols)
************************************************************************
*                                                                      *
* Object: to compute distances from a coordinate list                  *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
      Real*8 xyz(3,mCentr)
      Character*(LENIN) Lbls(mCentr)
      Real*8, Allocatable:: BST(:)
      Integer, Allocatable:: iBST(:,:)
*
      lu=6
      If (mCentr.le.Max_Center) Then
*
         Do i = 1, 2
            Write (Lu,*)
            If (i.eq.1) Then
               Fact = One
               Write (Lu,'(19X,A)')
     &            ' *************************************** '
               Write (Lu,'(19X,A)')
     &            ' *    InterNuclear Distances / Bohr    * '
               Write (Lu,'(19X,A)')
     &            ' *************************************** '
            Else
               Fact = Angstr
               Write (Lu,'(19X,A)')
     &            ' ******************************************* '
               Write (Lu,'(19X,A)')
     &            ' *    InterNuclear Distances / Angstrom    * '
               Write (Lu,'(19X,A)')
     &            ' ******************************************* '
            End If
            Do icc = 1, mCentr, iCols
               Write (Lu,*)
               If (iCols.EQ.6) Write (Lu,'( 9X,6(5X,I2,1X,A,2X))')
     &               (ic,Lbls(ic),ic=icc,Min(icc+5,mCentr))
               If (iCols.EQ.5) Write (Lu,'( 9X,5(5X,I2,1X,A,2X))')
     &               (ic,Lbls(ic),ic=icc,Min(icc+4,mCentr))

               Do jc = icc, mCentr
                  x1 = xyz(1,jc)
                  y1 = xyz(2,jc)
                  z1 = xyz(3,jc)
                  If (iCols.EQ.6) Write (Lu,'(I5,1X,A,1X,6(F10.6,6X))')
     &                  jc,Lbls(jc),
     &                  ( Fact * Sqrt( (xyz(1,ic)-x1)**2 +
     &                                 (xyz(2,ic)-y1)**2 +
     &                                 (xyz(3,ic)-z1)**2 ),
     &                  ic = icc, Min(jc,icc+5,mCentr))
                  If (iCols.EQ.5) Write (Lu,'(I5,1X,A,1X,5(F10.6,6X))')
     &                  jc,Lbls(jc),
     &                  ( Fact * Sqrt( (xyz(1,ic)-x1)**2 +
     &                                 (xyz(2,ic)-y1)**2 +
     &                                 (xyz(3,ic)-z1)**2 ),
     &                  ic = icc, Min(jc,icc+4,mCentr))
              End Do
            End Do
         End Do
         Return
      Else
*
         Write (Lu,*)
         Write (Lu,'(19X,A)')
     &            ' ************************************************* '
         Write (Lu,'(19X,A)')
     &            ' **** InterNuclear Distances / Bohr, Angstrom **** '
         Write (Lu,'(19X,A)')
     &            ' ************************************************* '
         Write (Lu,*)
         Write (Lu,'(A)') '     Atom centers     '//
     &                    '    Bohr        Angstrom'
*
CVV   Set .false. to get faster printing without sorting.
*
         If (.true.) then
            Thr_R=(3.0D0/Angstr)**2
            Thr_D=1D-4
            Call mma_allocate(BST,mCentr**2,Label='BST')
            Call mma_allocate(iBST,2,mCentr**2,Label='iBST')
            iiBST=0
            Do icc = 1, mCentr
               x1 = xyz(1,icc)
               y1 = xyz(2,icc)
               z1 = xyz(3,icc)
               Do jcc = 1, icc-1
                  x2 = xyz(1,jcc)
                  y2 = xyz(2,jcc)
                  z2 = xyz(3,jcc)
                  R = (x2-x1)**2+(y2-y1)**2+(z2-z1)**2
                  If (R.le.Thr_R) then
                     iiBST=iiBST+1
                     BST(iiBST)=R
                     iBST(1,iiBST)=icc
                     iBST(2,iiBST)=jcc
                  End If
               End Do
            End Do
*#define _DEBUG_
#ifdef _DEBUG_
            do ii=1,iiBST
                R=sqrt(BST(ii))
                Write (Lu,'(2(I5,1X,A4),2(F10.6,6X))')
     &               iBST(1,ii),Lbls(iBST(1,ii)),
     &               iBST(2,ii),Lbls(iBST(2,ii)),
     &               R,R*Angstr
            end do
#endif
            MoreToGo=1
            Do While (MoreToGo.eq.1)
*
*              Find the shortest distance between any atoms
*
               R=100D0
               do ii=1,iiBST
                  R=MIN(R,BST(ii))
               enddo
*
               If (R.gt.90D0) Exit
*
               moretogo=0
               isfirst=1
               do ii=1,iiBST
*
                  If (abs(R-BST(ii)).lt.Thr_D) Then
                     If (isfirst.eq.1) Then
                        RR=SQRT(R)
                        Write (Lu,'(2(I5,1X,A),2(F10.6,6X))')
     &                         iBST(1,ii),Lbls(iBST(1,ii)),
     &                         iBST(2,ii),Lbls(iBST(2,ii)),
     &                         RR,RR*Angstr
                        isfirst=0
                     Else
                        Write (Lu,'(2(I5,1X,A),2(F10.6,6X))')
     &                  iBST(1,ii),Lbls(iBST(1,ii)),
     &                  iBST(2,ii),Lbls(iBST(2,ii))
                     End If
                     BST(ii)=100D0 ! Effectively remove from the list
                     moretogo=1
                   End If
               End Do
            End Do
            Call mma_deallocate(iBST)
            Call mma_deallocate(BST)
*
         Else
*
            Thr_R=3.0D0
            Do icc = 1, mCentr
               x1 = xyz(1,icc)
               y1 = xyz(2,icc)
               z1 = xyz(3,icc)
               Do jcc = 1, icc-1
                  x2 = xyz(1,jcc)
                  y2 = xyz(2,jcc)
                  z2 = xyz(3,jcc)
                  R = Sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                  If (R*Angstr.le.Thr_R)
     &               Write (Lu,'(2(I5,1X,A),2(F10.6,6X))')
     &               icc,Lbls(icc), jcc,Lbls(jcc), R, R*Angstr
               End Do
            End Do
          End If
      End If
*
      Return
      End
