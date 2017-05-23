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
      Subroutine Fold_tMat(nSym,nBas,A,B)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Fold up symmetry blocked matrix A (UT-storage)                   *
*     into triangular matrices with scaled (by two) off-diagonals      *
*     The input and output matrices may be                             *
*     identical.                                                       *
*                                                                      *
*     calling arguments:                                               *
*     nSym    : input, integer                                         *
*               number of symmetry blocks                              *
*     nBas    : input, array of integers                               *
*               matrix dimension per symmetry block                    *
*     A       : input, array of real*8                                 *
*               Unfolded input matrix (UT-storage)                     *
*     B       : output, array of real*8                                *
*               Folded output matrix  (UT-storage)                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension nBas(*) , A(*) , B(*)

      Parameter (two = 2.0d0)

*************************************
      iit(j) = j*(j+1)/2
*************************************
      ijt(i,j) = i*(i-1)/2 + j
*************************************


      iOff = 0

      Do iSym=1,nSym

        Do j=1,nBas(iSym)

          Do i=j+1,nBas(iSym)

            B( iOff + ijt(i,j) ) = two*A( iOff + ijt(i,j) )

          End Do

          B( iOff + iit(j) ) =  A( iOff + iit(j) )

        End Do

        iOff = iOff + iit(nBas(iSym))

      End Do

c      ij=0
c      do isym=1,nsym
c       do i=1,nbas(isym)
c        do j=1,i-1
c         ij=ij+1
c         B(ij)=2*A(ij)
c        end do
c        ij=ij+1
c        B(ij)=A(ij)
c       end do
c      end do


      Return
      End
