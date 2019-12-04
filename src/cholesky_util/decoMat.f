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
* Copyright (C) 2016, Giovanni Li Manni                                *
************************************************************************
      Subroutine DecoMat(MAT,dimens,eigenvec,NumV,rc)
************* by G. Li Manni Stuttgart March 2016 *************
*
* MAT:      copy of the D1A matrix (1RDM) in AO basis.
*           It will be destroyed in eigen_molcas
*
* DIMENS:   dimension of the square MAT = nbas(isym)
*
* EIGENVEC: First it is used as scratch inside eigen_molcas. Next it will contains eigenvectors
*           of the eigen-decomposed matrix in Cholesky form. Instead of eigenvectors X such that MAT = XDX^t,
*           eigenvectors are presented as Y = X(D**0.5) such that MAT = YY^t.
*           Negative eigenvalues are set to zero and eigenvalue larger than 2.0d0 set to 2.0d0
*
* NUMV    : Number of non negative eigenvalues
      implicit none
      integer dimens,NumV,NumVnull,rc,i,j
      real*8 MAT(dimens,dimens),eigenvec(dimens,dimens)
      real*8 eigenval(dimens)

      Character*12 routine
      Parameter (routine = 'DecoNegatMat')

      Call qEnter(routine)
      rc = 0
      NumV = 0
      NumVnull = 0

      If (dimens .lt. 1) then
       rc= -1
       write(6,*) 'matrix size < 1'
       Go To 10
      end if

* Step 1: diagonalize MAT. MAT is destroyed and replaced by the eigenvectors values.
* IMPORTANT: At this stage eigenvec is used as scratch.
      call eigen_molcas(dimens,MAT,eigenval,eigenvec)
* Move MAT to eigenvec. where it belongs. I could not wait!
      call dcopy_(dimens**2,MAT,1,eigenvec,1)
* Set to zero negative eigenvalue and to TWO values larger than 2.0d0.
* Count only the positive ones (counter NumV)
      do j = 1, dimens
       if(eigenval(j).gt.1.0d-12) then
          NumV = NumV + 1
          if(eigenval(j).gt.2.0d0) eigenval(j) = 2.0d0
       else
        eigenval(j) = 0.0d0
       end if
      end do
* Sort eigenvalues in decreasing order of occupation number
      call IncrSort(eigenval,eigenvec,dimens)
* Compute sqrt(eigenvalues)
      do i = 1, dimens
          eigenval(i) = sqrt(eigenval(i))
      end do
* Generate Y = X(D**0.5)
      do j = 1, dimens
          do i = 1, dimens
           eigenvec(i,j) = eigenvec(i,j)*eigenval(j)
          end do
      end do
****************** Exit ****************
10    Continue
      Call qExit(routine)
      return
      end subroutine

      SUBROUTINE IncrSort(EVal,EVec,dimens)
      Implicit none
      integer dimens, i, k, j, l
      Real*8 EVal(dimens),EVec(dimens,dimens),swap
      Call qEnter('IncrSort')
      Do i = 1,dimens - 1
         k = i
         do j = i + 1, dimens
            If (EVal(j).gt.EVal(k)) k = j
         end do
         If (k.ne.i) Then
            Swap    = EVal(k)
            EVal(k) = EVal(i)
            EVal(i) = Swap
            do l = 1, dimens
               Swap      =   EVec(l,k)
               EVec(L,K) =   EVec(l,i)
               EVec(L,I) =   Swap
            end do
         End If
      End Do
      Call qExit('IncrSort')
      Return
      End Subroutine
