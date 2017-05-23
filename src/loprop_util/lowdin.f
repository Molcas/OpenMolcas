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
      Subroutine Lowdin(S,C,nDim)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
      ReaL*8 C(nDim,nDim),S(nDim,nDim)
*
      Call Allocate_Work(ip_S,nDim*(nDim+1)/2)
      Call Allocate_Work(ip_B,nDim**2)
*
      Do i = 1, nDim
         Do j = 1, i
            ioff = ip_S -1 + i*(i-1)/2 + j
            Work(iOff)=S(i,j)
         End Do
      End Do
      call dcopy_(nDim**2,Zero,0,Work(ip_B),1)
      call dcopy_(nDim,One,0,Work(ip_B),nDim+1)
*
      Call Lowdin_(S,Work(ip_S),C,nDim,nDim,Work(ip_B))
*
      Call Free_Work(ip_B)
      Call Free_Work(ip_S)
      Return
      End

      Subroutine Lowdin_(S,Eval,C,nDim,nDim2,Blk)
*
*     S: full-storage overlap matrix (it will be destroyed!)
*     C: on exit, the S^-1/2 matrix
*
      Implicit ReaL*8 (A-H,O-Z)
#include "real.fh"
      Real*8 C(nDim,nDim),S(nDim,nDim), Blk(nDim,nDim),
     &       Eval(nDim*(nDim+1)/2)
*
C     NDim2=NDim2
C
      Data DIAGTH,DANGER/1.0D-12,1.0D3/
C
C          diagth  threshold for matrix diagonalization used in
C                     subroutine jacobi.  in jacobi, this  constant
C                     is called doneth.
C          danger  criterion for deciding that the job should be
C                     aborted due to numerical problems caused by near
C                     linear dependencies in the basis set.
C                     all eigenvalues of the weighted overlap
C                     matrix must be greater than DIAGTH*DANGER.
C
      toosml=diagth*danger

C  diagonalize overlap matrix
C
      Call Jacob(Eval,Blk,nDim,nDim)
C
c     Call RecPrt('Blk',' ',Blk,nDim,nDim)
c     Call TriPrt('Eval',' ',Eval,nDim)

C  form the inverse sqrt of the overlap matrix of the vectors:
C  (avoid numerical problems of linear dependence (too small eigenvals)
C  by prescreening
C  the  eigenvalues)
C
      Do i=1,nDim
         eigenv=eval(i*(i+1)/2)
         if(eigenv.lt.toosml) go to 900
         eval(i*(i+1)/2)=ONE/sqrt(eigenv)
      End Do
*
* --- Compute S^-1/2
*
      do i=1,nDim
        do  j=1,i
          sij=zero
          do  k=1,nDim
             sij=sij+eval(k*(k+1)/2)*blk(i,k)*blk(j,k)
          enddo
          C(i,j)=sij
          C(j,i)=sij
        enddo
      enddo
C
c     Call RecPrt('C',' ',C,nDim,nDim)
C
      Return
  900 write(6,910) eigenv,toosmL
  910 format(/1X,'An eigenvalue of the overlap matrix of the ',
     *   'symmetrized Jacobi transf. ',
     *   'matrix of ',E13.5,' has been found.'/1X,
     *   'This is lower than the allowed threshold of ',E13.5)
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(S)
         Call Unused_integer(nDim2)
      End If
      End
