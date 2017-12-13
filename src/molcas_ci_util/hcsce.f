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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine HCSCE(N,H,S,C,E,M)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Solve the secular equations HC=SCE. This routine is part of the  *
*     Davidson diagonalization procedure used to find the lowest roots *
*     of the CI-Hamiltonian. Because of that it is crucial to use      *
*     Schmidt orthogonalization such that the latest vector (in chro-  *
*     nological order) remains unchanged and the previous are ortho-   *
*     gonalized relative to it. It is also important that the input    *
*     data remain unchanged.                                           *
*                                                                      *
*     calling arguments:                                               *
*     N       : Type integer, input.                                   *
*               Dimensions of the secular equations.                   *
*     H       : Type double precision real, input.                     *
*               Hamiltonian in matrix representation.                  *
*     S       : Type double precision real, input.                     *
*               Overlap matrix.                                        *
*     C       : Type double precision real, output.                    *
*               Matrix containing the eigenvectors.                    *
*     E       : Type double precision real, output.                    *
*               Vector containing the eigenvalues.                     *
*     M       : Type integer, output.                                  *
*               This is the number of linearly independent basis       *
*               vectors that span space the metric given by the        *
*               overlap matrix.                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8(A-H,O-Z)

      Dimension H(N*(N+1)/2),S(N*(N+1)/2),C(N,N),E(N)
      dimension WGronk(2)

#include "WrkSpc.fh"
#include "timers.fh"

*      Character*12 method

      Call qEnter('HCSCE')
      Call Timing(Longines_1,Swatch,Swatch,Swatch)

* PAM 2009: On input, M=max possible orthonormal solutions to HC=SCE
* Save it.
      MMAX=M

*     allocate temporary work space
      Call GetMem('Temp1','Allo','Real',lw1,N*N)
      Call GetMem('Temp2','Allo','Real',lw2,N*N)
      Call GetMem('Temp3','Allo','Real',lw3,N*N)
      Call GetMem('Temp4','Allo','Real',lw4,N)

*     make local copies of H and S
      Call Square(S,Work(lw1),1,N,N)
      Call Square(H,Work(lw2),1,N,N)

*     Schmidt orthogonalization
      call dcopy_(N*N,0.0d0,0,C,1)
      call dcopy_(N,1.0d0,0,C,N+1)

*      write(6,*)' HCSCE calling Schmidt.'
*      Call Schmidt(N,Work(lw1),C,Work(lw4),M)
*      write(6,*)' HCSCE back from Schmidt. M=',M
* PAM 2009: It seems that no provision is made for the case that
* the returned M, = nr of ON vectors produced, is smaller than N?!
* Also, the whole thing looks very inefficient. But I just make
* some provisional changes now (090216).
*      write(6,*)' HCSCE check eigenvalues. N=',N
*      Call eigv(N,Work(lw1))
*      write(6,*)' HCSCE calling NewGS.'
      Call NewGS(N,Work(lw1),C,Work(lw4),M)
*      write(6,*)' HCSCE back from NewGS. M=',M
* Possibly in very difficult cases, NewGS produced too many vectors:
       M=MIN(M,MMAX)

*     transform H to an orthogonal basis
* PAM 2009: Rewritten, use only M orthogonal vectors
      Call DGEMM_('N','N',
     &            N,M,N,
     &            1.0d0,Work(lw2),N,
     &            C,N,
     &            0.0d0,Work(lw3),N)
      Call DGEMM_('T','N',
     &            M,M,N,
     &            1.0d0,C,N,
     &            Work(lw3),N,
     &            0.0d0,Work(lw2),M)

* PAM 2009: Replace by DSYEV call.
*      method='Householder'
*      method='Jacobi'

*     diagonalize and extract eigenvalues
*      If ( method.eq.'Jacobi' ) then
*        Do i=1,N
*           Do j=1,i
*              Work(j+i*(i-1)/2+lw2-1)=Work(j+(i-1)*N+lw2-1)
*           End Do
*        End Do
*        Call Jacob(Work(lw2),C,N,N)
*        Call JacOrd(Work(lw2),C,N,N)
*        Do i=1,N
*           E(i)=Work(i*(i+1)/2+lw2-1)
*        End Do
*      End If
*      If ( method.eq.'Householder' ) then
*        Call Eigen_Molcas(N,Work(lw2),E,Work(lw4))
*        Call DGEMM_('N','N',
*     &              N,N,N,
*     &              1.0d0,C,N,
*     &              Work(lw2),N,
*     &              0.0d0,Work(lw3),N)
*        call dcopy_(N*N,Work(lw3),1,C,1)
*      End If
* PAM 2009 Equivalent, DSYEV, note now use just M, not all N:
      INFO=0
      call dsyev_('V','L',M,Work(lw2),M,E,WGRONK,-1,INFO )
      NSCRATCH=INT(WGRONK(1))
      CALL GETMEM('SCRATCH','ALLO','REAL',LSCRATCH,NSCRATCH)
      call dsyev_('V','L',M,WORK(lw2),M,E,WORK(LSCRATCH),
     &            NSCRATCH,INFO)
      CALL GETMEM('SCRATCH','FREE','REAL',LSCRATCH,NSCRATCH)
      Call DGEMM_('N','N',
     &              N,M,M,
     &              1.0d0,C,N,
     &              Work(lw2),M,
     &              0.0d0,Work(lw3),N)
      call dcopy_(N*M,Work(lw3),1,C,1)

*     deallocate temporary work space
      Call GetMem('Temp4','Free','Real',lw4,N)
      Call GetMem('Temp3','Free','Real',lw3,N*N)
      Call GetMem('Temp2','Free','Real',lw2,N*N)
      Call GetMem('Temp1','Free','Real',lw1,N*N)

      Call Timing(Longines_2,Swatch,Swatch,Swatch)
      Longines_2 = Longines_2 - Longines_1
      Longines_3 = Longines_3 + Longines_2
      Call qExit('HCSCE')

      Return
      End

      Subroutine NewGS(N,S,C,Temp,M)
      Implicit Real*8(A-H,O-Z)
      Dimension S(N,N),C(N,N),Temp(N)
* PAM 2009, New Gram-Schmidt 090216
#include "warnings.fh"

      M=0
      do i=1,N
       X=S(i,i)
       If(X.lt.1.0D-6) Goto 90
       Y=1.0D0/sqrt(X)
       call dcopy_(N,0.0D0,0,C(1,M+1),1)
       C(i,M+1)=Y
       call dcopy_(N,S(1,i),1,Temp,1)
       Call DSCAL_(N,Y,Temp,1)

       Loop=0
  10   Continue
       Loop=Loop+1
       Do k=1,M
        ovl=DDOT_(N,Temp,1,C(1,k),1)
        call daxpy_(N,-ovl,C(1,k),1,C(1,M+1),1)
       End Do
       Call dGeMV_('N',N,N,1.0D0,S,N,C(1,M+1),1,0.0D0,Temp,1)
       xn2=DDOT_(N,Temp,1,C(1,M+1),1)

       If(xn2.lt.1.0D-6) Goto 90

       Y=1.0D0/sqrt(xn2)
       Call DSCAL_(N,Y,C(1,M+1),1)
       Call dGeMV_('N',N,N,1.0D0,S,N,C(1,M+1),1,0.0D0,Temp,1)
       If(Loop.eq.1 .and. Y.gt.100.0D0) Goto 10
*MGD issues with many states : to be very safe, test the result
       isfail=0
       Do  k=1,M
        ovl=DDOT_(N,Temp,1,C(1,k),1)
        if (abs(ovl).gt.1.0d-4) isfail=1
        if (abs(ovl).gt.1.0d-4) write(6,*) 'failed for', M,k,ovl
       End Do
       If (isfail.eq.0) M=M+1

  90   Continue
      End do
      Return
      End
      Subroutine Schmidt(N,S,C,Temp,M)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Schmidt orthogonalization such that the latest vector (in chro-  *
*     nological order) remains unchanged and the previous are ortho-   *
*     gonalized relative to it.                                        *
*                                                                      *
*     calling arguments:                                               *
*     N       : Type integer, input.                                   *
*               Dimensions of the overlap matrix.                      *
*     S       : Type double precision real, input.                     *
*               Overlap matrix.                                        *
*     C       : Type double precision real, output.                    *
*               Matrix containing the transformation vectors.          *
*     Temp    : Type double precision real, input/output.              *
*               Scratch area.                                          *
*     M       : Type integer, output.                                  *
*               This is the number of linearly independent basis       *
*               vectors that span space the metric given by the        *
*               overlap matrix.                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8(A-H,O-Z)

      Dimension S(N,N),C(N,N),Temp(N)

      Logical forward,backward

*     Call qEnter('Schmidt')

      forward  = .true.
      backward = .not.forward

      M=0
      Do i=1,N
         Do j=1,N
            C(j,i)=0.0d0
         End Do
         C(i,i)=1.0d0/Sqrt(S(i,i))
      End Do

*     forward orthonormalization
*     (vectors 1 remains unchanged)
      If ( forward ) then
        Do i=1,N
           Alpha=C(i,i)
           Do j=1,N
              Temp(j)=S(j,i)*Alpha
           End Do
           Do j=1,i-1
              Sum=0.0d0
               Do k=1,i
                 Sum=Sum+C(k,j)*Temp(k)
              End Do
              Do k=1,i
                 C(k,i)=C(k,i)-Sum*C(k,j)
              End Do
           End Do
           Sum=0.0d0
           Do k=1,i
              Sum=Sum+C(k,i)*Temp(k)
           End Do
           If ( Sum.gt.1.0d-9 ) then
              M=M+1
              Alpha=1.0d0/Sqrt(Sum)
              Do k=1,i
                 C(k,i)=C(k,i)*Alpha
              End Do
           Else
              Do k=1,i
                 C(k,i)=0.0d0
              End Do
           End If
        End Do
      End If

*     backward orthonormalization
*     (vectors N remains unchanged)
      If ( backward ) then
        Do i=N,1,-1
           Alpha=C(i,i)
           Do j=1,N
              Temp(j)=S(j,i)*Alpha
           End Do
           Do j=N,i+1,-1
              Sum=0.0d0
              Do k=j,N
                 Sum=Sum+C(k,j)*Temp(k)
              End Do
              Do k=j,N
                 C(k,i)=C(k,i)-Sum*C(k,j)
              End Do
           End Do
           Sum=0.0d0
           Do k=i,N
              Sum=Sum+C(k,i)*Temp(k)
           End Do
           If ( Sum.gt.1.0d-9 ) then
              M=M+1
              Alpha=1.0d0/Sqrt(Sum)
              Do k=i,N
                 C(k,i)=C(k,i)*Alpha
              End Do
           Else
              Do k=i,N
                 C(k,i)=0.0d0
              End Do
           End If
        End Do
      End If

*     Call qExit('Schmidt')

      Return
      End
*      subroutine eigv(N,S)
*      implicit real*8 (a-h,o-z)
*      dimension S(N,N)
*      dimension WGronk(2)

*
*      CALL GETMEM('LSCP','ALLO','REAL',LSCP,N**2)
*      CALL DCOPY_(N**2,S,1,WORK(LSCP),1)
*      CALL GETMEM('LEIG','ALLO','REAL',LEIG,N)
*      INFO=0
*      call dsyev_('V','L',N,WORK(LSCP),N,WORK(LEIG),WGRONK,-1,INFO )
*      NSCRATCH=INT(WGRONK(1))
*      CALL GETMEM('SCRATCH','ALLO','REAL',LSCRATCH,NSCRATCH)
*      call dsyev_('V','L',N,WORK(LSCP),N,WORK(LEIG),WORK(LSCRATCH),
*     &            NSCRATCH,INFO)
*      CALL GETMEM('SCRATCH','FREE','REAL',LSCRATCH,NSCRATCH)
*      write(6,*)' Largest eigenvalues of S matrix:'
*      write(6,'(5g16.8)')(Work(leig+i),i=MAX(0,N-10),N-1)
*      CALL GETMEM('LEIG','FREE','REAL',LEIG,N)
*      CALL GETMEM('LSCP','FREE','REAL',LSCP,N**2)
*
*      return
*      end
