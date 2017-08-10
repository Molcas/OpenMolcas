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
      Subroutine Find_Dipole_Center(q_A,q_B,Dipole_A,Dipole_B,qn_A,qn_B,
     &                              R_A,R_B,EC_A,EC_B,EC_AB,t,iPlot)
      Implicit Real*8 (A-H,O-Z)
      Dimension EC_A(3),EC_B(3),EC_AB(3)
      Real*8 Multipole_Expansion
      External Multipole_Expansion,x_Golden2
      Logical Absolute
      Parameter (Absolute = .False.)
      Parameter (E_Threshold = 1.0D-12, R_Threshold = 1.0D-12)
#include "real.fh"
*
      R = R_B - R_A
*
      nSteps = 100
      Delta  = R/(DBLE(nSteps)+1.0D0)
      R_Test = R_A + Delta
      If (iPlot .eq. 1) Then
       write(6, *)'Electronic contributions: q_A, q_B = ',q_a,q_b
         DO i = 1,nSteps
            R_Test = R_A + Delta*DBLE(i)
            E = Multipole_Expansion(q_A,q_B,Dipole_A,Dipole_B,
     &                              R_A,R_B,R_Test,Absolute)
            Write(6,'(1X,A,F6.3,1X,F20.12)') 'R, E = ',R_test,E
            Call xFlush(6)
         End Do
         write(6, *)'Nuclear contributions: q_A, q_B = ',qn_A,qn_B
         DO i = 1,nSteps
            R_Test = R_A + Delta*DBLE(i)
            E = Multipole_Expansion(qn_A,qn_B,Zero,Zero,R_A,R_B,
     &                              R_Test,Absolute)
            Write(6,'(1X,A,F6.3,1X,F20.12)') 'R, E = ',R_test,E
            Call xFlush(6)
         End Do
       write(6, *)'Total contributions: q_A, q_B = ',q_a+qn_a,q_b+qn_b
         DO i = 1,nSteps
            R_Test = R_A + Delta*DBLE(i)
            E = Multipole_Expansion(q_A+qn_A,q_B+qn_B,Dipole_A,Dipole_B,
     &                              R_A,R_B,R_Test,Absolute)
            Write(6,'(1X,A,F6.3,1X,F20.12)') 'R, E = ',R_test,E
            Call xFlush(6)
         End Do
      End If
*
      ax = (R_A+R_B)*.5 + Delta
      bx = (R_A+R_B)*.5 - Delta
* First make an initial bracketing of the minima
      Call n_mnBrak2(ax,bx,cx,fa,fb,fc,multipole_expansion,q_a,q_b,
     &             dipole_a,dipole_b,r_a,r_b)
* Use Golden to find the minima
      golden = x_golden2(ax,bx,cx,multipole_expansion,r_threshold,
     &             e_threshold,r_best,q_a,q_b,
     &             dipole_a,dipole_b,r_a,r_b)
      t = (r_best-half*R)/R
      write(6,'(A,3F18.10)')
     &       't_el , r_best, golden = ',t, r_best, golden
      call xflush(6)
*
*
*
      ax = (R_A+R_B)*.5 + Delta
      bx = (R_A+R_B)*.5 - Delta
C      ax = -0.245
C      bx = -0.263
* First make an initial bracketing of the minima
      Call n_mnBrak2(ax,bx,cx,fa,fb,fc,multipole_expansion,qn_a,qn_b,
     &             zero,zero,r_a,r_b)
* Use Golden to find the minima
      golden = x_golden2(ax,bx,cx,multipole_expansion,r_threshold,
     &             e_threshold,r_best2,qn_a,qn_b,
     &             zero,zero,r_a,r_b)
      t = (r_best2-half*R)/R
      write(6,'(A,3F18.10)')
     &       't_nuc, r_best, golden = ',t, r_best2, golden
      call xflush(6)
      r_best = (abs(q_a+q_b)*R_best+abs(qn_a+qn_b)*R_best2)
     &               /(abs(q_a+q_b)+abs(qn_a+qn_b))
      t = (r_best-half*R)/R
      write(6,'(A,3F18.10)')
     &       't_fit, r_best, golden = ',t, r_best, golden
      call xflush(6)
*
*
*
C??      ax = (R_A+R_B)*.5 + Delta
C??      bx = (R_A+R_B)*.5 - Delta
C??C      ax = -0.245
C??C      bx = -0.263
C??* First make an initial bracketing of the minima
C??      Call n_mnBrak2(ax,bx,cx,fa,fb,fc,multipole_expansion,q_a+qn_a,
C??     &             q_b+qn_b,dipole_a,dipole_b,r_a,r_b)
C??* Use Golden to find the minima
C??      golden = x_golden2(ax,bx,cx,multipole_expansion,r_threshold,
C??     &             e_threshold,r_best,q_a+qn_a,q_b+qn_b,
C??     &             dipole_a,dipole_b,r_a,r_b)
C??      t = (r_best-half*R)/R
C??      write(6,'(A,3F18.10)')
C??     &       't_sum, r_best, golden = ',t, r_best, golden
C??      call xflush(6)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(EC_A)
         Call Unused_real_array(EC_B)
         Call Unused_real_array(EC_AB)
      End If
      End
