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
      Subroutine KE_Covalent(N,lant,t,u,OPT, HCOV )
c this function computes the covalent CF Hamiltonian ofr a given Lanthanide
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer N,OPT,lant
      Real(kind=8) ::  t,u
      Real(kind=8) ::  WCG ! Clebsh_Gordan Coefficeints
      Complex(kind=8) ::  HCOV(N,N)
c local variables
      Integer i,j,JLn,ms1,ns1
      Real(kind=8) ::  HCOV1(N,N)
      Real(kind=8) ::  test1
      external WCG
#include "stdalloc.fh"
#include "jcoeff.fh"

      HCOV1=0.0_wp
      HCOV=(0.0_wp,0.0_wp)
      JLn=N-1

      Do i=1,N
      ms1=-(N-1)+2*(i-1)
         Do j=1,N
         ns1=-(N-1)+2*(j-1)

      If(OPT.eq.1) Then !FULL calculation

            Do iLS=1,4
               Do iJ=1,17
                  Do iK=0,6,2
                     Do ika=-4,4,4
      test1=0.0_wp
      test1=WCG(JLn, JLn, 2*iK, 0,  JLn, JLn )
      If(test1.eq.0.0_wp) Go To 107
      HCOV1(i,j)=HCOV1(i,j)
     &                    + t*t/(u+dE(lant,iLS,iJ))
     &                    * Jx(lant,iLS,iJ, iK,ika,0,0)
     &                    * WCG(JLn, ns1, 2*iK, 2*ika, JLn, ms1 )
     &                    / WCG(JLn, JLn, 2*iK,     0, JLn, JLn )
 107  continue
                     End Do
                  End Do
               End Do
            End Do

      Else If(OPT.eq.2) Then ! 1/U approximation

            Do iLS=1,4
               Do iJ=1,17
                  Do iK=0,6,2
                     Do ika=-4,4,4
      test1=0.0_wp
      test1=WCG(JLn, JLn, 2*iK, 0,  JLn, JLn )
      If(test1.eq.0.0_wp) Go To 108
      HCOV1(i,j)=HCOV1(i,j) + t*t/u
     &                    * Jx(lant,iLS,iJ, iK,ika,0,0)
     &                    * WCG(JLn, ns1, 2*iK, 2*ika, JLn, ms1 )
     &                    / WCG(JLn, JLn, 2*iK,     0, JLn, JLn )
 108  continue
                     End Do
                  End Do
               End Do
            End Do
      End If
      ! kind=8, complex double precision
      HCOV(i,j)=cmplx(HCOV1(i,j),0.0_wp, wp)
      End Do !j
      End Do !i

      Call mma_deallocate(Jx)
      Return
      End
