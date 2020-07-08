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
      Subroutine KE_exchange(N1,N2,lant,t,u,OPT,HEXC)
c this function computes the exchange+covalent contributions to Hamiltonian of a given Lanthanide
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)           :: N1, N2, OPT, lant
      Real(kind=8), intent(in)     :: t, u
      Complex(kind=8), intent(out) :: HEXC(N1,N1,N2,N2)
c local variables
      Integer       :: i1,j1,i2,j2,JLn,SR,ms1,ns1,ms2,ns2
      Real(kind=8) :: HEXC1(N1,N1,N2,N2)
      Real(kind=8) :: test1, test2
      Real(kind=8) :: Jfinal(0:7, -5:5, 0:1, -1:1)
      Real(kind=8) :: WCG ! Clebsh_Gordan Coefficeints
      External      :: WCG
#include "stdalloc.fh"
#include "jcoeff.fh"

      HEXC1=0.0_wp
      HEXC=(0.0_wp,0.0_wp)
      JLn=N1-1  !nexch(iLn)-1
      SR =N2-1  !nexch(iRad)-1
      Do i1=1,N1
      ms1=-(N1-1)+2*(i1-1)
      Do j1=1,N1
      ns1=-(N1-1)+2*(j1-1)
      Do i2=1,N2
      ms2=-(N2-1)+2*(i2-1)
      Do j2=1,N2
      ns2=-(N2-1)+2*(j2-1)


      If(OPT.eq.1) Then ! FULL calculation

      Do iLS=1,4
         Do iJ=1,17
            Do iK=0,7
               Do ika=-5,5
                  Do iP=0,1
                     Do iph=-1,1
                     test1=0.0_wp
                     test1=WCG(JLn, JLn, 2*iK, 0,  JLn, JLn )
     &                    *WCG( SR,  SR, 2*iP, 0,   SR,  SR )
                     test2=0.0_wp
                     test2=WCG(JLn,  ns1, 2*iK, 2*ika, JLn,  ms1 )
     &                    *WCG( SR,  ns2, 2*iP, 2*iph,  SR,  ms2 )
                     If(test1.eq.0.0_wp) Go To 104
                     If(test2.eq.0.0_wp) Go To 104

                   HEXC1(i1,j1,i2,j2) = HEXC1(i1,j1,i2,j2)
     &                           + t*t/(u+dE(lant,iLS,iJ))
     &                           * Jx(lant,iLS,iJ, iK,ika,iP,iph)
     &                           * WCG(JLn, ns1, 2*iK, 2*ika, JLn, ms1 )
     &                           * WCG( SR, ns2, 2*iP, 2*iph,  SR, ms2 )
     &                           / WCG(JLn, JLn, 2*iK,     0, JLn, JLn )
     &                           / WCG( SR,  SR, 2*iP,     0,  SR,  SR )

 104                 continue
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do


      If ( (i1.eq.1).and.(j1.eq.1).and.(i2.eq.1).and.(j2.eq.1) ) Then
      Jfinal=0.0_wp
      Write(6,'(A)') 'Jfinal(iK,ika,iP,iph):'
      Do iK=0,7
         Do ika=-5,5
            Do iP=0,1
               Do iph=-1,1

                  Do iLS=1,4
                     Do iJ=1,17

                 Jfinal(iK,ika,iP,iph) = Jfinal(iK,ika,iP,iph)
     &                           + t*t/(u+dE(lant,iLS,iJ))
     &                           * Jx(lant,iLS,iJ, iK,ika,iP,iph)

                     End Do
                  End Do
      If ( ABS(Jfinal(iK,ika,iP,iph)).gt.1.d-13 ) Then
         Write(6,'(4i4,4x,2E24.14)') iK,ika,iP,iph,
     &                              Jfinal(iK,ika,iP,iph),
     &                              Jfinal(iK,ika,iP,iph)*0.123984193_wp
      End If
               End Do
            End Do
         End Do
      End Do
      End If

      Else If(OPT.eq.2) Then ! 1/U approximation

      Do iLS=1,4
         Do iJ=1,17
            Do iK=0,7
               Do ika=-5,5
                  Do iP=0,1
                     Do iph=-1,1
                     test1=0.0_wp
                     test1=WCG(JLn, JLn, 2*iK, 0,  JLn, JLn )
     &                    *WCG( SR,  SR, 2*iP, 0,   SR,  SR )
                     test2=0.0_wp
                     test2=WCG(JLn,  ns1, 2*iK, 2*ika, JLn,  ms1 )
     &                    *WCG( SR,  ns2, 2*iP, 2*iph,  SR,  ms2 )
                     If(test1.eq.0.0_wp) Go To 105
                     If(test2.eq.0.0_wp) Go To 105

                   HEXC1(i1,j1,i2,j2) = HEXC1(i1,j1,i2,j2)
     &                           + t*t/u
     &                           * Jx(lant,iLS,iJ, iK,ika,iP,iph)
     &                           * WCG(JLn, ns1, 2*iK, 2*ika, JLn, ms1 )
     &                           * WCG( SR, ns2, 2*iP, 2*iph,  SR, ms2 )
     &                           / WCG(JLn, JLn, 2*iK,     0, JLn, JLn )
     &                           / WCG( SR,  SR, 2*iP,     0,  SR,  SR )
 105                 continue
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

      End If !( OPT )

      ! kind=8, complex double precision
      HEXC(i1,j1,i2,j2) = cmplx(HEXC1(i1,j1,i2,j2),0.0_wp, wp)

      End Do !j2
      End Do !i2
      End Do !j1
      End Do !i1

      Call mma_deallocate(Jx)
      Return
      End

