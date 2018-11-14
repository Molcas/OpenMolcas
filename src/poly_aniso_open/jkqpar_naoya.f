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
      Subroutine JKQPar_Naoya(N1,N2,HEXCH,Jpar)
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer, intent(in)           :: N1, N2
      Complex(kind=wp), intent(in)  :: HEXCH(N1,N1,N2,N2)
      Complex(kind=wp), intent(out) :: Jpar( (N1-1), (-N1+1):(N1-1),
     &                                       (N2-1), (-N2+1):(N2-1) )
      ! local variables
      Integer          :: is1,is2,js1,js2,ms1,ms2,ns1,ns2,
     &                    k1,k2,q1,q2,j1,j2
      Complex(kind=wp) :: QMAT(N1,N1,N2,N2), trace
      Real(kind=wp)    :: WCG, OPER, FACT
      External         :: WCG
      Logical DBG
      DBG=.false.

c we need to project now the HEXCH: in products of ITOs
c  HEXCH = SUM(rank1,proj1,rank2,proj2)=
c         { B(rank1,proj1,rank2,proj2)* O1(rank1,proj1) * O2(rank2,proj2) }
c Naoya definition
c eq.40 in DoI:10.1103/PhysRevB.91.174438
      Jpar=(0.0_wp,0.0_wp)
      J1=N1-1 ! i.e. Double of the dimension of the spin on site 1
      J2=N2-1 ! i.e. Double of the dimension of the spin on site 2
      Do k1=0,J1
        Do q1=-k1,k1
          Do k2=0,J2
            Do q2=-k2,k2
              If( mod(k1,2).ne.1) Go To 104 ! the rank of individual spins must be even
              If( mod(k2,2).ne.1) Go To 104 ! the rank of individual spins must be even
              ! If the total rank is odd, Then it is a local ZFS contribution; ==> to be Done later
              QMAT=0.0_wp
              ! compute the qmat:
              ! projections q and q' are with opposite sign:
              !  -is1 and  -is2
              Do is1=1,N1
                ms1=2*is1-N1-1 ! spin projection on site 1
                Do js1=1,N1
                  ns1=2*js1-N1-1 ! spin projection on site 1
                  Do is2=1,N2
                    ms2=2*is2-N2-1 ! spin projection on site 2
                    Do js2=1,N2
                      ns2=2*js2-N2-1 ! spin projection on site 2
                      If ( WCG(J1,J1,2*k1,0,J1,J1)==0) Go To 103
                      If ( WCG(J2,J2,2*k2,0,J2,J2)==0) Go To 103

                      QMAT(is1,js1,is2,js2)=
     &                               WCG( J1, ns1, 2*k1,-2*q1, J1, ms1 )
     &                             * WCG( J2, ns2, 2*k2,-2*q2, J2, ms2 )
     &                             / WCG( J1,  J1, 2*k1,    0, J1,  J1 )
     &                             / WCG( J2,  J2, 2*k2,    0, J2,  J2 )

       If ( DBG .and. (ABS(QMAT(is1,js1,is2,js2)).gt.1.d-12) ) Then
          Write(6,'(8(A,i3),A,2F20.14)')
     &                       ' QMAT(',k1,',', q1,',', k2,',', q2,'|||',
     &                               is1,',',js1,',',is2,',',js2,') = ',
     &                         QMAT(is1,js1,is2,js2)
       End If

 103                Continue
                    End Do
                  End Do
                End Do
              End Do
              ! compute the trace Tr[QMAT * HEXCH]
              trace=(0.0_wp,0.0_wp)
              Do is1=1,N1
                Do js1=1,N1
                  Do is2=1,N2
                    Do js2=1,N2
                      trace = trace +  QMAT(is1,js1,is2,js2)
     &                              * HEXCH(js1,is1,js2,is2)
                    End Do
                  End Do
                End Do
              End Do

              If ( DBG .AND. (ABS(trace).gt.1.d-12) ) Then
                  Write(6,'(4(A,i3),A,2F20.14)')
     &                    'trace(',k1,',',q1,',',k2,',',q2,') = ',trace
              End If

              OPER=0.0_wp
              FACT=0.0_wp
              OPER= WCG(J1,J1,2*k1,0,J1,J1)*WCG(J1,J1,2*k1,0,J1,J1)
     &             *WCG(J2,J2,2*k2,0,J2,J2)*WCG(J2,J2,2*k2,0,J2,J2)
              FACT=dble(((-1)**(q1+q2))*((2*k1+1)*(2*k2+1)))/dble(N1*N2)

              Jpar(k1,q1,k2,q2)=cmplx(FACT*OPER,0.0_wp,wp)*trace

              If( DBG.and. (ABS(Jpar(k1,q1,k2,q2)).gt.0.5d-13) ) Then
                  Write(6,'(4(A,i3),A,2F20.14)')
     &                    '    J(',k1,',',q1,',',k2,',',q2,') = ',
     &                         Jpar(k1,q1,k2,q2)
       End If

 104        Continue
            End Do
          End Do
        End Do
      End Do

      Return
      End Subroutine JKQPar_Naoya
