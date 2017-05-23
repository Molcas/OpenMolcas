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
C
C----------------------------------------------------------------------|
C
      subroutine AODKHEXP(NBasis,M,MX,M0,Ep,E0,EL,ES,OL,
     &                    W,TA,TB,A1,A2,O2,E2,F2,WS)
      Implicit Real*8(A-H,O-Z)
C
C     Do the Douglas-Kroll-Hess Arbitrary order term
C
C     Storage requirements:
C     Ep,E0 - NBasis
C     EL,ES,OL - NBasis*NBasis
C     W,TA,TB,A1,A2 - NBasis*NBasis
C     O2,E2,F2 - Nbasis*Nbasis*M (M=order of DKH)
C
      Dimension Ep(*),E0(*),EL(NBasis,*),ES(NBasis,*),OL(NBasis,*),
     &  W(NBasis,*),A1(NBasis,*),A2(NBasis,*),TA(NBasis,*),TB(NBasis,*),
     &  O2(NBasis,NBasis,*),E2(NBasis,NBasis,*),F2(NBasis,NBasis,*),
     &  WS(Nbasis,NBasis,*)
      Logical Ifodd
      Save Zero, One, Half
      Data Zero/0.d0/, One/1.d0/, Half/0.5d0/
C
C     Copy 1st order DKH effective potential
C     E2(upper-left) O2(upper-right) F2(lower-right)
C
      Do 20 J = 1, NBasis
        Do 20 I = 1, NBasis
          E2(I,J,1) = EL(I,J)
          F2(I,J,1) = ES(I,J)
   20     O2(I,J,1) = OL(I,J)
C
C     Apply [M/2] times Douglas-Kroll-Hess Unitary transformations
C
      Do 30 Iut = 1, M/2
        Do 35 I = 1, NBasis
          Do 35 J = 1, NBasis
            W(J,I) = O2(J,I,Iut)/(Ep(I) + Ep(J))
   35       If(Iut.le.MX) WS(J,I,Iut*2-1) = W(J,I)
C      ! Apply W_{Iut} to O/E_{Ks}
        Do 40 Ks = M-Iut, 1, -1
C         ! W_{1} only apply to O/E_{1}
          If(Iut.eq.1.and.Ks.ge.2) Goto 40
          Do 50 Ioe = 1, 2
C           ! O_{k,k<Iut} was eliminated
            If(Ioe.eq.1.and.Ks.lt.Iut) Goto 50
            Ifodd = Ioe.eq.1
            Do 55 I = 1, NBasis
              Do 55 J = 1, NBasis
C               ! Copy O/E_{Ks} terms to temp arrays
                If(Ifodd)then
                  TA(J,I) = O2(J,I,Ks)
                Else
                  TA(J,I) = E2(J,I,Ks)
                  TB(J,I) = F2(J,I,Ks)
                Endif
   55           Continue
            Do 60 K = Ks, M, Iut
C             ! skip terms do not contribute to final DKH Hamiltonian (even,upper-left)
              If(K+Iut+Iut.gt.M.and.(K+Iut.gt.M.or..not.Ifodd)) Goto 60
              KK = (K - Ks)/Iut + 1
              If(Ioe.eq.1.and.Ks.eq.Iut)then
C               ! see Eq.(74) of JCP130(2009)044102
                If(KK.eq.1)Then
                  Cof = Half
                Else
                  Cof = Dble(KK)/(KK*KK - One)
                Endif
              Else
C               ! see Eq.(71) of JCP130(2009)044102
                Cof = One/KK
              Endif
              If(Ifodd)then
C               ! skip terms do not contribute to final DKH Hamiltonian (even,upper-left)
                If(K+Iut+Iut+Iut.le.M)then
                  Call DGEMM_('T','N',NBasis,NBasis,NBasis,Cof,
     $                       W,NBasis,TA,NBasis,Zero,A2,NBasis)
                Endif
                Call DGEMM_('N','T',NBasis,NBasis,NBasis,Cof,
     $                     W,NBasis,TA,NBasis,Zero,A1,NBasis)
C               ! ( 0  W)(0  O)   (0  O)( 0  W)   ( WO'+(WO')'     0       )
C               ! (-W' 0)(O' 0) - (O' 0)(-W' 0) = (    0       -W'O-(W'O)' )
C               !  where W'=W^{\dag} O'=O^{\dag}
                Do 62 I = 1, NBasis
                  Do 62 J = 1, NBasis
                    If(K+Iut+Iut+Iut.le.M)then
                      TB(J,I) =-A2(J,I) - A2(I,J)
                      F2(J,I,K+Iut) = F2(J,I,K+Iut) + TB(J,I)
                    Endif
                    TA(J,I) = A1(J,I) + A1(I,J)
   62               E2(J,I,K+Iut) = E2(J,I,K+Iut) + TA(J,I)
              Else
                Call DGEMM_('N','N',NBasis,NBasis,NBasis,Cof,
     $                     W,NBasis,TB,NBasis,Zero,A1,NBasis)
                Call DGEMM_('N','N',NBasis,NBasis,NBasis,Cof,
     $                     TA,NBasis,W,NBasis,Zero,A2,NBasis)
C               ! ( 0  W)(E 0)   (E 0)( 0  W)   (    0     WF-EW )
C               ! (-W' 0)(0 F) - (0 F)(-W' 0) = ( (WF-EW)'   0   )
C               !  where W'=W^{\dag}
                Do 64 I = 1, NBasis
                  Do 64 J = 1, NBasis
                    TA(J,I) = A1(J,I) - A2(J,I)
   64               O2(J,I,K+Iut) = O2(J,I,K+Iut) + TA(J,I)
              Endif
              Ifodd = .not.Ifodd
   60       Continue
   50     Continue
   40   Continue
   30 Continue
C
C     Sum all even terms to Douglas-Kroll-Hess Hamiltonian
C
      Do 69 I = 1, NBasis
   69   EL(I,I) = EL(I,I) + E0(I)
      Do 70 I = 2, M0
        Do 70 J = 1, NBasis
          Do 70 K = 1, Nbasis
   70       EL(K,J) = EL(K,J) + E2(K,J,I)
C
      Return
      End
