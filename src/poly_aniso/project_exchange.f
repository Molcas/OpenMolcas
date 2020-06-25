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
      Subroutine project_exch(N1,N2,S1,S2,M1,M2,E1,E2,HEXCH,Jpar,Jc)
c  this function determines the local pseuDospins and rotates the hamiltonian
c  to the local pseuDospin basis
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer N1,N2
!     spin-orbit energies on each site
      Real(kind=8) ::  E1(N1), E2(N2)
!     spin matrices on each site
      Complex(kind=8) ::  S1(3,N1,N1),S2(3,N2,N2)
!     magnetic moment matrices on each site
      Complex(kind=8) ::  M1(3,N1,N1),M2(3,N2,N2)
      Complex(kind=8) ::  HEXCH(N1,N1,N2,N2) ! exchange hamiltonian
      Complex(kind=8) ::  HEXCH2(N1,N1,N2,N2) ! exchange hamiltonian
      Complex(kind=8) ::  HEXCH3(N1,N1,N2,N2) ! exchange hamiltonian
c--  local variables --
      Integer ns1,ns2,is1,is2,iprint !,i1,j1,i2,j2
c      Integer k1,k2,q1,q2,js1,js2,ms1,ms2
      Real(kind=8) ::  Ethr,gtens(2,3), maxes(2,3,3),Jc(3,3)
      Complex(kind=8) ::  Z1(N1,N1), Z2(N2,N2)
      Complex(kind=8) ::  SR1(3,N1,N1),MR1(3,N1,N1)
      Complex(kind=8) ::  SR2(3,N2,N2),MR2(3,N2,N2)
      Complex(kind=8) ::  TMP(N1,N1)
c      Complex(kind=8) ::  DIP_O1(N1,N1)
c      Complex(kind=8) ::  DIP_W1(N1,N1)
c      Complex(kind=8) ::  DIP_O2(N2,N2)
c      Complex(kind=8) ::  DIP_W2(N2,N2)
c      Complex(kind=8) ::  SP_MOW1,SP_MOW2
c      Complex(kind=8) ::  QMAT(N1,N1,N2,N2) !,trace
c      Real(kind=8) ::  WCG ! Clebsh_Gordan Coefficeints
c      logical DBG
c      external WCG
      Complex(kind=8) ::  Jpar(N1-1,-N1+1:N1-1,N2-1,-N2+1:N2-1)
      Call qEnter('PA_projexch')
c      DBG=.false.

c      Write(6,'(A)') 'J parameters in the initial ab intio basis:'
c      Jpar=(0.0_wp,0.0_wp)
c      Call JKQPar(N1,N2,HEXCH,Jpar)
c      Jc=0.0_wp
c      Call tensor2cart(1,1,Jpar(1,-1:1,1,-1:1),Jc)


c determine the pseuDospin on each site (Z1 and Z2):
!     threshold for determination of the local pseuDospin main anisotropy axis
      Ethr=0.2_wp
      ns1=0
      ns2=0
      Do is1=1,N1
        If(E1(is1).lt.Ethr) Then
          ns1=ns1+1
        End If
      End Do
      Do is1=1,N2
        If(E2(is1).lt.Ethr) Then
          ns2=ns2+1
        End If
      End Do
      Write(6,'(A,i3)') 'size of local pseudospin, site 1  =',ns1
      Write(6,'(A,i3)') 'size of local pseudospin, site 2  =',ns2
      Call atens(  M1(1:3,1:ns1,1:ns1), ns1, gtens(1,1:3),
     &             maxes(1,1:3,1:3), 2 )
      Call atens(  M2(1:3,1:ns2,1:ns2), ns2, gtens(2,1:3),
     &             maxes(2,1:3,1:3), 2 )
c rotate the magnetic moment to the coordinate system of main magnetic axes on Ln
      SR1=(0.0_wp,0.0_wp)
      MR1=(0.0_wp,0.0_wp)
      SR2=(0.0_wp,0.0_wp)
      MR2=(0.0_wp,0.0_wp)
      Call rotmom2(   S1, N1, maxes(1,:,:), SR1 )
      Call rotmom2(   M1, N1, maxes(1,:,:), MR1 )
      Call rotmom2(   S2, N2, maxes(2,:,:), SR2 )
      Call rotmom2(   M2, N2, maxes(2,:,:), MR2 )
      Z1=(0.0_wp,0.0_wp)
      Z2=(0.0_wp,0.0_wp)
      iprint=1
      Call pseudospin(MR1,N1,Z1,3,1,iprint)
      Call pseudospin(MR2,N2,Z2,3,1,iprint)
c reWrite the exchange matrix in the basis of local pseuDospins:
      HEXCH2=(0.0_wp,0.0_wp)
      HEXCH3=(0.0_wp,0.0_wp)
      Do is1=1,N2
        Do is2=1,N2
          TMP(:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('C','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &                Z1(1:N1,1:N1), N1,
     &                HEXCH(1:N1,1:N1,is1,is2), N1, (0.0_wp,0.0_wp),
     &                TMP(1:N1,1:N1), N1 )
          Call ZGEMM_('N','N',N1,  N1, N1, (1.0_wp,0.0_wp),
     &                TMP(1:N1,1:N1), N1,
     &                 Z1(1:N1,1:N1), N1, (0.0_wp,0.0_wp),
     &                HEXCH2(1:N1,1:N1,is1,is2), N1 )
        End Do
      End Do
      Do is1=1,N1
        Do is2=1,N1
          TMP(:,:)=(0.0_wp,0.0_wp)
          Call ZGEMM_('C','N',N2,  N2, N2, (1.0_wp,0.0_wp),
     &                Z2(1:N2,1:N2), N2,
     &                HEXCH2(is1,is2,1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &                TMP(1:N2,1:N2), N2 )
          Call ZGEMM_('N','N',N2,  N2, N2, (1.0_wp,0.0_wp),
     &                TMP(1:N2,1:N2), N2,
     &                 Z2(1:N2,1:N2), N2, (0.0_wp,0.0_wp),
     &                HEXCH3(is1,is2,1:N2,1:N2), N2 )
        End Do
      End Do
      Jpar=(0.0_wp,0.0_wp)
      Call JKQPar(N1,N2,HEXCH3,Jpar)
      Jc=0.0_wp
      Call tensor2cart(Jpar(1,-1:1,1,-1:1),Jc)

      Call qExit('PA_projexch')
      Return
      End
