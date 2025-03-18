!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
      subroutine CalcAXX(AXX,W)
      use Constants, only: Zero,Two,Four
      use input_mclr, only: nRoots
      Implicit None
!*****Input
      Real*8,DIMENSION((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2)::W
!*****Output
      Real*8,DIMENSION(((nRoots-1)*nRoots/2)**2)::AXX

!*****Auxiliary Quantities
      INTEGER K,L,M,N,IKL,IMN,IKL2,IMN2,IKK,ILL,IMM,INN,IC,nRTri
      Real*8  VKLMN,VLKNM,VKLNM,VLKMN

      nRTri=(nRoots-1)*nRoots/2
      DO K=1,nRoots
      DO L=1,K-1
       IKL=(K-1)*K/2+L
       IKK=(K+1)*K/2
       ILL=(L+1)*L/2
       IKL2=(K-2)*(K-1)/2+L
        Do M=1,nRoots
        Do N=1,M-1
         IMN=(M-1)*M/2+N
         IMM=(M+1)*M/2
         INN=(N+1)*N/2
         IMN2=(M-2)*(M-1)/2+N
         VKLMN=Zero
         VLKNM=Zero
         VLKMN=Zero
         VKLNM=Zero
         IF(L.eq.M) THEN
          If(N.lt.K) Then
           IC=(K-1)*K/2+N
          Else
           IC=(N-1)*N/2+K
          End If
          VKLMN=W(IC,IKK)+W(IC,INN)-Two*W(IC,ILL)-Four*W(IKL,IMN)
         END IF
         IF(K.eq.N) THEN
          If(M.lt.L) Then
           IC=(L-1)*L/2+M
          Else
           IC=(M-1)*M/2+L
          End If
          VLKNM=W(IC,ILL)+W(IC,IMM)-Two*W(IC,IKK)-Four*W(IKL,IMN)
         END IF
         IF(K.eq.M) THEN
          If(N.lt.L) Then
           IC=(L-1)*L/2+N
          Else
           IC=(N-1)*N/2+L
          End If
          VLKMN=W(IC,ILL)+W(IC,INN)-Two*W(IC,IKK)-Four*W(IKL,IMN)
         END IF
         IF(L.eq.N) THEN
          If(M.lt.K) Then
           IC=(K-1)*K/2+M
          Else
           IC=(M-1)*M/2+K
          End If
          VKLNM=W(IC,IKK)+W(IC,IMM)-Two*W(IC,ILL)-Four*W(IKL,IMN)
         END IF
         AXX((IKL2-1)*nRTri+IMN2)=VKLMN+VLKNM-VKLNM-VLKMN
        End Do
        End Do
       END DO
       END DO
      END SUBROUTINE CalcAXX
