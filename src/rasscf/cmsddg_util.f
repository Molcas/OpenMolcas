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
* Copyright (C) 2022, Jie J. Bao                                       *
************************************************************************
******************************************************************
* history:                                                       *
* Jie J. Bao, on Apr. 12, 2022, created this file.               *
******************************************************************

*This file contains subroutines that calculates the DDg matrix, Q_aa,
* the Gradient and Hessian of Q_aa wrt intermediate-state rotations.

************************************************************************
      Subroutine CalcDDg(DDg,GD,Dg,nDDg,nGD,lRoots2,NAC2)
******************************************************************
*     Gtuvx : two-electron integral, g_tuvx                      *
*     GD    : "generalized 1-e density matrix"                   *
*              GD^KL: transition density matrix from L to K      *
*              GD^KK: density matrix for state K                 *
*     Dg    :  sum_{vx}{GD^KL_vx * g_tuvx}                       *
*     In GDorbit and Dgorbit, the leading index is orbital index;*
*     In GDstate and Dgstate, the leading index is state index.  *
*                                                                *
*     DDg   :  sum_{tuvx}{GD^KL_tu * GD^MN_vx * g_tuvx}          *
*      namely, sum_{tu}{GD^KL_tu * Dg^MN_tu}                     *
******************************************************************
      INTEGER nDDg, nGD, lRoots2, NAC2
      Real*8 DDg(nDDg),GD(nGD),Dg(nGD)


      CALL DGEMM_('T','N',lRoots2,lRoots2,NAC2,
     &               1.0d0, Dg, NAC2, GD, NAC2,
     &               0.0d0, DDg, lRoots2)

      RETURN
      End Subroutine

************************************************************************


      Subroutine CalcHessCMS(Hess,DDg,nDDg,lRoots,nSPair)
      INTEGER nDDg,lRoots,nSPair
      Real*8 DDg(nDDg),Hess(nSPair**2)
      Real*8 Vklmn,Vlknm,Vklnm,Vlkmn
      INTEGER K,L,M,N,iKL,iMN,iLoc1,iLoc2,iLoc3,iLoc4,iLoc5,
     &        lRoots2,lRoots3,lRoots23

      lRoots2=lRoots**2
      lRoots3=lRoots2*lRoots
      lRoots23=lRoots2+lRoots3

      DO K=2,lRoots
      DO L=1,K-1
       iKL=(K-2)*(K-1)/2+L
       Do M=2,lRoots
       Do N=1,M-1
        iMN=(M-2)*(M-1)/2+N
        Vklmn=0.0d0
        Vlknm=0.0d0
        Vlkmn=0.0d0
        Vklnm=0.0d0
        iLoc5=K+(L-1)*lRoots+(M-1)*lRoots2+(N-1)*lRoots3
        IF(L.eq.M) THEN
         iLoc4=K+(N-1)*lRoots
         iLoc1=iLoc4+(K-1)*lRoots23
         iLoc2=iLoc4+(N-1)*lRoots23
         iLoc3=iLoc4+(L-1)*lRoots23
         Vklmn=DDg(iLoc1)+DDg(iLoc2)-2.0d0*DDg(iLoc3)-4.0d0*DDg(iLoc5)
C         Vklmn=DDg(K,N,K,K)+DDg(K,N,N,N)-2.0d0*DDg(K,N,L,L)
C     &         -4.0d0*DDg(K,L,M,N)
        END IF
        IF(K.eq.N) THEN
         iLoc4=L+(M-1)*lRoots
         iLoc1=iLoc4+(L-1)*lRoots23
         iLoc2=iLoc4+(M-1)*lRoots23
         iLoc3=iLoc4+(K-1)*lRoots23
         Vlknm=DDg(iLoc1)+DDg(iLoc2)-2.0d0*DDg(iLoc3)-4.0d0*DDg(iLoc5)
C         Vlknm=DDg(L,M,L,L)+DDg(L,M,M,M)-2.0d0*DDg(L,M,K,K)
C     &         -4.0d0*DDg(K,L,M,N)
        END IF
        IF(K.eq.M) THEN
         iLoc4=L+(N-1)*lRoots
         iLoc1=iLoc4+(L-1)*lRoots23
         iLoc2=iLoc4+(N-1)*lRoots23
         iLoc3=iLoc4+(K-1)*lRoots23
         Vlkmn=DDg(iLoc1)+DDg(iLoc2)-2.0d0*DDg(iLoc3)-4.0d0*DDg(iLoc5)
C         Vlkmn=DDg(L,N,L,L)+DDg(L,N,N,N)-2.0d0*DDg(L,N,K,K)
C     &         -4.0d0*DDg(K,L,M,N)
        END IF
        IF(L.eq.N) THEN
         iLoc4=K+(M-1)*lRoots
         iLoc1=iLoc4+(K-1)*lRoots23
         iLoc2=iLoc4+(M-1)*lRoots23
         iLoc3=iLoc4+(L-1)*lRoots23
         Vklnm=DDg(iLoc1)+DDg(iLoc2)-2.0d0*DDg(iLoc3)-4.0d0*DDg(iLoc5)
C         Vklnm=DDg(K,M,K,K)+DDg(K,M,M,M)-2.0d0*DDg(K,M,L,L)
C     &         -4.0d0*DDg(K,L,M,N)
        END IF
        Hess((iKL-1)*nSPair+iMN)=Vklmn+Vlknm-Vklnm-Vlkmn
       End Do
       End Do
      END DO
      END DO


      RETURN
      End Subroutine

************************************************************************

      Subroutine CalcGradCMS(Grad,DDg,nDDg,lRoots,nSPair)
      INTEGER nDDg,lRoots,nSPair
      Real*8 Grad(nSPair),DDg(nDDg)

      INTEGER K,L,iKL,lRoots2,lRoots3,iLoc1,iLoc2

      lRoots2=lRoots**2
      lRoots3=lRoots*lRoots2

      DO K=2,lRoots
       Do L=1,K-1
        iLoc1=K+(K-1)*lRoots+(K-1)*lRoots2+(L-1)*lRoots3
        iLoc2=L+(L-1)*lRoots+(K-1)*lRoots2+(L-1)*lRoots3
        iKL=(K-1)*(K-2)/2+L
        Grad(iKL)=DDg(iLoc1)-DDg(iLoc2)
       End Do
      END DO
      CALL DSCal_(nSPair,2.0d0,Grad,1)
      RETURN
      End Subroutine
************************************************************************


************************************************************************
      Subroutine CalcQaa(Qaa,DDg,lRoots,nDDg)
      INTEGER lRoots,nDDg
      Real*8 DDg(nDDg)
      Real*8 Qaa

      INTEGER iState,iLoc,Int1,lRoots2

      lRoots2=lRoots**2
      Int1=(lRoots2+1)*(lRoots+1)
      Qaa=0.0d0
      DO iState=1,lRoots
       iLoc=(iState-1)*Int1+1
       Qaa=Qaa+DDg(iLoc)
      END DO
      Qaa=Qaa/2.0d0

      RETURN
      End Subroutine
