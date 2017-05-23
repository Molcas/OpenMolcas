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
      Subroutine LevMarquart(iPotte,nPick,ipPick,nEPP,ipEPCo,Coo
     &                      ,dMullig,lMax,A,iAtom,jAtom,chP
     &                      ,Thrs1,Thrs2,nThrs,Chi2B,iPrint,AboveMul)
      Implicit real*8 (a-h,o-z)

#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension Coo(3),A(2),B(2),dMullig((lMax*(lMax**2+6*lMax+11)+6)/6)
      Dimension AlfMat(4),AlfMatI(4),rStore(nPick),rinvStore(nPick)
      Dimension ARaw(2,2),BRaw(2),dA(2),Pout(nPick)
      Dimension xStore(nPick),yStore(nPick),zStore(nPick)

      Logical AboveMul(2)
      Logical lStop1,lStop2,lStop3,lStop4
      Logical lScreen1,lScreen2,lScreen3,lScreen4

      Character*60 UtChar

*
*-- Set iteration count to zero here at the top of the
*   Levenberg-Marquart thing. Initiate lambda and the exponents.
*   The exponents are not God-given, but these rather low values
*   have been found to work nicely in most cases. Set highest
*   and lowest allowed values of exponents. Set highest and lowest
*   allowed values of exponent steps.
*
      Iter=0
      nStep=0
      dLambda=1d-3
      A(1)=0.5d0
      A(2)=0.5d0
      dUpper=3.0d0
      dLower=1d-1
      ddUpper=0.7d0
      ddLower=-0.7d0

*
*-- Start loop, set zeros.
*
9901  Continue
      Iter=Iter+1
      Chi2=0.0d0
      Chi2B=0.0d0
      Do i=1,2
        BRaw(i)=0.0d0
        dA(i)=0.0d0
        Do j=1,2
          ARaw(i,j)=0.0d0
        Enddo
      Enddo

*
*-- Compute distances and store them.
*
      Do iP=1,nPick
        ind=iWork(ipPick+iP-1)
        x=Work(ipEPCo+(ind-1)*3+0)-Coo(1)
        y=Work(ipEPCo+(ind-1)*3+1)-Coo(2)
        z=Work(ipEPCo+(ind-1)*3+2)-Coo(3)
        xStore(iP)=x
        yStore(iP)=y
        zStore(iP)=z
        rStore(iP)=sqrt(x**2+y**2+z**2)
        rinvStore(iP)=1.0d0/rStore(iP)
      Enddo

*
*-- Compute the electric potential and the difference between true and
*   modelled potential. Also compute first derivatives. Formaly, there
*   are second derivatives as well, but they are neglected for
*   stability reasons. For this relatively simple non-linear problem,
*   it is probably not meaningful to code a more intelligent
*   routine that includes the second derivatives with some
*   stabalizing modifications.
*
      DerMax=0d0
      Do iP=1,nPick
        r=rStore(iP)
        rinv=rinvStore(iP)
        x=xStore(iP)
        y=yStore(iP)
        z=zStore(iP)
        Potte=ElPot(r,rinv,x,y,z,dMullig,lMax,A,chP,.true.,.true.)
        Diffo=Work(iPotte+iP-1)-Potte
        Der1=dMullig(1)*(1.0d0+2.0d0*A(1)*r)*exp(-2.0d0*A(1)*r)
        Der2=(dMullig(2)*x+dMullig(3)*y+dMullig(4)*z)
     &      *(A(2)**2+2.0d0*A(2)**3*r)*exp(-2.0d0*A(2)*r)
        If(abs(Der1).gt.DerMax)Dermax=abs(Der1)
        If(abs(Der2).gt.DerMax)Dermax=abs(Der2)
*
*---- Assemble the alpha-matrix, the beta-vector and the chi2 scalar.
*
        Chi2=Chi2+Diffo**2
        BRaw(1)=BRaw(1)+Diffo*Der1
        BRaw(2)=BRaw(2)+Diffo*Der2
        ARaw(1,1)=ARaw(1,1)+Der1*Der1
        ARaw(2,1)=ARaw(2,1)+Der2*Der1
        ARaw(2,2)=ARaw(2,2)+Der2*Der2
      Enddo
      Chi2=Chi2/dble(nPick)

*
*-- Check if numerical and statistical well-behaved calculations are
*   expected or not.
*
      If(nPick.lt.5) then
        Write(6,*)
        Write(6,'(A)')' Error in numerical fit for exponent!'
        Write(6,'(A)')'    Too few electric potential points have'
        Write(6,'(A,2I4)')'    been sampled for (At1,At2):'
     &                              ,iAtom,jAtom
        Write(6,*)
        Call Quit(_RC_GENERAL_ERROR_)
      Endif
      If(DerMax.lt.1d-4) then
        Write(6,*)
        Write(6,'(A)')' Error in numerical fit for exponent!'
        Write(6,'(A)')'    To small derivative for parameter'
        Write(6,'(A,2I4)')'    in centre (At1,At2):'
     &                              ,iAtom,jAtom
        Write(6,'(A,E12.4)')'    Maximal derivative:',DerMax
        Write(6,*)
        Write(6,'(A)')'    Either include closer points, give a'
        Write(6,'(A)')'    better initial estimate, or increase'
        Write(6,'(A)')'    multipole magnitude threshold.'
        Write(6,*)
        Call Quit(_RC_GENERAL_ERROR_)
      Endif

*
*-- Solve the linear system to get step, dA.
*
      Call SolveA(AlfMat,AlfMatI,dLambda,dMullig,lMax
     &           ,ARaw,BRaw,dA,iPrint,AboveMul,ddUpper,ddLower)

*
*-- Make a screening of dA: if the maximal or minimal exponent
*   has been reached, but the optimization still pushes on, screen
*   them by putting them as zero, and modifying abovemul.
*
      Do i=1,2
        If(AboveMul(i)) then
          lScreen1=(A(i)+1d-8).gt.dUpper
          lScreen2=(dA(i)+1d-8).gt.ddUpper
          lScreen3=(A(i)-1d-8).lt.dLower
          lScreen4=(dA(i)-1d-8).lt.ddLower
          If((lScreen1.and.lScreen2).or.(lScreen3.and.lScreen4)) then
            dA(i)=0.0d0
            AboveMul(i)=.false.
          Endif
        Endif
      Enddo
      If(iPrint.ge.7) then
        Call RecPrt('deltaA',' ',dA,2,1)
      Endif

*
*-- Construct trial parameters and compute the chi2-scalar for
*   the trial parameters. To keep things within reasonable bounds
*   do not exceed certain values.
*
      B(1)=A(1)+dA(1)
      B(2)=A(2)+dA(2)
      If(B(1).lt.dLower)B(1)=dLower
      If(B(2).lt.dLower)B(2)=dLower
      If(B(1).gt.dUpper)B(1)=dUpper
      If(B(2).gt.dUpper)B(2)=dUpper
      Do iP=1,nPick
        r=rStore(iP)
        rinv=rinvStore(iP)
        x=xStore(iP)
        y=yStore(iP)
        z=zStore(iP)
        Potte=ElPot(r,rinv,x,y,z,dMullig,lMax,B,chP,.true.,.true.)
        Diffo=Work(iPotte+iP-1)-Potte
        Chi2B=Chi2B+Diffo**2
      Enddo
      Chi2B=Chi2B/dble(nPick)

*
*-- Take appropriate action given the different Stop-criterion. They
*   measure the following: (1) The error should not change too much
*   between steps; observe that it is not meaningful to have a too
*   tight threshold in this regard due to the statistical nature of
*   the problem. (2) The modifier parameter has to stabalize and not
*   be far out in the linear region, rather in the second order regime.
*   (3) The last step should be a decrease. (4) A certain number of
*   steps should preceed that decreases the error; this threshold
*   'overlaps' some with second threshold. Halt the optimization if
*   no convergence is reached after 40 iterations. This is by far
*   a generous limit.
*
      lStop1=abs(Chi2-Chi2B).lt.Thrs1
      lStop2=dLambda.lt.Thrs2
      lStop3=Chi2.gt.Chi2B
      lStop4=nStep.ge.nThrs
      If(iPrint.ge.6) then
        Write(6,*)
        Write(6,790)Iter
        Write(6,794)'Chi2','Chi2, trial','Lambda'
        Write(6,791)Chi2,Chi2B,dLambda
        Write(6,795)'Exponents:','Charge','Dipole'
        Write(6,792)A(1),A(2)
      Endif
      If(lStop1.and.lStop2.and.lStop3.and.lStop4) then
        A(1)=B(1)
        A(2)=B(2)
        Go To 9902
      Endif
      If(Iter.gt.40) then
        Write(6,*)
        Write(6,*)'No convergence in the Levenberg-Marquart section.'
        Call Quit(_RC_NOT_CONVERGED_)
      Endif
      If(Chi2B.ge.Chi2) then
        dLambda=dLambda*10.0d0
        Go To 9901
      Else
        dLambda=dLambda*0.1d0
        nStep=nStep+1
        A(1)=B(1)
        A(2)=B(2)
        Go To 9901
      Endif

9902  Continue

*
*-- Optional printing when convergence is reached.
*
      If(iPrint.ge.5) then
        Write(6,*)
        Write(6,*)
        Write(6,793)Iter
        Write(6,794)'Chi2','Chi2, trial','Lambda'
        Write(6,791)Chi2,Chi2B,dLambda
        Write(6,795)'Exponents:','Charge','Dipole'
        Write(6,792)A(1),A(2)
        Do iP=1,nPick
          r=rStore(iP)
          rinv=rinvStore(iP)
          x=xStore(iP)
          y=yStore(iP)
          z=zStore(iP)
          Pout(iP)=ElPot(r,rinv,x,y,z,dMullig,lMax,A,chP,.true.,.true.)
        Enddo
        Write(UtChar,'(A,2I3)')'Approximate partial density potential,'
     &                       //' centre',iAtom,jAtom
        Call RecPrt(UtChar,' ',Pout,nPick,1)
        Call RecPrt('Distance to points',' ',rStore,nPick,1)
      Endif

790   Format('Levenberg-Marquart optimization. Iteration:',I3)
794   Format('   ',A,'          ',A,'   ',A)
791   Format(3(E13.6,' '))
795   Format(' ',A,' ',A,'      ',A)
792   Format('       ',2F12.6)
793   Format('Convergence reached in iteration ',I2)

*
*-- And yes, there is an end here as well.
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nEPP)
      End


*
*-- The electric potential with diffuse s- and p-functions. No
*   higher than d-functions (with non-zero trace).
*
      real*8 Function ElPot(r,rinv,x,y,z,dMullig,lMax,A,chP
     &                               ,lDOrNot1,lDOrNot2)
      Implicit real*8 (a-h,o-z)

#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension A(2),dMullig((lMax*(lMax**2+6*lMax+11)+6)/6)
      Dimension dL2(6),dL3(10),dL4(15),dL5(21)

      Logical lDOrNot1,lDOrNot2

      ElPot=0.0d0
      If(lMax.ge.0) then
        If(lDOrNot1) then
          dCh=chP*rinv
          dCh=dCh+dMullig(1)*rinv*(1.0d0-(1.0d0+A(1)*r)
     &                          *exp(-2.0d0*A(1)*r))
          ElPot=dCh
        Else
          dL2(1)=dMullig(1)+chP
          ElPot=ElPot+ElPointPot(rinv,x,y,z,0,dL2)
        Endif
      Endif
      If(lMax.ge.1) then
        If(lDOrNot2) then
          ar=A(2)*r
          dM=(x*dMullig(2)+y*dMullig(3)+z*dMullig(4))
     &       *rinv**3*(1.0d0-(1.0d0+2.0d0*ar+2.0d0*ar**2+ar**3)
     &       *exp(-2.0d0*ar))
          ElPot=ElPot+dM
        Else
          dL2(1)=dMullig(2)
          dL2(2)=dMullig(3)
          dL2(3)=dMullig(4)
          ElPot=ElPot+ElPointPot(rinv,x,y,z,1,dL2)
        Endif
      Endif
      If(lMax.ge.2) then
        dL2(1)=dMullig(5)
        dL2(2)=dMullig(6)
        dL2(3)=dMullig(7)
        dL2(4)=dMullig(8)
        dL2(5)=dMullig(9)
        dL2(6)=dMullig(10)
        ElPot=ElPot+ElPointPot(rinv,x,y,z,2,dL2)
      Endif
      If(lMax.ge.3) then
        dL3(1)=dMullig(11)
        dL3(2)=dMullig(12)
        dL3(3)=dMullig(13)
        dL3(4)=dMullig(14)
        dL3(5)=dMullig(15)
        dL3(6)=dMullig(16)
        dL3(7)=dMullig(17)
        dL3(8)=dMullig(18)
        dL3(9)=dMullig(19)
        dL3(10)=dMullig(20)
        ElPot=ElPot+ElPointPot(rinv,x,y,z,3,dL3)
      Endif
      If(lMax.ge.4) then
        dL4(1)=dMullig(21)
        dL4(2)=dMullig(22)
        dL4(3)=dMullig(23)
        dL4(4)=dMullig(24)
        dL4(5)=dMullig(25)
        dL4(6)=dMullig(26)
        dL4(7)=dMullig(27)
        dL4(8)=dMullig(28)
        dL4(9)=dMullig(29)
        dL4(10)=dMullig(30)
        dL4(11)=dMullig(31)
        dL4(12)=dMullig(32)
        dL4(13)=dMullig(33)
        dL4(14)=dMullig(34)
        dL4(15)=dMullig(35)
        ElPot=ElPot+ElpointPot(rinv,x,y,z,4,dL4)
      Endif
      If(lMax.ge.5) then
        dL5(1)=dMullig(36)
        dL5(2)=dMullig(37)
        dL5(3)=dMullig(38)
        dL5(4)=dMullig(39)
        dL5(5)=dMullig(40)
        dL5(6)=dMullig(41)
        dL5(7)=dMullig(42)
        dL5(8)=dMullig(43)
        dL5(9)=dMullig(44)
        dL5(10)=dMullig(45)
        dL5(11)=dMullig(46)
        dL5(12)=dMullig(47)
        dL5(13)=dMullig(48)
        dL5(14)=dMullig(49)
        dL5(15)=dMullig(50)
        dL5(16)=dMullig(51)
        dL5(17)=dMullig(52)
        dL5(18)=dMullig(53)
        dL5(19)=dMullig(54)
        dL5(20)=dMullig(55)
        dL5(21)=dMullig(56)
        ElPot=ElPot+ElpointPot(rinv,x,y,z,5,dL5)
      Endif
      If(lMax.ge.6) then
        Write(6,*)
        Write(6,*)'Oops! You hit the roof with respect to angular'
     &          //' momentum. Lower that, or do some programming.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif

      Return
      End
