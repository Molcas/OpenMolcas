!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Do_PI2(D1mo,nd1mo,TabMO,mAO,mGrid,nMOs,                &
     &   P2_ontop,nP2_ontop,RhoI,RhoA,mRho,Do_Grad,                     &
     &   P2MOCube,MOs,MOx,MOy,MOz)
!***********************************************************************
!                                                                      *
! Object: Calculation P2 ontop density and its derivatives             *
!                                                                      *
! Called from: Do_batch                                                *
!                                                                      *
! Calling    : FZero                                                   *
!                                                                      *
!    INPUT:                                                            *
!   D1mo     = one-body density matrix in MO basis                     *
!   nd1mo    = size of D1mo                                            *
!   TabMO    = MO values computed on grid                              *
!   nMOs     = number of MO basis                                      *
!   mAO      = number of derivatives of AO...                          *
!   mGrid    = number of grid points                                   *
!                                                                      *
!***********************************************************************
      use nq_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
!#include "stdalloc.fh"
      Real*8 D1mo(nd1mo),TabMO(mAO,mGrid,nMOs),                         &
     &       P2_ontop(nP2_ontop,mGrid)
      Real*8 RhoI(mRho,mGrid)
      Real*8 RhoA(mRho,mGrid)
      Logical Do_Grad

      INTEGER IOff1,IOff2
      REAL*8,DIMENSION(mGrid*NASHT)::P2MOCube,MOs,MOx,MOy,MOz
      Real*8 ddot_
      External DDot_
!***********************************************************************
!                                                                      *
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!                                                                      *
!***********************************************************************
      If (nP2_ontop.eq.4) Then
         If (mAO.ne.4.or.mRho.ne.4) Then
           Call WarningMessage(2,' Somthings wrong in dim. in p2cs')
           Call Abend()
       End If
      Else If (nP2_ontop.eq.6) Then
         If (mAO.ne.10.or.mRho.ne.6) Then
           Call WarningMessage(2,' Somthings wrong in dim. in p2cs')
           Call Abend()
      End If
      End If
!
      Call FZero(P2_ontop,mGrid*nP2_ontop)
      jOffA_ = 0
      jOffB_ = 0
      Do iIrrep = 0, mIrrep-1
         iOff_Ash(iIrrep)=jOffA_
         iOff_Bas(iIrrep)=jOffB_
         iOff_BasAct(iIrrep)=jOffB_ + nIsh(iIrrep) + nFro(iIrrep)
         jOffA_=jOffA_+nAsh(iIrrep)
         jOffB_=jOffB_+mBas(iIrrep)
      End Do
!***********************************************************************
!                                                                      *
!   P(1,...) - P_2                                                     *
!   P(2,...), P(3,...), P(4,...) - grad P_2                            *
!   Not implemented:                                                   *
!   P(5,...) - grad^2 P_2                                              *
!   P(6,...) - additional part grad^2 P_2 for CS functional            *
!   P(5) and P(6) removed                                              *
!***********************************************************************

!***********************************************************************
!          Inactive part:                                              *
!***********************************************************************
      NumIsh = 0
      NumAsh = 0
      Do iIrrep=0, mIrrep-1
         NumIsh = NumIsh + nISh(iIrrep)
         NumAsh = NumAsh + nAsh(iIrrep)
      End Do
!
      Do iGrid = 1, mGrid
      Do iIrrep=0, mIrrep-1
!      Write(6,*) " Symm:",iIrrep
       Do i_=1,nISh(iIrrep) + nFro(iIrrep)
        i = iOff_Bas(iIrrep) + i_
!
!         Write(6,*)  " do_p2: Inact-Inact:", iIrrep,i,
!     &               TabMO(1,iGrid,i)
!
          RhoI(1,iGrid) = RhoI(1,iGrid) +                               &
     &                    TabMO(1,iGrid,i)*TabMO(1,iGrid,i)
!      write(6,'(A15,2I3,2G15.8)')'iGrid,i,MO,RhoI',
!     &   iGrid,i,TabMO(1,iGrid,i), RhoI(1,iGrid)
!        if (Functional_type.eq.GGA_type.or.Do_Grad) then
        if (Functional_type.eq.GGA_type) then
          RhoI(2,iGrid) = RhoI(2,iGrid) +                               &
     &                    TabMO(1,iGrid,i)*TabMO(2,iGrid,i)
          RhoI(3,iGrid) = RhoI(3,iGrid) +                               &
     &                    TabMO(1,iGrid,i)*TabMO(3,iGrid,i)
          RhoI(4,iGrid) = RhoI(4,iGrid) +                               &
     &                    TabMO(1,iGrid,i)*TabMO(4,iGrid,i)
        end if
!
       End Do         ! i_
      End Do         ! iIrrep
      End Do         ! iGrid
!
      If (NumIsh.ne.0) Then
      Do iGrid = 1, mGrid
         P2_ontop(1,iGrid) = RhoI(1,iGrid)*RhoI(1,iGrid)
!         Write(6,'(A15,I3,1G28.20)')'iGrid,P2(1)=',
!     *          iGrid,P2_ontop(1,iGrid)
!
!        if (Functional_type.eq.GGA_type.or.Do_Grad) then
        if (Functional_type.eq.GGA_type) then
            P2_ontop(2,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(2,iGrid)
!            Write(6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
            P2_ontop(3,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(3,iGrid)
!            Write(6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
            P2_ontop(4,iGrid) = 4.0d0*RhoI(1,iGrid)*RhoI(4,iGrid)
!            Write(6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
         end if
        if (Functional_type.eq.LDA_type.and.Do_Grad) then
!Here I must
!1. transform the 2-body density matrix to AO

!2. Loop over effective gradients
!3. Calculate P2_ontop_d(eff_Grad,iGrid)


        end if
      End Do
      End If


!***********************************************************************
!          Active-Inactive part:                                       *
!***********************************************************************
      If (NumIsh.ne.0.and.NumAsh.ne.0) Then
       Do kIrrep = 0, mIrrep-1
         Do k_ = 1, nASh(kIrrep)
            k= k_ + iOff_BasAct(kIrrep)
            Do lIrrep = 0, mIrrep-1
               Do l_ = 1, nAsh(lIrrep)
                  l= l_ + iOff_BasAct(lIrrep)
                  kl=iTri(k_ + iOff_Ash(kIrrep) ,                       &
     &                    l_ + iOff_Ash(lIrrep) )
                  Do iGrid = 1, mGrid
                     RhoA(1,iGrid) = RhoA(1,iGrid) +                    &
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(1,iGrid,l)
!       Write(6,'(A35,3I3,3G15.8)') 'iGrid,k,l,D1mo(kl),Tab(k),Tab(l)=',
!     &    iGrid,k,l,D1mo(kl),TabMO(1,iGrid,k),TabMO(1,iGrid,l)
!                     if (Functional_type.eq.GGA_type.or.Do_Grad) Then
                     if (Functional_type.eq.GGA_type) Then
                     RhoA(2,iGrid) = RhoA(2,iGrid) +                    &
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(2,iGrid,l)
                     RhoA(3,iGrid) = RhoA(3,iGrid) +                    &
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(3,iGrid,l)
!        write(6,*) 'RhoA(4,iGrid) bf =', RhoA(4,iGrid)
                     RhoA(4,iGrid) = RhoA(4,iGrid) +                    &
     &                   D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
!       Write(6,*) 'D1mo(kl),Tab(1,k),Tab(1,l)=',
!     &    D1mo(kl)*TabMO(1,iGrid,k)*TabMO(4,iGrid,l)
                     end if
                  End Do     ! iGrid
               End Do      ! l_
            End Do       ! lIrrep
         End Do        ! k_
       End Do        ! kIrrep
!
       Do iGrid = 1, mGrid
               P2_ontop(1,iGrid) = P2_ontop(1,iGrid) +                  &
     &                           RhoI(1,iGrid)*RhoA(1,iGrid)
        if (Functional_type.eq.GGA_type) Then
               P2_ontop(2,iGrid) = P2_ontop(2,iGrid) +                  &
     &                     2.0d0*RhoI(2,iGrid)*RhoA(1,iGrid) +          &
     &                     2.0d0*RhoI(1,iGrid)*RhoA(2,iGrid)
!            Write(6,'(A,1f28.20)') 'P2(2)   =',P2_ontop(2,iGrid)
               P2_ontop(3,iGrid) = P2_ontop(3,iGrid) +                  &
     &                     2.0d0*RhoI(3,iGrid)*RhoA(1,iGrid) +          &
     &                     2.0d0*RhoI(1,iGrid)*RhoA(3,iGrid)
!            Write(6,'(A,1f28.20)') 'P2(3)   =',P2_ontop(3,iGrid)
               P2_ontop(4,iGrid) = P2_ontop(4,iGrid) +                  &
     &                     2.0d0*RhoI(4,iGrid)*RhoA(1,iGrid) +          &
     &                     2.0d0*RhoI(1,iGrid)*RhoA(4,iGrid)
!            Write(6,'(A,1f28.20)') 'P2(4)   =',P2_ontop(4,iGrid)
        end if
       End Do ! loop over grid points
      End If ! if Inactive
!***********************************************************************
!
!          Active-Active part:
!
!***********************************************************************


      IF (NumAsh.eq.0) RETURN

!      write(6,*) 'P2MOCube in do_pi2'
!      CALL RecPrt(' ',' ',P2MOCube,NASHT,mGrid)
!
!      write(6,*) 'MOs array in do_pi2'
!      CALL RecPrt(' ',' ',MOs,NASHT,mGrid)

      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT
       Do kIrrep=0,mIrrep-1
        IOff2=IOff1+iOff_Ash(kIrrep)+1
        P2_ontop(1,iGrid)=P2_ontop(1,iGrid)+                            &
     &  ddot_(nAsh(kIrrep),MOs(IOff2),1,P2MOCube(IOff2),1)
       End Do
      END DO

      IF(Functional_type.eq.GGA_type) THEN
       DO iGrid=1,mGrid
        IOff1=(iGrid-1)*NASHT
        Do kIrrep=0,mIrrep-1
         IOff2=IOff1+iOff_Ash(kIrrep)+1
         P2_ontop(2,iGrid)=P2_ontop(2,iGrid)+                           &
     &4.0d0*ddot_(nAsh(kIrrep),MOx(IOff2),1,P2MOCube(IOff2),1)
         P2_ontop(3,iGrid)=P2_ontop(3,iGrid)+                           &
     &4.0d0*ddot_(nAsh(kIrrep),MOy(IOff2),1,P2MOCube(IOff2),1)
         P2_ontop(4,iGrid)=P2_ontop(4,iGrid)+                           &
     &4.0d0*ddot_(nAsh(kIrrep),MOz(IOff2),1,P2MOCube(IOff2),1)
        End Do
       END DO
      END IF

!      write(6,*) 'On-top density new code'
!C      write(6,'(10(F9.6,1X))')(P2_Ontop(1,iGrid),iGrid=1,mGrid)
!      write(6,*)(P2_Ontop(1,iGrid),iGrid=1,mGrid)
      RETURN
      END



