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

subroutine AODKHEXP(NBasis,M,MX,M0,Ep,E0,EL,ES,OL,W,TA,TB,A1,A2,O2,E2,F2,WS)
! Do the Douglas-Kroll-Hess Arbitrary order term
!
! Storage requirements:
! Ep,E0 - NBasis
! EL,ES,OL - NBasis*NBasis
! W,TA,TB,A1,A2 - NBasis*NBasis
! O2,E2,F2 - Nbasis*Nbasis*M (M=order of DKH)

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBasis, M, MX, M0
real(kind=wp), intent(in) :: Ep(NBasis), E0(NBasis), ES(NBasis,NBasis), OL(NBasis,NBasis)
real(kind=wp), intent(inout) :: EL(NBasis,NBasis), TA(NBasis,NBasis), TB(NBasis,NBasis), A1(NBasis,NBasis), A2(NBasis,NBasis), &
                                O2(NBasis,NBasis,M), E2(NBasis,NBasis,M), F2(NBasis,NBasis,M), WS(NBasis,NBasis,M)
real(kind=wp), intent(out) :: W(NBasis,NBasis)
integer(kind=iwp) :: I, Ioe, Iut, J, K, KK, Ks
real(kind=wp) :: Cof
logical(kind=iwp) :: Ifodd

! Copy 1st order DKH effective potential
! E2(upper-left) O2(upper-right) F2(lower-right)

do J=1,NBasis
  do I=1,NBasis
    E2(I,J,1) = EL(I,J)
    F2(I,J,1) = ES(I,J)
    O2(I,J,1) = OL(I,J)
  end do
end do

! Apply [M/2] times Douglas-Kroll-Hess Unitary transformations

do Iut=1,M/2
  do I=1,NBasis
    do J=1,NBasis
      W(J,I) = O2(J,I,Iut)/(Ep(I)+Ep(J))
      if (Iut <= MX) WS(J,I,Iut*2-1) = W(J,I)
    end do
  end do
  ! Apply W_{Iut} to O/E_{Ks}
  do Ks=M-Iut,1,-1
    ! W_{1} only apply to O/E_{1}
    if (Iut == 1 .and. Ks >= 2) goto 40
    do Ioe=1,2
      ! O_{k,k<Iut} was eliminated
      if (Ioe == 1 .and. Ks < Iut) goto 50
      Ifodd = Ioe == 1
      do I=1,NBasis
        do J=1,NBasis
          ! Copy O/E_{Ks} terms to temp arrays
          if (Ifodd) then
            TA(J,I) = O2(J,I,Ks)
          else
            TA(J,I) = E2(J,I,Ks)
            TB(J,I) = F2(J,I,Ks)
          end if
        end do
      end do
      do K=Ks,M,Iut
        ! skip terms do not contribute to final DKH Hamiltonian (even,upper-left)
        if (K+Iut+Iut > M .and. (K+Iut > M .or. .not. Ifodd)) goto 60
        KK = (K-Ks)/Iut+1
        if (Ioe == 1 .and. Ks == Iut) then
          ! see Eq.(74) of JCP130(2009)044102
          if (KK == 1) then
            Cof = Half
          else
            Cof = real(KK,kind=wp)/(KK*KK-One)
          end if
        else
          ! see Eq.(71) of JCP130(2009)044102
          Cof = One/KK
        end if
        if (Ifodd) then
          ! skip terms do not contribute to final DKH Hamiltonian (even,upper-left)
          if (K+Iut+Iut+Iut <= M) then
            call DGEMM_('T','N',NBasis,NBasis,NBasis,Cof,W,NBasis,TA,NBasis,Zero,A2,NBasis)
          end if
          call DGEMM_('N','T',NBasis,NBasis,NBasis,Cof,W,NBasis,TA,NBasis,Zero,A1,NBasis)
          ! ( 0  W)(0  O)   (0  O)( 0  W)   ( WO'+(WO')'     0       )
          ! (-W' 0)(O' 0) - (O' 0)(-W' 0) = (    0       -W'O-(W'O)' )
          !  where W'=W^{\dag} O'=O^{\dag}
          do I=1,NBasis
            do J=1,NBasis
              if (K+Iut+Iut+Iut <= M) then
                TB(J,I) = -A2(J,I)-A2(I,J)
                F2(J,I,K+Iut) = F2(J,I,K+Iut)+TB(J,I)
              end if
              TA(J,I) = A1(J,I)+A1(I,J)
              E2(J,I,K+Iut) = E2(J,I,K+Iut)+TA(J,I)
            end do
          end do
        else
          call DGEMM_('N','N',NBasis,NBasis,NBasis,Cof,W,NBasis,TB,NBasis,Zero,A1,NBasis)
          call DGEMM_('N','N',NBasis,NBasis,NBasis,Cof,TA,NBasis,W,NBasis,Zero,A2,NBasis)
          ! ( 0  W)(E 0)   (E 0)( 0  W)   (    0     WF-EW )
          ! (-W' 0)(0 F) - (0 F)(-W' 0) = ( (WF-EW)'   0   )
          !  where W'=W^{\dag}
          do I=1,NBasis
            do J=1,NBasis
              TA(J,I) = A1(J,I)-A2(J,I)
              O2(J,I,K+Iut) = O2(J,I,K+Iut)+TA(J,I)
            end do
          end do
        end if
        Ifodd = .not. Ifodd
60      continue
      end do
50    continue
    end do
40  continue
  end do
end do

! Sum all even terms to Douglas-Kroll-Hess Hamiltonian

do I=1,NBasis
  EL(I,I) = EL(I,I)+E0(I)
end do
do I=2,M0
  do J=1,NBasis
    do K=1,Nbasis
      EL(K,J) = EL(K,J)+E2(K,J,I)
    end do
  end do
end do

return

end subroutine AODKHEXP
