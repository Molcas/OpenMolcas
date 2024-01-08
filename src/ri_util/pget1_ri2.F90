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
! Copyright (C) 1992,2007, Roland Lindh                                *
!               2009, Francesco Aquilante                              *
!***********************************************************************

subroutine PGet1_RI2(PAO,ijkl,nPAO,iCmp,iAO,iAOst,jBas,lBas,kOp,ExFac,CoulFac,PMax,V_K,U_K,mV_K,Z_p_K,nSA)
!***********************************************************************
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density.                        *
!                                                                      *
!          (Only for use with C1 point group symmetry)                 *
!                                                                      *
!          The indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!                                                                      *
!             Modified for RI-DFT, March 2007                          *
!                                                                      *
!             Modified for RI-HF/CAS, Dec 2009 (F. Aquilante)          *
!***********************************************************************

use Symmetry_Info, only: Mul
use Basis_Info, only: nBas
use SOAO_Info, only: iAOtSO
use pso_stuff, only: A_PT2, DMdiag, Gamma_On, lPSO, lSA, nnP, nPos
use RI_glob, only: A, AMP2, CijK, iAdrCVec, iMP2prpt, iUHF, LuAVector, LuCVector, nAuxVe, nChOrb, nIJ1, nKvec, tavec
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ijkl, nPAO, iCmp(4), iAO(4), iAOst(4), jBas, lBas, kOp(4), mV_K, nSA
real(kind=wp), intent(out) :: PAO(ijkl,nPAO), PMax
real(kind=wp), intent(in) :: ExFac, CoulFac, V_K(mV_K,nSA), U_K(mV_K), Z_p_K(nnP(0),mV_K,*)
integer(kind=iwp) :: i, i2, i4, iAdrA, iAdrJ, iAdrL, iE, iOffA, iPAO, iS, iSO, iSO2, iSym, j, jAOj, jik, jil, jp, jSO, jSOj, jSym, &
                     k, kl, kSym, l, lAOl, lSO, lSOl, lSym, lTot, n, nijkl, nik, nik1, nik2
real(kind=wp) :: Cpu, Cpu1, Cpu2, Fac, Factor, temp, temp2, tempJ_mp2, tempK_mp2, tmp, Wall, Wall1, Wall2
real(kind=wp), pointer :: CiKj(:,:), CiKl(:), V2(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
do i=1,nSA
  call RecPrt('PGet1_RI2: V_k',' ',V_k(1,i),1,mV_k)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! DeSymP will treat up to eight fold degeneracy due to permutational
! symmetry of shell quadruplets. We will have to compensate for that
! here since we only have shell doublets.

call CWTime(Cpu1,Wall1)

if (min(lBas,jBas) == 0) return

Fac = Quart
PMax = Zero
iPAO = 0
iOffA = nBas(0)

!                                                                      *
!***********************************************************************
!                                                                      *

if (ExFac == Zero) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Pure DFT

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)

    do i4=1,iCmp(4)

      lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

      iPAO = iPAO+1
      nijkl = 0

      do lAOl=0,lBas-1
        lSOl = lSO+lAOl-iOffA
        do jAOj=0,jBas-1

          jSOj = jSO+jAOj-iOffA
          nijkl = nijkl+1

          ! Coulomb contribution
          temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
          !temp = Zero

          PMax = max(PMax,abs(temp))
          PAO(nijkl,iPAO) = Fac*temp

        end do
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((iMP2prpt /= 2) .and. (.not. lPSO) .and. (iUHF == 0)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Hybrid DFT and HF

  iSO = 1

  jSym = 1
  kSym = jSym
  iSym = 1
  lSym = Mul(jSym,iSym)

  nik = nIJ1(iSym,kSym,iSO)

  n = nik*jBas
  iS = 1
  iE = n
  CiKj(1:n,1:1) => CijK(iS:iE)
  n = nik*lBas
  iS = iE+1
  iE = iE+n
  CiKl(1:n) => CijK(iS:iE)

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)

    ! Pick up the MO transformed fitting coefficients, C_ik^J
    jSOj = jSO-iOffA
    iAdrJ = nik*(jSOj-1)+iAdrCVec(jSym,iSym,1)
    call dDaFile(LuCVector(jSym,1),2,CikJ(:,1),nik*jBas,iAdrJ)

    do i4=1,iCmp(4)
      lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

      iPAO = iPAO+1
      nijkl = 0

      if (lSO /= jSO) then
        lSOl = lSO-iOffA
        iAdrL = nik*(lSOl-1)+iAdrCVec(jSym,iSym,1)
        call dDaFile(LuCVector(jSym,1),2,CiKl,nik*lBas,iAdrL)

        V2(1:) => CiKl(1:)
      else
        V2(1:) => CiKj(1:,1)
      end if

      A(1:jBas*lBas) = Zero
      call DGEMM_('T','N',jBas,lBas,nik,One,CiKj,nik,V2,nik,Zero,A,jBas)

      do lAOl=0,lBas-1
        lSOl = lSO+lAOl-iOffA
        do jAOj=0,jBas-1
          jSOj = jSO+jAOj-iOffA
          nijkl = nijkl+1

          temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
          temp = temp-ExFac*A(nijkl)

          PMax = max(PMax,abs(temp))
          PAO(nijkl,iPAO) = Fac*temp
        end do
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((iMP2prpt /= 2) .and. (.not. lPSO) .and. (iUHF == 1)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Hybrid UDFT and UHF

  jSym = 1
  kSym = jSym
  iSym = 1
  lSym = Mul(jSym,iSym)
  nik1 = nIJ1(iSym,kSym,1)
  nik2 = nIJ1(iSym,kSym,2)
  nik = max(nik1,nik2)

  n = nik*jBas
  iS = 1
  iE = n*2
  CiKj(1:n,1:2) => CijK(iS:iE)

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    jSOj = jSO-iOffA

    ! Pick up the MO transformed fitting coefficients, C_ik^J
    if (nik1 /= 0) then
      iAdrJ = nik1*(jSOj-1)+iAdrCVec(jSym,iSym,1)
      call dDaFile(LuCVector(jSym,1),2,CiKj(:,1),nik1*jBas,iAdrJ)
    end if

    if (nik2 /= 0) then
      iAdrJ = nik2*(jSOj-1)+iAdrCVec(jSym,iSym,2)
      call dDaFile(LuCVector(jSym,2),2,CikJ(:,2),nik2*jBas,iAdrJ)
    end if

    do i4=1,iCmp(4)
      lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

      iPAO = iPAO+1
      nijkl = 0

      Factor = Zero
      A(1:jBas*lBas) = Zero
      do iSO=1,nKVec
        nik = nIJ1(iSym,kSym,iSO)

        CiKl(1:nik*lBas) => CijK(iE+1:iE+nik*lBas)

        if (nik == 0) cycle

        if (lSO /= jSO) then
          lSOl = lSO-iOffA
          iAdrL = nik*(lSOl-1)+iAdrCVec(jSym,iSym,iSO)
          call dDaFile(LuCVector(jSym,iSO),2,CiKl,nik*lBas,iAdrL)
          V2(1:) => CiKl(1:)
        else
          V2(1:) => CiKj(1:,iSO)
        end if

        call DGEMM_('T','N',jBas,lBas,nik,One,CikJ(:,iSO),nik,V2,nik,Factor,A,jBas)
        Factor = One
      end do

      do lAOl=0,lBas-1
        lSOl = lSO+lAOl-iOffA
        do jAOj=0,jBas-1
          jSOj = jSO+jAOj-iOffA
          nijkl = nijkl+1

          temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
          temp = temp-Two*ExFac*A(nijkl)

          PMax = max(PMax,abs(temp))
          PAO(nijkl,iPAO) = Fac*temp
        end do
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((iMP2prpt /= 2) .and. lPSO .and. (.not. LSA)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! CASSCF

  iSO = 1

  jSym = 1
  kSym = jSym
  iSym = 1
  lSym = Mul(jSym,iSym)

  nik = nIJ1(iSym,kSym,iSO)
  iS = 1
  iE = nik*jBas
  CiKj(1:nik*jBas,1:1) => CijK(iS:iE)
  iS = iE+1
  iE = iE+nik*lBas
  CiKl(1:nik*lBas) => CijK(iS:iE)

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)

    ! Pick up the MO transformed fitting coefficients, C_ik^J
    jSOj = jSO-iOffA
    iAdrJ = nik*(jSOj-1)+iAdrCVec(jSym,iSym,1)
    call dDaFile(LuCVector(jSym,1),2,CiKj(:,1),nik*jBas,iAdrJ)

    do i4=1,iCmp(4)
      lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

      iPAO = iPAO+1
      nijkl = 0

      if (lSO /= jSO) then
        lSOl = lSO-iOffA
        iAdrL = nik*(lSOl-1)+iAdrCVec(jSym,iSym,1)
        call dDaFile(LuCVector(jSym,1),2,CiKl,nik*lBas,iAdrL)

        V2(1:) => CiKl(1:)
      else
        V2(1:) => CiKj(1:,1)
      end if

      A(1:jBas*lBas) = Zero
      call DGEMM_('T','N',jBas,lBas,nik,One,CiKj(:,1),nik,V2,nik,Zero,A,jBas)

      do lAOl=0,lBas-1
        lSOl = lSO+lAOl-iOffA
        do jAOj=0,jBas-1
          jSOj = jSO+jAOj-iOffA
          nijkl = nijkl+1

          temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)
          temp = temp-ExFac*A(nijkl)

          ! Active space contribution
          temp2 = Zero
          do jp=1,nnP(0)
            temp2 = temp2+sign(One,DMdiag(jp,1))*Z_p_K(jp,jSOj,1)*Z_p_K(jp,lSOl,1)
          end do
          temp = temp+temp2

          PMax = max(PMax,abs(temp))
          PAO(nijkl,iPAO) = Fac*temp
        end do
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((iMP2prpt /= 2) .and. lPSO .and. lSA) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !   SA-CASSCF

  jSym = 1
  kSym = jSym
  iSym = 1
  lSym = Mul(jSym,iSym)

  nik1 = nIJ1(iSym,kSym,1)
  nik2 = nIJ1(iSym,kSym,2)
  nik = max(nik1,nik2)

  iS = 1
  iE = 2*nik*jBas
  CiKj(1:nik*jBas,1:2) => CijK(iS:iE)

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
    jSOj = jSO-iOffA

    ! Pick up the MO transformed fitting coefficients, C_ik^J
    if (nik1 /= 0) then
      iAdrJ = nik1*(jSOj-1)+iAdrCVec(jSym,iSym,1)
      call dDaFile(LuCVector(jSym,1),2,CikJ(:,1),nik1*jBas,iAdrJ)
    end if

    if (nik2 /= 0) then
      iAdrJ = nik2*(jSOj-1)+iAdrCVec(jSym,iSym,2)

      call dDaFile(LuCVector(jSym,2),2,CikJ(:,2),nik2*jBas,iAdrJ)
    end if

    do i4=1,iCmp(4)
      lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)

      iPAO = iPAO+1
      nijkl = 0

      Factor = Zero
      A(1:jBas*lBas) = Zero

      do iSO=1,nKVec
        nik = nIJ1(iSym,kSym,iSO)

        CiKl(1:nik*lBas) => CijK(iE+1:iE+nik*lBas)

        if (nik == 0) cycle

        if (lSO /= jSO) then
          lSOl = lSO-iOffA
          iAdrL = nik*(lSOl-1)+iAdrCVec(jSym,iSym,iSO)
          call dDaFile(LuCVector(jSym,iSO),2,CikL,nik*lBas,iAdrL)
          V2(1:) => CiKl(1:)
        else
          V2(1:) => CiKj(1:,iSO)
        end if

        ! Here one should keep track of negative eigenvalues of the densities

        iSO2 = iSO+2

        do l=1,lBas
          do k=1,jBas

            tmp = Zero

            do i=1,nChOrb(0,iSO)
              do j=1,nChOrb(0,iSO2)

                jik = j+nChOrb(0,iSO2)*(i-1)+nik*(k-1)
                jil = j+nChOrb(0,iSO2)*(i-1)+nik*(l-1)
                if (j <= npos(0,iSO)) then
                  tmp = tmp+CiKj(jik,iSO)*V2(jil)
                else
                  tmp = tmp-CiKj(jik,iSO)*V2(jil)
                end if
              end do
            end do

            kl = k+jBas*(l-1)
            A(kl) = Factor*A(kl)+tmp

          end do
        end do
        Factor = One

      end do

      do lAOl=0,lBas-1
        lSOl = lSO+lAOl-iOffA
        do jAOj=0,jBas-1
          jSOj = jSO+jAOj-iOffA
          nijkl = nijkl+1

          temp = CoulFac*(V_K(lSOl,1)*V_K(jSOj,2)+V_K(lSOl,2)*V_K(jSOj,1)+V_K(lSOl,3)*V_K(jSOj,4)+V_K(lSOl,4)*V_K(jSOj,3)+ &
                          V_K(lSOl,1)*V_K(jSOj,5)+V_K(lSOl,5)*V_K(jSOj,1))
          temp = temp-ExFac*A(nijkl)

          ! Active space contribution
          temp2 = Zero
          do jp=1,nnP(0)
            temp2 = temp2+sign(One,DMdiag(jp,1))*Z_p_K(jp,jSOj,1)*Z_p_K(jp,lSOl,1)+ &
                    sign(Two,DMdiag(jp,2))*(Z_p_K(jp,jSOj,2)*Z_p_K(jp,lSOl,3)+Z_p_K(jp,jSOj,3)*Z_p_K(jp,lSOl,2))
          end do
          temp = temp+temp2

          if (Gamma_On) temp = temp+A_PT2(lSOl,jSOj) ! For CASPT2

          PMax = max(PMax,abs(temp))
          PAO(nijkl,iPAO) = Fac*temp
        end do
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! MP2

  iSO = 1

  jSym = 1
  kSym = jSym
  iSym = 1
  lSym = Mul(jSym,iSym)

  nik = nIJ1(iSym,lSym,iSO)

  iS = 1
  iE = nik*max(jBas,lBas)*2
  CiKj(1:nik*max(jBas,lBas),1:2) => CijK(iS:iE)

  do i2=1,iCmp(2)
    jSO = iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)

    jSOj = jSO-iOffA
    iAdrJ = nik*(jSOj-1)+iAdrCVec(jSym,iSym,1)
    call dDaFile(LuCVector(jSym,1),2,CiKj(:,1),nik*jBas,iAdrJ)

    do i4=1,iCmp(4)
      lSO = iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
      iPAO = iPAO+1
      nijkl = 0

      if (lSO /= jSO) then
        lSOl = lSO-iOffA
        iAdrL = nik*(lSOl-1)+iAdrCVec(jSym,iSym,1)
        call dDaFile(LuCVector(jSym,1),2,CiKj(:,2),nik*lBas,iAdrL)

        V2(1:) => CiKj(1:,2)
      else
        V2(1:) => CiKj(1:,1)
      end if

      A(1:jBas*lBas) = Zero
      call DGEMM_('T','N',jBas,lBas,nik,One,CiKj(:,1),nik,V2,nik,Zero,A,jBas)

      do lAOl=0,lBas-1
        lSOl = lSO+lAOl-iOffA

        ! While the I/O here has been moved outside the
        ! inner loop this needs to be reconsidered and
        ! improved such that it can be moved out yet
        ! another loop (or more.)

        lTot = jBas
        iAdrA = nAuxVe*(lSOl-1)+(jSO-iOffA)
        call dDaFile(LuAVector(1),2,AMP2(:,1),lTot,iAdrA)
        iAdrA = nAuxVe*(lSOl-1)+(jSO-iOffA)
        call dDaFile(LuAVector(2),2,AMP2(:,2),lTot,iAdrA)

        do jAOj=0,jBas-1
          jSOj = jSO+jAOj-iOffA
          nijkl = nijkl+1

          temp = CoulFac*V_K(jSOj,1)*V_K(lSOl,1)+CoulFac*V_K(jSOj,1)*U_K(lSOl)+CoulFac*U_K(jSOj)*V_K(lSOl,1)-ExFac*A(nijkl)

          tempJ_mp2 = AMP2(1+jAOj,2)
          temp = temp+tempJ_mp2*CoulFac

          tempK_mp2 = AMP2(1+jAOj,1)
          temp = temp-ExFac*half*tempK_mp2

          PMax = max(PMax,abs(temp))
          PAO(nijkl,iPAO) = Fac*temp
        end do
      end do
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if

nullify(CiKj,CiKl,V2)
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPAO /= nPAO) then
  write(u6,*) ' Error in PGet1_RI2!'
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet1_RI2:PAO ',' ',PAO,ijkl,nPAO)
#endif
call CWTime(Cpu2,Wall2)
Cpu = Cpu2-Cpu1
Wall = Wall2-Wall1
tavec(1) = tavec(1)+Cpu
tavec(2) = tavec(2)+Wall

return

end subroutine PGet1_RI2
