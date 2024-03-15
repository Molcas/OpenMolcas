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
!***********************************************************************

subroutine PGet2_RI2(iCmp,jBas,lBas,iAO,iAOst,nijkl,PSO,nPSO,ExFac,CoulFac,PMax,V_K,mV_K,Z_p_K,nSA,nZ_p_k)
!***********************************************************************
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density matrix.                 *
!                                                                      *
!          The indices have been scrambled before calling this routine.*
!          Hence we must take special care in order to regain the      *
!          canonical order.                                            *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!                                                                      *
!             Modified to RI-DFT, March 2007                           *
!***********************************************************************

use Basis_Info, only: nBas, nBas_Aux
use SOAO_Info, only: iAOtSO
use pso_stuff, only: DMdiag, lPSO, lSA, nnP
use Symmetry_Info, only: Mul, nIrrep
use RI_glob, only: A, CijK, iAdrCVec, iMP2prpt, iUHF, LuCVector, nIJ1, tavec
use Constants, only: Zero, One, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCmp(4), jBas, lBas, iAO(4), iAOst(4), nijkl, nPSO, mV_K, nSA, nZ_p_k
real(kind=wp), intent(out) :: PSO(nijkl,nPSO), PMax
real(kind=wp), intent(in) :: ExFac, CoulFac, V_K(mV_K,nSA), Z_p_K(nZ_p_k,*)
integer(kind=iwp) :: CumnnP(0:7), CumnnP2(0:7), i, i2, i4, iAdrJ, iAdrL, iE, iS, iSO, iSym, j, j2, j4, jAOj, jp, jpSOj, jpSOl, js, &
                     jSO, jSOj, jSym(0:7), kSym, lAOl, ls, lSO, lSOl, lSym(0:7), MemSO2, mijkl, nB, nik, njSym, nlSym
real(kind=wp) :: Cpu, Cpu1, Cpu2, Fac, Fact, temp, temp2, Wall, Wall1, Wall2
real(kind=wp), pointer :: CiKj(:), CiKl(:), V2(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('V_K',' ',V_K,1,mV_K)
#endif

call CWTime(Cpu1,Wall1)
!                                                                      *
!***********************************************************************
!                                                                      *
PMax = Zero
iSO = 1

PSO(:,:) = Zero

if (lPSO) then
  CumnnP(0) = 0
  CumnnP2(0) = 0
  do i=1,nIrrep-1
    nB = nBas_Aux(i-1)
    if (i == 1) nB = nB-1
    CumnnP(i) = CumnnP(i-1)+nnP(i-1)
    CumnnP2(i) = CumnnP2(i-1)+nnP(i-1)*nB
  end do
end if

!                                                                      *
!***********************************************************************
!                                                                      *
Fac = Quart
MemSO2 = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Pure DFT

if (ExFac == Zero) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do i2=1,iCmp(2)
    njSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) then
        jSym(njSym) = j
        njSym = njSym+1
      end if
    end do

    do i4=1,iCmp(4)
      nlSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(4)+i4,j) > 0) then
          lSym(nlSym) = j
          nlSym = nlSym+1
        end if
      end do
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Loop over irreps which are spanned by the basis function.

      do js=0,njSym-1
        j2 = jSym(js)
        jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)

        do ls=0,nlSym-1
          j4 = lSym(ls)
          if (j2 /= j4) cycle
          lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

          MemSO2 = MemSO2+1
          if (j2 /= 0) cycle

          mijkl = 0
          do lAOl=0,lBas-1
            lSOl = lSO+lAOl-nBas(j4)
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-nBas(j2)
              mijkl = mijkl+1

              ! Coulomb contribution
              if (j2 == 0) then
                ! j4 == 0 also
                temp = V_K(jSOj,1)*V_K(lSOl,1)*CoulFac
                !temp = Zero
              else
                temp = Zero
              end if

              PMax = max(PMax,abs(Temp))
              PSO(mijkl,MemSO2) = Fac*temp

            end do
          end do

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

  do i2=1,iCmp(2)
    njSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) then
        jSym(njSym) = j
        njSym = njSym+1
      end if
    end do

    do i4=1,iCmp(4)
      nlSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(4)+i4,j) > 0) then
          lSym(nlSym) = j
          nlSym = nlSym+1
        end if
      end do
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Loop over irreps which are spanned by the basis function.

      do js=0,njSym-1
        j2 = jSym(js)
        jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)

        do ls=0,nlSym-1
          j4 = lSym(ls)
          if (j2 /= j4) cycle
          lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

          MemSO2 = MemSO2+1

          A(1:jBas*lBas) = Zero

          do iSym=1,nIrrep
            kSym = Mul(j2+1,iSym)
            nik = nIJ1(iSym,kSym,iSO)

            if (nik == 0) cycle

            iS = 1
            iE = nik*jBas
            CiKj(1:nik*jBas) => CijK(iS:iE)

            jSOj = jSO-nBas(j2)
            iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
            call dDaFile(LuCVector(j2+1,iSO),2,CikJ,nik*jBas,iAdrJ)

            if (lSO /= jSO) then
              iS = iE+1
              iE = iE+nik*lBas
              CiKl(1:nik*lBas) => CijK(iS:iE)

              lSOl = lSO-nBas(j4)
              iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
              call dDaFile(LuCVector(j4+1,iSO),2,CiKl,nik*lBas,iAdrL)
              V2(1:) => CiKl(1:)
            else
              V2(1:) => CiKj(1:)
            end if

            Fact = One
            if (iSym /= kSym) Fact = Half
            call DGEMM_('T','N',jBas,lBas,nik,Fact,CikJ,nik,V2,nik,One,A,jBas)

          end do

          mijkl = 0
          do lAOl=0,lBas-1
            lSOl = lSO+lAOl-nBas(j4)
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-nBas(j2)
              mijkl = mijkl+1

              ! Coulomb contribution
              if (j2 == 0) then
                ! j4 == 0 also
                temp = V_K(jSOj,1)*V_K(lSOl,1)*CoulFac
                !temp = Zero
              else
                temp = Zero
              end if

              ! Exchange contribution
              temp = temp-ExFac*A(mijkl)

              PMax = max(PMax,abs(Temp))
              PSO(mijkl,MemSO2) = Fac*temp

            end do
          end do

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

  write(u6,*) 'Pget2_RI2: UDFT/UHF not implemented yet.'
  call Abend()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((iMP2prpt /= 2) .and. lPSO .and. (.not. LSA)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! CASSCF

  do i2=1,iCmp(2)
    njSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) then
        jSym(njSym) = j
        njSym = njSym+1
      end if
    end do

    do i4=1,iCmp(4)
      nlSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(4)+i4,j) > 0) then
          lSym(nlSym) = j
          nlSym = nlSym+1
        end if
      end do
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Loop over irreps which are spanned by the basis function.

      do js=0,njSym-1
        j2 = jSym(js)
        jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)

        do ls=0,nlSym-1
          j4 = lSym(ls)
          if (j2 /= j4) cycle
          lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

          MemSO2 = MemSO2+1

          A(1:jBas*lBas) = Zero

          do iSym=1,nIrrep
            kSym = Mul(j2+1,iSym)
            nik = nIJ1(iSym,kSym,iSO)

            if (nik == 0) cycle

            iS = 1
            iE = nik*jBas
            CiKj(1:nik*jBas) => CijK(iS:iE)

            jSOj = jSO-nBas(j2)
            iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
            call dDaFile(LuCVector(j2+1,iSO),2,CiKj,nik*jBas,iAdrJ)

            if (lSO /= jSO) then
              iS = iE+1
              iE = iE+nik*lBas
              CiKl(1:nik*lBas) => CijK(iS:iE)

              lSOl = lSO-nBas(j4)
              iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
              call dDaFile(LuCVector(j4+1,iSO),2,CiKl,nik*lBas,iAdrL)
              V2(1:) => CiKl(1:)
            else
              V2(1:) => CiKj(1:)
            end if

            Fact = One
            if (iSym /= kSym) Fact = Half
            call DGEMM_('T','N',jBas,lBas,nik,Fact,CikJ,nik,V2,nik,One,A,jBas)

          end do

          mijkl = 0
          do lAOl=0,lBas-1
            lSOl = lSO+lAOl-nBas(j4)
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-nBas(j2)
              mijkl = mijkl+1

              ! Coulomb contribution
              if (j2 == 0) then
                ! j4 == 0 also
                temp = V_K(jSOj,1)*V_K(lSOl,1)*CoulFac
                !temp = Zero
              else
                temp = Zero
              end if

              ! Exchange contribution
              temp = temp-ExFac*A(mijkl)

              temp2 = Zero
              jpSOj = CumnnP2(j2)+(jSOj-1)*nnP(j2)
              jpSOl = CumnnP2(j2)+(lSOl-1)*nnP(j2)
              do jp=1,nnP(j2)
                temp2 = temp2+sign(One,DMdiag(CumnnP(j2)+jp,1))*Z_p_K(jpSOj+jp,1)*Z_p_K(jpSOl+jp,1)
              end do
              temp = temp+temp2

              PMax = max(PMax,abs(Temp))
              PSO(mijkl,MemSO2) = Fac*temp

            end do
          end do

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
  ! SA-CASSCF

  write(u6,*) 'Pget2_ri2: SA-CASSCF not implemented yet'
  call Abend()

  do i2=1,iCmp(2)
    njSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) then
        jSym(njSym) = j
        njSym = njSym+1
      end if
    end do

    do i4=1,iCmp(4)
      nlSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(4)+i4,j) > 0) then
          lSym(nlSym) = j
          nlSym = nlSym+1
        end if
      end do
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Loop over irreps which are spanned by the basis function.

      do js=0,njSym-1
        j2 = jSym(js)
        jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)

        do ls=0,nlSym-1
          j4 = lSym(ls)
          if (j2 /= j4) cycle
          lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

          MemSO2 = MemSO2+1

          A(1:jBas*lBas) = Zero

          do iSym=1,nIrrep
            kSym = Mul(j2+1,iSym)
            nik = nIJ1(iSym,kSym,iSO)

            if (nik == 0) cycle

            iS = 1
            iE = nik*jBas
            CiKj(1:nik*jBas) => CijK(iS:iE)

            jSOj = jSO-nBas(j2)
            iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
            call dDaFile(LuCVector(j2+1,iSO),2,CiKj,nik*jBas,iAdrJ)

            if (lSO /= jSO) then
              iS = iE+1
              iE = iE+nik*lBas
              CiKl(1:nik*lBas) => CijK(iS:iE)

              lSOl = lSO-nBas(j4)
              iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
              call dDaFile(LuCVector(j4+1,iSO),2,CiKl,nik*lBas,iAdrL)
              V2(1:) => CiKl(1:)
            else
              V2(1:) => CiKj(1:)
            end if

            Fact = One
            if (iSym /= kSym) Fact = Half
            call DGEMM_('T','N',jBas,lBas,nik,Fact,CiKJ,nik,V2,nik,One,A,jBas)

          end do

          mijkl = 0
          do lAOl=0,lBas-1
            lSOl = lSO+lAOl-nBas(j4)
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-nBas(j2)
              mijkl = mijkl+1

              ! Coulomb contribution
              if (j2 == 0) then
                ! j4 == 0 also
                temp = CoulFac*(V_K(lSOl,1)*V_K(jSOj,2)+V_K(lSOl,2)*V_K(jSOj,1)+V_K(lSOl,3)*V_K(jSOj,4)+V_K(lSOl,4)*V_K(jSOj,3))
                !temp = Zero
              else
                temp = Zero
              end if

              ! Exchange contribution
              temp = temp-ExFac*A(mijkl)

              PMax = max(PMax,abs(Temp))
              PSO(mijkl,MemSO2) = Fac*temp

            end do
          end do

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

  write(u6,*) 'Pget2_ri2: MP2 not implemented yet'
  call Abend()

  do i2=1,iCmp(2)
    njSym = 0
    do j=0,nIrrep-1
      if (iAOtSO(iAO(2)+i2,j) > 0) then
        jSym(njSym) = j
        njSym = njSym+1
      end if
    end do

    do i4=1,iCmp(4)
      nlSym = 0
      do j=0,nIrrep-1
        if (iAOtSO(iAO(4)+i4,j) > 0) then
          lSym(nlSym) = j
          nlSym = nlSym+1
        end if
      end do
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Loop over irreps which are spanned by the basis function.

      do js=0,njSym-1
        j2 = jSym(js)
        jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)

        do ls=0,nlSym-1
          j4 = lSym(ls)
          if (j2 /= j4) cycle
          lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)

          MemSO2 = MemSO2+1

          A(1:jBas*lBas) = Zero

          do iSym=1,nIrrep
            kSym = Mul(j2+1,iSym)
            nik = nIJ1(iSym,kSym,iSO)

            if (nik == 0) cycle

            iS = 1
            iE = nik*jBas
            CiKj(1:nik*jBas) => CijK(iS:iE)

            jSOj = jSO-nBas(j2)
            iAdrJ = nik*(jSOj-1)+iAdrCVec(j2+1,iSym,iSO)
            call dDaFile(LuCVector(j2+1,iSO),2,CiKj,nik*jBas,iAdrJ)

            if (lSO /= jSO) then
              iS = iE+1
              iE = iE+nik*lBas
              CiKl(1:nik*lBas) => CijK(iS:iE)

              lSOl = lSO-nBas(j4)
              iAdrL = nik*(lSOl-1)+iAdrCVec(j4+1,iSym,iSO)
              call dDaFile(LuCVector(j4+1,iSO),2,CiKl,nik*lBas,iAdrL)
              V2(1:) => CiKl(1:)
            else
              V2(1:) => CiKj(1:)
            end if

            Fact = One
            if (iSym /= kSym) Fact = Half
            call DGEMM_('T','N',jBas,lBas,nik,Fact,CiKj,nik,V2,nik,One,A,jBas)

          end do

          mijkl = 0
          do lAOl=0,lBas-1
            lSOl = lSO+lAOl-nBas(j4)
            do jAOj=0,jBas-1
              jSOj = jSO+jAOj-nBas(j2)
              mijkl = mijkl+1

              ! Coulomb contribution
              if (j2 == 0) then
                ! j4 == 0 also
                temp = V_K(jSOj,1)*V_K(lSOl,1)*CoulFac
                !temp = Zero
              else
                temp = Zero
              end if

              ! Exchange contribution
              temp = temp-ExFac*A(mijkl)

              PMax = max(PMax,abs(Temp))
              PSO(mijkl,MemSO2) = Fac*temp

            end do
          end do

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
if (nPSO /= MemSO2) then
  write(u6,*) ' PGet2: nPSO /= MemSO2'
  write(u6,*) nPSO,MemSO2
  call Abend()
end if

#ifdef _DEBUGPRINT_
call RecPrt(' In PGet2_RI2:PSO ',' ',PSO,nijkl,nPSO)
#endif

call CWTime(Cpu2,Wall2)
Cpu = Cpu2-Cpu1
Wall = Wall2-Wall1
tavec(1) = tavec(1)+Cpu
tavec(2) = tavec(2)+Wall

return

end subroutine PGet2_RI2
