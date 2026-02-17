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

subroutine SXHAM(D,P,PA,FP,SXN,F1,F2,DIA,G,H,HDIAG,DF,DDIAG)
! RASSCF VERSION IBM-3090: SX-SECTION
!
! Objective: Compute the auxiliary matrices G,H,F1,F2,and SXN
! defined as:
!
! G(pq)=sum(r,s)(2*P(pqrs)-D(pq)*D(rs))*FP(rs)
! H(pq)=G(pq)+DF(pq)+DF(qp)  where
! DF(pq)=sum(r)D(pr)*FP(rq)
! F1 is the occupied part of the SX Fock matrix FP
! F2 is the active-external part of the same matrix
! DIA is the occupied part of the density matrix(squared)
! SXN contains normalization factors for the SX states
! HDIAG is the diagonal of the SX Hamiltonian
! All matrices are blocked according to symmetry
! Each block having the dimension
! G, DIA, AND F1:     (NIO+NAO)**2
! F2:                 (NAO+NEO)**2
! H:                  NAO*NAE
! SXN:                (NIO+NAO)*(NAO+NEO)
! HDIAG:              NROOT+(NIO+NAO)*(NAO+NEO)
!
! These matrices are used to simplify the calculation
! of the sigma vector (see SIGVEC for further details).
!
! Called from SXCTL
!
! Subroutine calls: none
!
! ********** IBM-3090 RELEASE 89 01 23 **********

use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: Ener, ExFac, ICICP, IROOT, ISCF, ITER, IXSYM, IZROT, LVSHFT, NROOT, NSXS, SXSHFT
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: NASH, NBAS, NFRO, NISH, NORB, NSSH, NSYM
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: D(*), P(*), PA(*), FP(*), SXN(*), F1(*), F2(*), DIA(*), G(*), H(*), HDIAG(*), DF(*), DDIAG(*)
integer(kind=iwp) :: I, IASHI, IASHJ, IPQ, IPRLEV, IQP, IROOT1, ISTAE, ISTBM, ISTD, ISTFP, ISTFPJ, ISTH, ISTIA, ISTZ, ISYM, IX, &
                     IX1, JSYM, JUSTONE, NAE, NAO, NAOJ, NEO, NIA, NIO, NIOJ, NO, NP, nP2Act, NPR, NQ, NR, NRR, NT, NTT, NTU, &
                     NTUT, NTUTU, NTUVX, NU, NUT, NV, NVT, NVX, NVXF, NVXT, NX, NXT
real(kind=wp) :: DPP, DRR, DTU, FAC, FACD, FPP, GTU, HDMIN, HPP, P2Act(1), PRPR, SXNRM2, XLEV
real(kind=wp), parameter :: THRA = 1.0e-6_wp
#include "warnings.h"

! -- THRA: THRESHOLD FOR WARNING, ACTIVE OCC NO CLOSE TO 0 OR 2.

! Local print level (if any)
IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering SXHAM'

! Loop over all symmetry blocks

ISTBM = 0
ISTD = 0
ISTFP = 0
ISTIA = 0
ISTAE = 0
ISTH = 0
IASHI = 0
ISTZ = 0
IX1 = 0
do ISYM=1,NSYM
  IX = IX1+NFRO(ISYM)
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  !NR1O = NRS1(ISYM)
  !NR2O = NRS2(ISYM)
  !NR3O = NRS3(ISYM)
  !NAO = NR1O+NR2O+NR3O
  NIA = NIO+NAO
  NEO = NSSH(ISYM)
  NAE = NAO+NEO
  NO = NORB(ISYM)
  if ((NIA == 0) .or. (NAE == 0)) GO TO 98

  ! Form the diagonal of the density matrix for all orbitals

  justone = 0
  do NP=1,NO
    if (NP <= NIO) DDIAG(NP) = Two
    if (NP > NIA) DDIAG(NP) = Zero
    if ((NP > NIO) .and. (NP <= NIA)) then
      DDIAG(NP) = D(ISTD+nTri_Elem(NP-NIO))
      if ((DDIAG(NP) < THRA) .and. (ISCF == 0)) then
        if (JustOne == 0) then
          call WarningMessage(1,'Problems in orbital optimization.')
          JustOne = JustOne+1
        end if
        write(u6,'(1x,a,i2,a,i4,a,es14.6)') ' Warning: In symmetry ',ISYM,', orbital p=',NP, &
                                            ' has diagonal density matrix element D(p,p) close to zero. D(p,p)=',DDIAG(NP)
      end if
      if ((Two-DDIAG(NP) < THRA) .and. (ISCF == 0)) then
        if (JustOne == 0) then
          call WarningMessage(1,'Problems in orbital optimization.')
          JustOne = JustOne+1
        end if
        write(u6,'(1x,a,i2,a,i4,a,es14.6)') ' Warning: In symmetry ',ISYM,', orbital p=',NP, &
                                            ' has diagonal density matrix element D(p,p) close to two. (2 - D(p,p))=',Two-DDIAG(NP)
      end if
    end if
  end do

  ! Compute the normalization constants:
  ! SXN(pq)=4*(P(prrp)-P(prpr))+D(rr)+D(pp)

  NPR = ISTBM
  ! NP loops over active and secondary indices.
  do NP=NIO+1,NO
    ! NR loops over inactive and active indices.
    do NR=1,NIA
      NPR = NPR+1
      SXN(NPR) = One
      if (NR >= NP) GO TO 95
      DRR = DDIAG(NR)
      DPP = DDIAG(NP)
      PRPR = -Two*DPP
      if (NP > NIA) PRPR = Zero
      if ((NR > NIO) .and. (NP <= NIA)) then
        NT = NP-NIO+IASHI
        NU = NR-NIO+IASHI
        NTU = iTri(NT,NU)
        NTUTU = nTri_Elem(NTU)
        PRPR = -Four*PA(NTUTU)
      end if
      SXNRM2 = PRPR+DRR+DPP
      if (SXNRM2 < -1.0e-12_wp) then
        call WarningMessage(1,'Negative norm occured in SXHAM.')
        write(u6,*) 'SXHAM Error: Negative SXNRM2.'
        write(u6,*) ' Symmetry block ISYM:',ISYM
        write(u6,*) '      Orbitals NP,NR:',NP,NR
        write(u6,*) '                PRPR:',PRPR
        write(u6,*) '                 DRR:',DRR
        write(u6,*) '                 DPP:',DPP
        write(u6,*) '      Norm**2 SXNRM2:',SXNRM2
        write(u6,*)
        write(u6,*) ' The squared norm of a trial vector has been'
        write(u6,*) ' computed to be negative.'
        write(u6,*) ' This is possible only for some severe malfunction'
        write(u6,*) ' of the rasscf program. Please issue a bug report.'
        write(u6,*)
        call Quit(_RC_GENERAL_ERROR_)
      end if

      if (SXNRM2 < 1.0e-12_wp) then
        SXN(NPR) = Zero
      else
        SXN(NPR) = One/sqrt(SXNRM2)
      end if
95    continue
    end do
  end do

  ! Form the occupied (F1) and active-external (F2) part of FP
  ! and the occupied part DIA of the density matrix D

  IPQ = ISTFP
  DIA(ISTIA+1:ISTIA+NIA**2) = Zero
  do NP=1,NO
    do NQ=1,NP
      IPQ = IPQ+1
      if ((NP <= NIA) .and. (NQ <= NIA)) then
        F1(ISTIA+NIA*(NP-1)+NQ) = FP(IPQ)
        F1(ISTIA+NIA*(NQ-1)+NP) = FP(IPQ)
      end if
      if ((NP > NIO) .and. (NQ > NIO)) then
        F2(NAE*(NP-NIO-1)+NQ-NIO+ISTAE) = FP(IPQ)
        F2(NAE*(NQ-NIO-1)+NP-NIO+ISTAE) = FP(IPQ)
      end if
      if ((NP <= NIO) .and. (NP == NQ)) DIA(ISTIA+NIA*(NP-1)+NP) = Two
      if (((NP > NIO) .and. (NP <= NIA)) .and. NQ > NIO) then
        NT = NP-NIO
        NU = NQ-NIO
        NTU = iTri(NT,NU)+ISTD
        DIA(ISTIA+NIA*(NP-1)+NQ) = D(NTU)
        DIA(ISTIA+NIA*(NQ-1)+NP) = D(NTU)
      end if
    end do
  end do
  if (IPRLEV >= DEBUG) then
    write(u6,*) 'DIA in SXHAM:'
    write(u6,*) (DIA(ISTIA+i),i=1,NIA**2)
  end if

  ! Form the matrix DF (row index active, column index all orbitals)

  if (NAO /= 0) then
    call DGEMM_('N','N',NAO,NAE,NAO,One,DIA(ISTIA+NIO*NIA+NIO+1),NIA,F2(ISTAE+1),NAE,Zero,DF(NAO*NIO+1),NAO)
    if (NIO /= 0) call DGEMM_('N','N',NAO,NIO,NAO,One,DIA(ISTIA+NIO*NIA+NIO+1),NIA,F1(ISTIA+NIO+1),NIA,Zero,DF,NAO)
  end if

  ! compute the G matrix:
  ! G(ij)=-2*FP(ij)
  ! G(ti)=G(it)=-DF(ti)
  ! G(tu)=sum(vx)(2P(tuvx)-D(tu)D(vx))F(vx)

  if (ExFac /= One) then
    call Get_Temp('nP2Act  ',P2Act,1)
    nP2Act = int(P2Act(1))
    call Get_Temp('P2_RAW  ',P,nP2Act)
  end if
  do NP=1,NIA
    do NQ=1,NP
      IPQ = ISTIA+NIA*(NP-1)+NQ
      IQP = ISTIA+NIA*(NQ-1)+NP
      if (NP <= NIO) then
        G(IPQ) = -Two*F1(IPQ)
        G(IQP) = G(IPQ)
      else if ((NP > NIO) .and. (NQ <= NIO)) then
        G(IPQ) = -DF(NAO*(NQ-1)+NP-NIO)
        G(IQP) = G(IPQ)
      else if (NQ > NIO) then
        DTU = DIA(ISTIA+NIA*(NP-1)+NQ)
        GTU = Zero
        NTT = NP-NIO+IASHI
        NUT = NQ-NIO+IASHI
        NTUT = iTri(NTT,NUT)

        IASHJ = 0
        ISTFPJ = 0
        NVX = 0
        do JSYM=1,NSYM
          NAOJ = NASH(JSYM)
          NIOJ = NISH(JSYM)
          if (NAOJ == 0) GO TO 23
          do NV=1,NAOJ
            NVT = NV+IASHJ
            do NX=1,NV
              NXT = NX+IASHJ
              NVX = NVX+1
              NVXT = iTri(NVT,NXT)
              NTUVX = iTri(NTUT,NVXT)
              FAC = Two
              FACD = Two
              if (NV == NX) FACD = One
              if (NTUT < NVXT) then
                if ((NVT /= NXT) .and. (NTT == NUT)) FAC = Four
                if ((NVT == NXT) .and. (NTT /= NUT)) FAC = One
              end if
              NVXF = ISTFPJ+iTri(NV+NIOJ,NX+NIOJ)
              GTU = GTU+FP(NVXF)*(FAC*P(NTUVX)-FACD*DTU*D(NVX))
            end do
          end do
          IASHJ = IASHJ+NAOJ
23        ISTFPJ = ISTFPJ+nTri_Elem(NORB(JSYM))
        end do
        G(IPQ) = GTU
        G(IQP) = GTU
      end if
    end do
  end do
  if (ExFac /= One) call Get_Temp('P2_KS   ',P,nP2Act)

  ! FORM THE H MATRIX (only the nae*nao block is needed)
  ! H(tu)=G(tu)+DF(tu)+DF(ut)
  ! H(at)=DF(ta)
  ! H(ab) and H(ta) are not included.

  if (NAO /= 0) then
    IPQ = ISTH
    do NP=1,NAO
      do NQ=1,NAE
        IPQ = IPQ+1
        if (NQ <= NAO) then
          H(IPQ) = G(ISTIA+NIA*(NP+NIO-1)+NQ+NIO)+DF(NAO*(NP+NIO-1)+NQ)+DF(NAO*(NQ+NIO-1)+NP)
        else
          H(IPQ) = DF(NAO*(NQ+NIO-1)+NP)
        end if
      end do
    end do
  end if

  ! Diagonal elements of the SX Hamiltonian

  NPR = ISTBM
  do NP=1,NAE
    FPP = F2(ISTAE+NAE*(NP-1)+NP)
    HPP = Zero
    DPP = Zero
    if (NP <= NAO) then
      HPP = H(ISTH+NAE*(NP-1)+NP)
      DPP = DIA(ISTIA+NIA*(NP+NIO-1)+NP+NIO)
    end if
    NRR = ISTIA+1
    do NR=1,NIA
      NPR = NPR+1
      HDIAG(NPR+NROOT) = (G(NRR)-HPP+DIA(NRR)*FPP+DPP*F1(NRR))*SXN(NPR)**2
      NRR = NRR+NIA+1

      ! Make this matrix element large for forbidden rotations

      if ((NP <= NAO) .and. (NR > NIO)) then
        NT = NP
        NU = NR-NIO
        if (NT <= NU) then
          HDIAG(NPR+NROOT) = 1.0e32_wp
          SXN(NPR) = Zero
        else
          NTU = ISTZ+nTri_Elem(NT-2)+NU
          if (IZROT(NTU) /= 0) then
            HDIAG(NPR+NROOT) = 1.0e32_wp
            SXN(NPR) = Zero
            !write(u6,*) 'SXHAM. IZROT == 1 for NP,NR=',NP,NR
          end if
        end if
      end if

      ! Make HDIAG large for rotations forbidden by IXSYM input

      if (IXSYM(NR+IX) /= IXSYM(NP+NIO+IX)) then
        !write(u6,*) 'SXHAM. IXSYM forbids NP,NR=',NP,NR
        HDIAG(NPR+NROOT) = 1.0e32_wp
        SXN(NPR) = Zero
      end if

    end do
  end do
  !write(u6,*) 'SXHAM: The IXSYM array='
  !write(u6,'(1x,40i2)') (ixsym(i),i=1,ntot)

  ! Test print of all matrices

  if (IPRLEV >= DEBUG) then
    write(u6,1000) ISYM
    write(u6,1100) (DDIAG(I),I=1,NO)
    write(u6,1200) (SXN(I),I=ISTBM+1,ISTBM+NIA*NAE)
    write(u6,1300) (DF(I),I=1,NAO*NO)
    write(u6,1400) (F1(I),I=ISTIA+1,ISTIA+NIA**2)
    write(u6,1500) (F2(I),I=ISTAE+1,ISTAE+NAE**2)
    write(u6,1600) (G(I),I=ISTIA+1,ISTIA+NIA**2)
    write(u6,1700) (H(I),I=ISTH+1,ISTH+NAO*NAE)
    write(u6,1800) (HDIAG(I),I=ISTBM+1+NROOT,ISTBM+NAE*NIA+NROOT)
  end if

  IASHI = IASHI+NAO
  ISTD = ISTD+nTri_Elem(NAO)
98 ISTIA = ISTIA+NIA**2
  ISTAE = ISTAE+NAE**2
  ISTFP = ISTFP+nTri_Elem(NO)
  ISTBM = ISTBM+NIA*NAE
  ISTH = ISTH+NAO*NAE
  ISTZ = ISTZ+nTri_Elem(NAO-1)
  IX1 = IX1+NBAS(ISYM)

end do

! Level shift section

XLEV = LVSHFT
SXSHFT = Zero
HDMIN = 1.0e32_wp
do I=NROOT+1,NROOT+NSXS
  if (HDIAG(I) < HDMIN) HDMIN = HDIAG(I)
end do
SXSHFT = max(XLEV-HDMIN,Zero)

! no level shift if input (alpha) is zero

if (LVSHFT == Zero) SXSHFT = Zero
do I=NROOT+1,NROOT+NSXS
  HDIAG(I) = HDIAG(I)+SXSHFT
end do
if (IPRLEV >= DEBUG) write(u6,1900) HDMIN,SXSHFT

! Add diagonal elements for the reference space (CI-states)

HDIAG(1) = Zero
if (ICICP /= 0) then
  IROOT1 = IROOT(1)
  do I=1,NROOT
    HDIAG(I) = ENER(I,ITER)-ENER(IROOT1,ITER)
  end do
end if

return

1000 format(/1X,'Matrices in SXHAM for symmetry block',I2)
1100 format(/1X,'Diagonal of the full density matrix'/(1X,5G18.6))
1200 format(/1X,'Normalization constants for SX-states'/(1X,5G18.6))
1300 format(/1X,'The matrix sum(D(pr)F(rq)) -DF'/(1X,5G18.6))
1400 format(/1X,'The occupied part of the Fock matrix'/(1X,5G18.6))
1500 format(/1X,'The act-ext part of the Fock matrix'/(1X,5G18.6))
1600 format(/1X,'The G-matrix'/(1X,5G18.6))
1700 format(/1X,'The H-matrix'/(1X,5G18.6))
1800 format(/1X,'Diagonal of the SX Hamiltonian'/(1X,5G18.6))
1900 format(1X,'Lowest diagonal element has the value',F12.6,' A level shift of',F12.6,' has been used')

end subroutine SXHAM
