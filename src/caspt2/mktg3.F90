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

subroutine MKTG3(LSYM1,LSYM2,CI1,CI2,OVL,TG1,TG2,NTG3,TG3)
! Procedure for computing 1-body, 2-body, and 3-body transition
! density elements with active indices only.

! In: Wave functions CI1, with symmetry LSYM1, and CI2, with
!  symmetry LSYM2.
!
! Out: Transition density matrices, denoted here TG1, TG2 and TG3.
! Storage: TG1 and TG2 are simple two- and four-index arrays, and
! includes also such zeroes that are implied by symmetry.
! But TG3 is quite large, and while it is stored with zeroes, it
! is made more compact by the following addressing:

! <Psi1|E_tuvxyz|Psi2> is stored in TG3(ITG3) where
!    ITG3= nTri3_Elem(i-1) + nTri_Elem(j-1) + k
!     i  = max(tu,vx,yz)
!     j  = mid(tu,vx,yz)
!     k  = min(tu,vx,yz)
! tu stands for the pair index tu= t + NASHT*(u-1), etc., and t is
! the usual active orbital number, when they are enumerated across
! all the symmetries (The "absolute" active index).

use Index_Functions, only: nTri_Elem, nTri3_Elem
use Symmetry_Info, only: Mul
use sguga, only: CIS, EXS, L2ACT, SGS
use caspt2_module, only: IASYM, ISCF, NACTEL, NASHT
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs, MyRank
#endif
use caspt2_module, only: MxCI
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LSYM1, LSYM2, NTG3
real(kind=wp), intent(in) :: CI1(MXCI), CI2(MXCI)
real(kind=wp), intent(out) :: OVL, TG1(NASHT,NASHT), TG2(NASHT,NASHT,NASHT,NASHT), TG3(NTG3)
integer(kind=iwp) :: IL, IND1, IND2, IND3, IP, IP1, IP1END, IP1STA, IP2, IP3, IP3END, IP3STA, IS1, IS2, IS3, ISSG1, ISSG2, ISTAU, &
                     IT, IT1, IT2, IT3, ITG3, ITS, IU, IU1, IU2, IU3, IUS, IV, IVS, IX, IXS, IY, IYS, IZ, IZS, JL, jtuvxyz, L, &
                     LFROM, LSGM1, LSGM2, LTAU, LTO, NCI1, nLev, NTAU, NTG3WRK, NTUBUF, NVECS, NYZBUF
real(kind=wp) :: OCC, VAL
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iTask
logical(kind=iwp) :: Poor_Par
#endif
integer(kind=iwp), allocatable :: P2LEV(:,:)
real(kind=wp), allocatable :: TG3WRK(:)
real(kind=wp), external :: DDot_

nLev = SGS%nLev

! Put in zeroes. Recognize special cases:
OVL = One
if (NASHT == 0) return
if (LSYM1 /= LSYM2) OVL = Zero
TG1(:,:) = Zero
TG2(:,:,:,:) = Zero
TG3(:) = Zero

if (NACTEL == 0) return

if (ISCF == 0) then

  ! Here, for regular CAS or RAS cases.

  ! Special pair index allows true RAS cases to be handled:
  call mma_allocate(P2LEV,2,NASHT**2,Label='P2LEV')
  IP = 0
  ! First, IL < JL pairs.
  do IL=1,NLEV-1
    P2LEV(1,IP+1:IP+NLEV-IL) = IL
    P2LEV(2,IP+1:IP+NLEV-IL) = [(JL,JL=IL+1,NLEV)]
    IP = IP+NLEV-IL
  end do
  ! Then, IL = JL pairs.
  P2LEV(1,IP+1:IP+NLEV) = [(IL,IL=1,NLEV)]
  P2LEV(2,IP+1:IP+NLEV) = [(IL,IL=1,NLEV)]
  IP = IP+NLEV
  ! Last, IL > JL pairs.
  do IL=2,NLEV
    P2LEV(1,IP+1:IP+IL-1) = IL
    P2LEV(2,IP+1:IP+IL-1) = [(JL,JL=1,IL-1)]
    IP = IP+IL-1
  end do
  ! If now any matrix element E(t1u1)E(t2u2)..E(tnun) is arranged
  ! such that the pair indices are non-decreasing, then the matrix
  ! element can be correctly computed by performing explicit
  ! excitations within the RAS space.
  ! But we also need the 'usual' pair index in order to use the
  ! packed addressing.

  NCI1 = CIS%NCSF(LSYM1)
  ! Overlap:
  if (LSYM1 == LSYM2) OVL = DDOT_(NCI1,CI1,1,CI2,1)
  ! Allocate as many vectors as possible:
  ! Wishful thinking:
  NVECS = 2*NASHT**2+1
  ! But what is really available?
  call mma_MaxDBLE(NTG3WRK)
  NTG3WRK = min(MXCI*NVECS,NTG3WRK)
  NVECS = NTG3WRK/MXCI
  NTG3WRK = NVECS*MXCI
  ! Find optimal subdivision of available vectors:
  NYZBUF = nint(real(NVECS-1,kind=wp)/real(NASHT,kind=wp))
  NYZBUF = max(1,NYZBUF)
  NTUBUF = min(NASHT**2,NVECS-1-NYZBUF)
  NYZBUF = NVECS-1-NTUBUF
  ! Insufficient memory?
  if (NTUBUF <= 0) then
    write(u6,*) ' Too little memory left for MKTG3.'
    write(u6,*) ' Need at least 3 vectors of length MXCI=',MXCI
    call ABEND()
  end if
  if (NTUBUF <= (NASHT**2)/5) then
    write(u6,*) ' WARNING: MKTG3 will be inefficient owing to'
    write(u6,*) ' small memory.'
  end if
  call mma_allocate(TG3WRK,NTG3WRK,Label='TG#WRK')
  ! And divide it up:
  LSGM1 = 1
  LTAU = LSGM1+NTUBUF*MXCI
  LSGM2 = LTAU+MXCI

# ifdef _MOLCAS_MPP_
  !! enable poor parallelization, if applicable
  if (Is_Real_Par()) then
    POOR_PAR = .false.
    iTask = 0
    !if ((NTUBUF == NYZBUF) .and. (NTUBUF == NASHT**2)) POOR_PAR = .true.
  end if
# endif
  ! Sectioning loops over pair indices IP3 (ket side):
  do IP3STA=1,NASHT**2,NYZBUF
    IP3END = min(NASHT**2,IP3STA-1+NYZBUF)
    ! Compute a section of sigma vectors E(YZ)*PSI2 to memory:
    LTO = LSGM2
    do IP3=IP3STA,IP3END
      ! Translate to levels in the SGUGA coupling order:
      IL = P2LEV(1,IP3)
      JL = P2LEV(2,IP3)
      IY = L2ACT(IL)
      IZ = L2ACT(JL)
      IYS = IASYM(IY)
      IZS = IASYM(IZ)
      ISSG2 = Mul(Mul(IYS,IZS),LSYM2)
      TG3WRK(LTO:LTO+MXCI-1) = Zero
      ! LTO is first element of Sigma2 = E(YZ) Psi2
      call SG_Epq_Psi(SGS,CIS,EXS,IL,JL,One,LSYM2,CI2,TG3WRK(LTO))
      if (ISSG2 == LSYM1) TG1(IY,IZ) = DDOT_(NCI1,CI1,1,TG3WRK(LTO),1)
      LTO = LTO+MXCI
    end do
    ! Sectioning loops over pair indices IP1 (bra side):
    do IP1STA=IP3STA,NASHT**2,NTUBUF
      IP1END = min(NASHT**2,IP1STA-1+NTUBUF)
      ! Compute a section of sigma vectors E(UT)*PSI1 to memory:
      LTO = LSGM1
      do IP1=IP1STA,IP1END
        ! Translate to levels:
        JL = P2LEV(1,IP1)
        IL = P2LEV(2,IP1)
        IT = L2ACT(IL)
        IU = L2ACT(JL)
        ITS = IASYM(IT)
        IUS = IASYM(IU)
        ISSG1 = Mul(Mul(ITS,IUS),LSYM1)
        TG3WRK(LTO:LTO+MXCI-1) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,IL,JL,One,LSYM1,CI1,TG3WRK(LTO))
        LTO = LTO+MXCI
      end do
      ! Now compute as many elements as possible:
      LFROM = LSGM2
      do IP3=IP3STA,IP3END
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          iTask = iTask+1
          if (POOR_PAR .and. (mod(iTask,nProcs) /= MyRank)) then
            LFROM = LFROM+MXCI
            cycle
          end if
        end if
#       endif
        IY = L2ACT(P2LEV(1,IP3))
        IZ = L2ACT(P2LEV(2,IP3))
        ! LFROM will be start element of Sigma2=E(YZ) Psi2
        IYS = IASYM(IY)
        IZS = IASYM(IZ)
        ISSG2 = Mul(Mul(IYS,IZS),LSYM2)
        do IP2=IP3,IP1END
          IL = P2LEV(1,IP2)
          JL = P2LEV(2,IP2)
          IV = L2ACT(IL)
          IX = L2ACT(JL)
          IVS = IASYM(IV)
          IXS = IASYM(IX)
          ISTAU = Mul(Mul(IVS,IXS),ISSG2)
          NTAU = CIS%NCSF(ISTAU)
          TG3WRK(LTAU:LTAU+MXCI-1) = Zero
          ! LTAU  will be start element of Tau=E(VX) Sigma2=E(VX) E(YZ) Psi2
          call SG_Epq_Psi(SGS,CIS,EXS,IL,JL,One,ISSG2,TG3WRK(LFROM),TG3WRK(LTAU))
          if (ISTAU == LSYM1) TG2(IV,IX,IY,IZ) = DDOT_(NTAU,TG3WRK(LTAU),1,CI1,1)
          do IP1=max(IP2,IP1STA),IP1END
            IT = L2ACT(P2LEV(1,IP1))
            IU = L2ACT(P2LEV(2,IP1))
            ITS = IASYM(IT)
            IUS = IASYM(IU)
            ISSG1 = Mul(Mul(ITS,IUS),LSYM1)
            if (ISSG1 == ISTAU) then
              L = LSGM1+MXCI*(IP1-IP1STA)
              VAL = DDOT_(NTAU,TG3WRK(LTAU),1,TG3WRK(L),1)
              ! Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
              ! Code to put it in correct place:
              call get_tg3_index(IT,IU,IV,IX,IY,IZ,NASHT,jtuvxyz)
              TG3(JTUVXYZ) = VAL

              ! End of symmetry requirement IF-clause:
            end if
            ! End of IP1 loop.
          end do
          ! End of IP2 loop.
        end do
        LFROM = LFROM+MXCI
        ! End of IP3 loop.
      end do
      ! End of IP1STA sectioning loop
    end do
    ! End of IP3STA sectioning loop
  end do
  call mma_deallocate(TG3WRK)
  ! Now the computed elements of TG2 contain <PSI1|E(IT1,IU1)E(IT2,IU2)|PSI2>
  ! and TG3 contains <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
  ! Add here the necessary Kronecker deltas times 2-body matrix
  ! elements and lower, so we get a true normal-ordered density matrix
  ! element.

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par() .and. POOR_PAR) then
    call GADGOP(TG2,NASHT**4,'+')
    call GADGOP(TG3,NTG3,'+')
  end if
# endif

  ! First, the 2-particle density matrix:
  ! <PSI1|E(T,U,V,X)|PSI2>  = <PSI1|E(TU)E(VX)|PSI2> - D(V,U)*TG2(T,U,V,X)
  do IP1=1,NASHT**2
    IT = L2ACT(P2LEV(1,IP1))
    IU = L2ACT(P2LEV(2,IP1))
    do IP2=1,IP1
      IV = L2ACT(P2LEV(1,IP2))
      IX = L2ACT(P2LEV(2,IP2))
      if (IV == IU) TG2(IT,IU,IV,IX) = TG2(IT,IU,IV,IX)-TG1(IT,IX)
      TG2(IV,IX,IT,IU) = TG2(IT,IU,IV,IX)
    end do
  end do
  ! and then the 3-particle density matrix:
  ! <PSI1|E(T,U,V,X,Y,Z)|PSI2>  = <PSI1|E(TU)E(VX)E(YZ)|PSI2>
  ! -D(Y,X)*(TG2(T,U,V,Z)+D(V,U)*TG1(T,Z))
  ! -D(V,U)*TG2(T,X,Y,Z) C -D(Y,U)*TG2(V,X,T,Z)
  do IP1=1,NASHT**2
    IT = L2ACT(P2LEV(1,IP1))
    IU = L2ACT(P2LEV(2,IP1))
    ITS = IASYM(IT)
    IUS = IASYM(IU)
    IS1 = Mul(Mul(ITS,IUS),LSYM1)
    do IP2=1,IP1
      IV = L2ACT(P2LEV(1,IP2))
      IX = L2ACT(P2LEV(2,IP2))
      IVS = IASYM(IV)
      IXS = IASYM(IX)
      IS2 = Mul(Mul(IVS,IXS),IS1)
      do IP3=1,IP2
        IY = L2ACT(P2LEV(1,IP3))
        IZ = L2ACT(P2LEV(2,IP3))
        IYS = IASYM(IY)
        IZS = IASYM(IZ)
        IS3 = Mul(Mul(IYS,IZS),IS2)
        if (IS3 == LSYM2) then
          call get_tg3_index(IT,IU,IV,IX,IY,IZ,NASHT,jtuvxyz)
          VAL = TG3(JTUVXYZ)
          if (IY == IX) then
            VAL = VAL-TG2(IT,IU,IV,IZ)
            if (IV == IU) VAL = VAL-TG1(IT,IZ)
          end if
          if (IV == IU) VAL = VAL-TG2(IT,IX,IY,IZ)
          if (IY == IU) VAL = VAL-TG2(IV,IX,IT,IZ)
          TG3(JTUVXYZ) = VAL
        end if
      end do
    end do
  end do
  call mma_deallocate(P2LEV)

else

  ! -Special code for the closed-shell or hi-spin cases:
  ! ISCF=1 for closed-shell, =2 for hispin
  OCC = Two
  if (ISCF == 2) OCC = One
  do IT=1,NASHT
    TG1(IT,IT) = OCC
  end do
  if (NACTEL == 1) return
  do IT=1,NASHT
    do IU=1,NASHT
      TG2(IT,IT,IU,IU) = TG1(IT,IT)*TG1(IU,IU)
      if (IU == IT) then
        TG2(IT,IT,IU,IU) = TG2(IT,IT,IU,IU)-TG1(IT,IU)
      else
        TG2(IT,IU,IU,IT) = -TG1(IT,IT)
      end if
    end do
  end do
  if (NACTEL == 2) return
  do IT1=1,NLEV
    do IU1=1,NLEV
      IND1 = IT1+NASHT*(IU1-1)
      do IT2=1,NLEV
        do IU2=1,IU1
          IND2 = IT2+NASHT*(IU2-1)
          if (IND2 > IND1) cycle
          do IT3=1,NLEV
            do IU3=1,IU2
              IND3 = IT3+NASHT*(IU3-1)
              if (IND3 > IND2) cycle
              VAL = TG1(IT1,IU1)*TG1(IT2,IU2)*TG1(IT3,IU3)

              ! Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
              ! Add here the necessary Kronecker deltas times 2-body matrix
              ! elements and lower, so we get a true normal-ordered density matrix
              ! element.

              ! <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
              ! = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
              ! -D(T3,U2)*(TG2(T1,U1,T2,U3)+D(T2,U1)*TG1(T1,U3))
              ! -D(T2,U1)*TG2(T1,U2,T3,U3)
              ! -D(T3,U1)*TG2(T2,U2,T1,U3)

              if (IT3 == IU2) then
                VAL = VAL-TG2(IT1,IU1,IT2,IU3)
                if (IT2 == IU1) VAL = VAL-TG1(IT1,IU3)
              end if
              if (IT2 == IU1) VAL = VAL-TG2(IT1,IU2,IT3,IU3)
              if (IT3 == IU1) VAL = VAL-TG2(IT2,IU2,IT1,IU3)

              ! VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
              ITG3 = nTri3_Elem(IND1-1)+nTri_Elem(IND2-1)+IND3
              TG3(ITG3) = VAL

            end do
          end do
        end do
      end do
    end do
  end do

end if

end subroutine MKTG3
