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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine DERTG3(DOG3,LSYM1,LSYM2,NCONF,NASHT,CI1,CI2,OVL,DTG1,DTG2,NTG3,DTG3,CLAG1,CLAG2)
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
!    ITG3= ((i+1)*i*(i-1))/6 + (j*(j-1))/2 + k
!     i  = max(tu,vx,yz)
!     j  = mid(tu,vx,yz)
!     k  = min(tu,vx,yz)
! tu stands for the pair index tu= t + NASHT*(u-1), etc., and t is
! the usual active orbital number, when they are enumerated across
! all the symmetries (The ''absolute'' active index).

use Symmetry_Info, only: Mul
use sguga, only: SGS, L2ACT, CIS, EXS
use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
use definitions, only: iwp, wp, u6
use caspt2_module, only: NACTEL, ISCF, IASYM
use caspt2_module, only: MXCI
use Constants, only: Zero, One

implicit none
logical(kind=iwp), intent(in) :: DOG3
integer(kind=iwp), intent(in) :: LSYM1, LSYM2, NCONF, NASHT, NTG3
real(kind=wp), intent(in) :: CI1(NCONF), CI2(NCONF), OVL
real(kind=wp), intent(inout) :: DTG1(NASHT,NASHT), DTG2(NASHT,NASHT,NASHT,NASHT), DTG3(NTG3), CLAG1(NCONF), CLAG2(NCONF)
real(kind=wp), allocatable :: TG3WRK(:), BUF1(:), DTU(:,:), DYZ(:,:)
integer(kind=iwp), allocatable :: P2LEV(:)
integer(kind=iwp) :: nLev, LP2LEV1, LP2LEV2, IP, IL, JL, IP1, IT, IU, ITU, ITS, IUS, IS1, IP2, IV, IX, IVX, IVS, IXS, IS2, IP3, &
                     IY, IZ, IYS, IZS, IS3, IYZ, JTU, JVX, JYZ, JTUVXYZ, NCI1, NVECS, NTG3WRK, NYZBUF, NTUBUF, LSGM1, LTAU, LSGM2, &
                     IP3STA, IP3END, LTO, ISSG2, IP1STA, IP1END, ISSG1, LFROM, LFROMD, IM, JM, ISTAU, NTAU, ibuf
real(kind=wp) :: VAL

nLev = SGS%nLev

! Put in zeroes. Recognize special cases:
!OVL = One
if (NASHT == 0) return
if (NACTEL == 0) return
!if (LSYM1 /= LSYM2) OVL = Zero
!TG1(:,:) = Zero
!TG2(:,:,:,:) = Zero
!TG3(1:NTG3) = Zero

if (ISCF /= 0) then
  ! -Special code for the closed-shell or hi-spin cases:
  ! ISCF=1 for closed-shell, =2 for hispin
  write(u6,*) 'Here is the special case'
  write(u6,*) 'not yet'
  call abend()
end if
! Here, for regular CAS or RAS cases.

! Special pair index allows true RAS cases to be handled:
call mma_allocate(P2LEV,2*NASHT**2,Label='P2LEV')
LP2LEV1 = 1
LP2LEV2 = 1+NASHT**2
IP = 0
! First, IL < JL pairs.
do IL=1,NLEV-1
  do JL=IL+1,NLEV
    IP = IP+1
    P2LEV(LP2LEV1-1+IP) = IL
    P2LEV(LP2LEV2-1+IP) = JL
  end do
end do
! Then, IL = JL pairs.
do IL=1,NLEV
  IP = IP+1
  P2LEV(LP2LEV1-1+IP) = IL
  P2LEV(LP2LEV2-1+IP) = IL
end do
! Last, IL > JL pairs.
do IL=2,NLEV
  do JL=1,IL-1
    IP = IP+1
    P2LEV(LP2LEV1-1+IP) = IL
    P2LEV(LP2LEV2-1+IP) = JL
  end do
end do

! First, the 3-particle density matrix:
! <PSI1|E(T,U,V,X,Y,Z)|PSI2>  = <PSI1|E(TU)E(VX)E(YZ)|PSI2>
! -D(Y,X)*(TG2(T,U,V,Z)+D(V,U)*TG1(T,Z))
! -D(V,U)*TG2(T,X,Y,Z) C -D(Y,U)*TG2(V,X,T,Z)
if (DOG3) then
  do IP1=1,NASHT**2
    IT = L2ACT(P2LEV(LP2LEV1-1+IP1))
    IU = L2ACT(P2LEV(LP2LEV2-1+IP1))
    ITU = IT+NASHT*(IU-1)
    ITS = IASYM(IT)
    IUS = IASYM(IU)
    IS1 = Mul(Mul(ITS,IUS),LSYM1)
    do IP2=1,IP1
      IV = L2ACT(P2LEV(LP2LEV1-1+IP2))
      IX = L2ACT(P2LEV(LP2LEV2-1+IP2))
      IVX = IV+NASHT*(IX-1)
      IVS = IASYM(IV)
      IXS = IASYM(IX)
      IS2 = Mul(Mul(IVS,IXS),IS1)
      do IP3=1,IP2
        IY = L2ACT(P2LEV(LP2LEV1-1+IP3))
        IZ = L2ACT(P2LEV(LP2LEV2-1+IP3))
        IYS = IASYM(IY)
        IZS = IASYM(IZ)
        IS3 = Mul(Mul(IYS,IZS),IS2)
        if (IS3 == LSYM2) then
          IYZ = IY+NASHT*(IZ-1)
          if (ITU < IVX) then
            if (ITU >= IYZ) then
              JTU = IVX
              JVX = ITU
              JYZ = IYZ
            else if (IVX < IYZ) then
              JTU = IYZ
              JVX = IVX
              JYZ = ITU
            else
              JTU = IVX
              JVX = IYZ
              JYZ = ITU
            end if
          else
            if (ITU < IYZ) then
              JTU = IYZ
              JVX = ITU
              JYZ = IVX
            else if (IVX >= IYZ) then
              JTU = ITU
              JVX = IVX
              JYZ = IYZ
            else
              JTU = ITU
              JVX = IYZ
              JYZ = IVX
            end if
          end if
          JTUVXYZ = ((JTU+1)*JTU*(JTU-1))/6+(JVX*(JVX-1))/2+JYZ
          VAL = DTG3(JTUVXYZ)
          if (IY == IX) then
            DTG2(IT,IU,IV,IZ) = DTG2(IT,IU,IV,IZ)-VAL
            if (IV == IU)  DTG1(IT,IZ) = DTG1(IT,IZ)-VAL
          end if
          if (IV == IU) DTG2(IT,IX,IY,IZ) = DTG2(IT,IX,IY,IZ)-VAL
          if (IY == IU) DTG2(IV,IX,IT,IZ) = DTG2(IV,IX,IT,IZ)-VAL
        end if
      end do
    end do
  end do
end if

! Then, the 2-particle density matrix:
! <PSI1|E(T,U,V,X)|PSI2>  = <PSI1|E(TU)E(VX)|PSI2> - D(V,U)*TG2(T,U,V,X)
do IP1=1,NASHT**2
  IT = L2ACT(P2LEV(LP2LEV1-1+IP1))
  IU = L2ACT(P2LEV(LP2LEV2-1+IP1))
  do IP2=1,IP1
    IV = L2ACT(P2LEV(LP2LEV1-1+IP2))
    IX = L2ACT(P2LEV(LP2LEV2-1+IP2))
    if (IP1 /= IP2) then
      DTG2(IT,IU,IV,IX) = DTG2(IT,IU,IV,IX)+DTG2(IV,IX,IT,IU)
      DTG2(IV,IX,IT,IU) = Zero
    end if
    if (IV == IU) DTG1(IT,IX) = DTG1(IT,IX)-DTG2(IT,IU,IV,IX)
  end do
end do

! If now any matrix element E(t1u1)E(t2u2)..E(tnun) is arranged
! such that the pair indices are non-decreasing, then the matrix
! element can be correctly computed by performing explicit
! excitations within the RAS space.
! But we also need the 'usual' pair index in order to use the
! packed addressing.

NCI1 = CIS%NCSF(LSYM1)
! Overlap:
!if (LSYM1 == LSYM2) OVL = DDOT_(NCI1,CI1,1,CI2,1)
if (LSYM1 == LSYM2) then
  CLAG2(1:NCI1) = CLAG2(1:NCI1)+OVL*CI1(1:NCI1)
  CLAG1(1:NCI1) = CLAG1(1:NCI1)+OVL*CI2(1:NCI1)
end if
!write(u6,*) 'overlap = ',DDOT_(NCI1,CI1,1,CI2,1)
! Allocate as many vectors as possible:
! Wishful thinking:
NVECS = 2*NASHT**2+1
! But what is really available?
call mma_MaxDBLE(NTG3WRK)
NTG3WRK = NTG3WRK-3*MXCI ! for BUF1, DTU, and DYZ, allocated later
NTG3WRK = NTG3WRK/2
NTG3WRK = min(MXCI*NVECS,NTG3WRK)
NVECS = NTG3WRK/MXCI
NTG3WRK = NVECS*MXCI
! Find optimal subdivision of available vectors:
NYZBUF = nint(real(NVECS-1,kind=wp)/real(NASHT,kind=wp),kind=iwp)
NYZBUF = max(1,NYZBUF)
NTUBUF = min(NASHT**2,NVECS-1-NYZBUF)
NYZBUF = NVECS-1-NTUBUF
! Insufficient memory?
if (NTUBUF <= 0) then
  write(u6,*) ' Too little memory left for MKTG3.'
  write(u6,*) ' Need at least 6 vectors of length MXCI=',MXCI
  call ABEND()
end if
if (NTUBUF <= (NASHT**2)/5) then
  write(u6,*) ' WARNING: MKTG3 will be inefficient owing to'
  write(u6,*) ' small memory.'
end if
call mma_allocate(TG3WRK,NTG3WRK,Label='TG3WRK')
call mma_allocate(BUF1,MXCI,Label='BUF1')

call mma_allocate(DTU,MXCI,NTUBUF,Label='DTU')
call mma_allocate(DYZ,MXCI,NYZBUF,Label='DYZ')

! And divide it up:
!! LSGM1: NTUBUF vectors
!! LTAU : 1 vector
!! LSGM2: NYZBUF vectors
LSGM1 = 1
LTAU = LSGM1+NTUBUF*MXCI
LSGM2 = LTAU+MXCI

! Sectioning loops over pair indices IP3 (ket side):
do IP3STA=1,NASHT**2,NYZBUF
  IP3END = min(NASHT**2,IP3STA-1+NYZBUF)
  ! Compute a section of sigma vectors E(YZ)*PSI2 to memory:
  LTO = LSGM2
  do IP3=IP3STA,IP3END
    ! Translate to levels in the SGUGA coupling order:
    IL = P2LEV(LP2LEV1-1+IP3)
    JL = P2LEV(LP2LEV2-1+IP3)
    IY = L2ACT(IL)
    IZ = L2ACT(JL)
    IYS = IASYM(IY)
    IZS = IASYM(IZ)
    ISSG2 = Mul(Mul(IYS,IZS),LSYM2)
    TG3WRK(LTO:LTO+MXCI-1) = Zero
    ! LTO is first element of Sigma2 = E(YZ) Psi2
    call SG_Epq_Psi(SGS,CIS,EXS,IL,JL,One,LSYM2,CI2,TG3WRK(LTO))
    !! It is possible to calculate the contribution using
    !! DGEMV, but DAXPY seems to be faster than DGEMV
    if ((ISSG2 == LSYM1) .and. (DTG1(IY,IZ) /= Zero)) CLAG1(1:NCI1) = CLAG1(1:NCI1)+DTG1(IY,IZ)*TG3WRK(LTO:LTO+NCI1-1)
    LTO = LTO+MXCI
  end do

  DYZ(1:MXCI,1:NYZBUF) = Zero
  ! Sectioning loops over pair indices IP1 (bra side):
  do IP1STA=IP3STA,NASHT**2,NTUBUF
    IP1END = min(NASHT**2,IP1STA-1+NTUBUF)
    ! Compute a section of sigma vectors E(UT)*PSI1 to memory:
    LTO = LSGM1
    !! <Psi1|E(TU)
    do IP1=IP1STA,IP1END
      ! Translate to levels:
      JL = P2LEV(LP2LEV1-1+IP1)
      IL = P2LEV(LP2LEV2-1+IP1)
      IT = L2ACT(IL)
      IU = L2ACT(JL)
      ITS = IASYM(IT)
      IUS = IASYM(IU)
      ISSG1 = Mul(Mul(ITS,IUS),LSYM1)
      TG3WRK(LTO:LTO+MXCI-1) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,IL,JL,One,LSYM1,CI1,TG3WRK(LTO))
      if ((ISSG1 == LSYM1) .and. (DTG1(IU,IT) /= Zero) .and. (IP3STA == 1)) &
        CLAG2(1:NCI1) = CLAG2(1:NCI1)+DTG1(IU,IT)*TG3WRK(LTO:LTO+NCI1-1)
      LTO = LTO+MXCI
    end do

    DTU(1:MXCI,1:NTUBUF) = Zero
    ! Now compute as many elements as possible:
    LFROM = LSGM2
    LFROMD = 1
    do IP3=IP3STA,IP3END
      IY = L2ACT(P2LEV(LP2LEV1-1+IP3))
      IZ = L2ACT(P2LEV(LP2LEV2-1+IP3))
      ! LFROM will be start element of Sigma2=E(YZ) Psi2
      IYZ = IY+NASHT*(IZ-1)
      IYS = IASYM(IY)
      IZS = IASYM(IZ)
      ISSG2 = Mul(Mul(IYS,IZS),LSYM2)
      IM = P2LEV(LP2LEV1-1+IP3)
      JM = P2LEV(LP2LEV2-1+IP3)
      do IP2=IP3,IP1END
        IL = P2LEV(LP2LEV1-1+IP2)
        JL = P2LEV(LP2LEV2-1+IP2)
        IV = L2ACT(IL)
        IX = L2ACT(JL)
        IVX = IV+NASHT*(IX-1)
        IVS = IASYM(IV)
        IXS = IASYM(IX)
        ISTAU = Mul(Mul(IVS,IXS),ISSG2)
        NTAU = CIS%NCSF(ISTAU)
        TG3WRK(LTAU:LTAU+MXCI-1) = Zero
        ! LTAU  will be start element of Tau=E(VX) Sigma2=E(VX) E(YZ) Psi2
        !! LTAU = EvxEyz|Psi2>
        call SG_Epq_Psi(SGS,CIS,EXS,IL,JL,One,ISSG2,TG3WRK(LFROM),TG3WRK(LTAU))
        if ((ISTAU == LSYM1) .and. (DTG2(IV,IX,IY,IZ) /= Zero)) then
          !DTG2(IV,IX,IY,IZ) = DDOT_(NTAU,TG3WRK(LTAU),1,CI1,1)
          !! For left derivative: <I|Evx Eyz|Psi2>
          CLAG1(1:NTAU) = CLAG1(1:NTAU)+DTG2(IV,IX,IY,IZ)*TG3WRK(LTAU:LTAU+NTAU-1)
          !! For right derivative: <Psi1|Evx Eyz|I>
          if ((IP2 >= IP1STA) .and. (IP2 <= IP1END)) then
            ibuf = lsgm1+mxci*(ip2-ip1sta)
            DYZ(1:MXCI,LFROMD) = DYZ(1:MXCI,LFROMD)+DTG2(IV,IX,IY,IZ)*TG3WRK(IBUF:IBUF+MXCI-1)
          else
            call SG_Epq_Psi(SGS,CIS,EXS,JL,IL,DTG2(IV,IX,IY,IZ),ISSG2,CI1,DYZ(1,LFROMD))
          end if
          DTG2(IV,IX,IY,IZ) = Zero
        end if
        if (DOG3) then
          BUF1(1:MXCI) = Zero
          do IP1=max(IP2,IP1STA),IP1END
            IT = L2ACT(P2LEV(LP2LEV1-1+IP1))
            IU = L2ACT(P2LEV(LP2LEV2-1+IP1))
            ITS = IASYM(IT)
            IUS = IASYM(IU)
            ISSG1 = Mul(Mul(ITS,IUS),LSYM1)
            if (ISSG1 == ISTAU) then
              !L = LSGM1+MXCI*(IP1-IP1STA)
              !VAL = DDOT_(NTAU,TG3WRK(LTAU),1,TG3WRK(L),1)
              ITU = IT+NASHT*(IU-1)
              ! Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
              ! Code to put it in correct place:
              if (ITU < IVX) then
                if (ITU >= IYZ) then
                  JTU = IVX
                  JVX = ITU
                  JYZ = IYZ
                else if (IVX < IYZ) then
                  JTU = IYZ
                  JVX = IVX
                  JYZ = ITU
                else
                  JTU = IVX
                  JVX = IYZ
                  JYZ = ITU
                end if
              else
                if (ITU < IYZ) then
                  JTU = IYZ
                  JVX = ITU
                  JYZ = IVX
                else if (IVX >= IYZ) then
                  JTU = ITU
                  JVX = IVX
                  JYZ = IYZ
                else
                  JTU = ITU
                  JVX = IYZ
                  JYZ = IVX
                end if
              end if
              JTUVXYZ = ((JTU+1)*JTU*(JTU-1))/6+(JVX*(JVX-1))/2+JYZ
              if (DTG3(JTUVXYZ) /= Zero) then
                !! For left derivative: <I|Evx Eyz|Psi2> * Dtuvxyz
                !! I don't understand, but this is much faster than
                !! processing all possible vectors at once with DGEMM
                !! (and DGER) after finishing the IP1 loop below
                DTU(1:MXCI,1+IP1-IP1STA) = DTU(1:MXCI,1+IP1-IP1STA)+DTG3(JTUVXYZ)*TG3WRK(LTAU:LTAU+MXCI-1)
                !! For right derivative: <Psi1|Etu|I> * Dtuvxyz
                !! This is also (slightly) faster than DGEMV, apparently
                call DaXpY_(MXCI,DTG3(JTUVXYZ),TG3WRK(LSGM1+MXCI*(IP1-IP1STA)),1,BUF1,1)
              end if
              ! End of symmetry requirement IF-clause:
            end if
            ! End of IP1 loop.
          end do
          !! Second operator for the right derivative:
          !! <Psi1|Etu Evx|I> * Dtuvxyz
          call SG_Epq_Psi(SGS,CIS,EXS,JL,IL,One,ISTAU,BUF1,DYZ(1,LFROMD))
        end if !! End of DOG3 clause
        ! End of IP2 loop.
      end do
      LFROM = LFROM+MXCI
      LFROMD = LFROMD+1
      ! End of IP3 loop.
    end do

    LTO = 1
    !! <I|Etu Evx Eyz|Psi2> * Dtuvxyz
    do IP1=IP1STA,IP1END
      ! Translate to levels:
      IL = P2LEV(LP2LEV1-1+IP1)
      JL = P2LEV(LP2LEV2-1+IP1)
      IT = L2ACT(IL)
      IU = L2ACT(JL)
      ITS = IASYM(IT)
      IUS = IASYM(IU)
      ISSG1 = Mul(Mul(ITS,IUS),LSYM1)
      call SG_Epq_Psi(SGS,CIS,EXS,IL,JL,One,LSYM1,DTU(1,LTO),CLAG1)
      LTO = LTO+1
    end do
    ! End of IP1STA sectioning loop
  end do

  LTO = 1
  do IP3=IP3STA,IP3END
    IY = L2ACT(P2LEV(LP2LEV1-1+IP3))
    IZ = L2ACT(P2LEV(LP2LEV2-1+IP3))
    ! LFROM will be start element of Sigma2=E(YZ) Psi2
    IYZ = IY+NASHT*(IZ-1)
    IYS = IASYM(IY)
    IZS = IASYM(IZ)
    ISSG2 = Mul(Mul(IYS,IZS),LSYM2)
    IM = P2LEV(LP2LEV1-1+IP3)
    JM = P2LEV(LP2LEV2-1+IP3)
    ! LTO is first element of Sigma2 = E(YZ) Psi2
    call SG_Epq_Psi(SGS,CIS,EXS,JM,IM,One,LSYM2,DYZ(1,LTO),CLAG2)
    LTO = LTO+1
  end do
  ! End of IP3STA sectioning loop
end do

call mma_deallocate(TG3WRK)
call mma_deallocate(BUF1)
call mma_deallocate(DTU)
call mma_deallocate(DYZ)
call mma_deallocate(P2LEV)

end subroutine DERTG3
