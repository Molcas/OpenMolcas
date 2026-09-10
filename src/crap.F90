subroutine SG_Epq_Psi_Flat(SGS,CIS,EXS,IP,IQ,ISYCI)

use Symmetry_Info, only: Mul
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout), target :: EXS
integer(kind=iwp), intent(in) :: IP, IQ, ISYCI
integer(kind=iwp) :: I, IC, ICS, INDEO, IOC, IOLW, IOUW, IPPOW, IPSHFT, ISGSTA, ISTA, ISYDC, ISYDSG, ISYP, ISYPQ, ISYQ, ISYSGM, &
                     ISYUC, ISYUSG, J, JC, JSTA, LICP, LLW, LUW, MVSGM, NCP, NDWNC, NDWNSG, NS1, NTMP, NUPC, &
                     NUPSG, NCP1, NCP2, MV, MVX, nCSFs, LICP1, LICP2
real(kind=wp) :: X

! declarations to facilitate the reuse of sigma vectors if possible
integer(kind=iwp), save:: i_save_p=0, i_save_q=0
integer(kind=iwp), save:: i_save_p_sym=-1, i_save_q_sym=-1
logical(kind=iwp) :: Reuse_Sigma, Swap, Test1, Test2, Test3
real(kind=wp) :: CI_ID=Zero
integer(kind=iwp) ::  iOff, jOff, nEntry

!*********************************************************************
!  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
!  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
!  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
!  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
!  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
!*********************************************************************

nEntry=0
nCSFs=CIS%NCSF(ISYCI)
! SYMMETRY OF ORBITALS:
ISYP = SGS%ISm(IP)
ISYQ = SGS%ISm(IQ)
ISYPQ = Mul(ISYP,ISYQ)
! SYMMETRY OF SIGMA ARRAY:
ISYSGM = Mul(ISYPQ,ISYCI)


if (IQ/=IP) then

! Unified code for the case P/=Q
Swap = IQ > IP
Test1 = MERGE(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
Test2 = MERGE(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)

if (Test1) then

  ! EXCITING CASE, IQ<IP<=MIDLEV.
  INDEO = MERGE(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
  do MVSGM=1,CIS%nMidV
    do ISYUSG=1,SGS%nSym

      NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
      if (NS1 == 0) cycle
      ISYDSG = Mul(ISYUSG,ISYSGM)           ! compute the lower symmetry
      ISYDC = Mul(ISYPQ,ISYDSG)             ! compute the symmetry of the sigma vector
      NDWNC = CIS%NOW(2,ISYDC,MVSGM)
      if (NDWNC == 0) cycle
      ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM) ! get the off-set to the sigma vector block
      NUPSG = CIS%NOW(1,ISYUSG,MVSGM)           ! number of upper half-walks by symmetry and midvertex.
      IOC = CIS%IOCSF(ISYUSG,MVSGM,ISYCI)       ! get the off-set to the CI vector block
      NCP = MERGE(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVSGM),Swap)
      if (NCP == 0) cycle
      NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)

      ! get the off-set to the block of compressed coupling coefficients.
      LICP = MERGE(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVSGM),Swap)
      ! CASE IS: LOWER HALF, EXCITE:
      nEntry=nEntry+1

    end do
  end do
else if (Test2) then
  ! EXCITING CASE, MIDLEV<IQ<IP
  INDEO = MERGE(2*SGS%nLev+(IQ*(IQ-1))/2+IP,2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
  do MVSGM=1,CIS%nMidV
    do ISYUSG=1,SGS%nSym
      NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
      if (NS1 == 0) cycle
      ISYUC = Mul(ISYPQ,ISYUSG)
      NUPC = CIS%NOW(1,ISYUC,MVSGM)
      if (NUPC == 0) cycle
      ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
      NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
      ISYDSG = Mul(ISYUSG,ISYSGM)
      NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
      IOC = CIS%IOCSF(ISYUC,MVSGM,ISYCI)
      NCP = MERGE(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVSGM),Swap)
      if (NCP == 0) cycle
      LICP = MERGE(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVSGM),Swap)
      ! CASE IS: UPPER HALF, EXCITE:
      nEntry=nEntry+1
    end do
  end do
else
  ! EXCITING CASE, IQ<=MIDLEV<IP

  iOff=1  ! Initialize pointer to the SGTMP scratch area
  Reuse_Sigma = .false.
  ! Check if there is a potential oppertunity to reuse the vector for the first part of the
  ! coupling coefficient. That is we are repeating the first index.
  Test3=MERGE((i_save_q == IQ) .and. (i_save_p_sym == ISYP), &
              (i_save_p == IP) .and. (i_save_q_sym == ISYQ), Swap)
  ! If we have opted in for reuse, and there is the oppertunity for reuse, then
  ! check that the CI vector has not changed since last time. If to we set Reuse_Sigma=.TRUE.
  if (EXS%Reuse_SGTMP .and.  Test3) Reuse_Sigma = CI_ID == sum(CI(1:nCSFs))+real(nCSFs,kind=wp)

  do MVSGM=1,CIS%nMidV
    do MV = 1, 2
      MVX = MERGE(EXS%MVR(MVSGM,MV),EXS%MVL(MVSGM,MV),SWAP)
      if (MVX == 0) cycle

      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        ISYUC = MERGE(Mul(ISYQ,ISYUSG),Mul(ISYP,ISYUSG),Swap)
        ISYDC = MERGE(Mul(ISYP,ISYDSG),Mul(ISYQ,ISYDSG),Swap)

        NUPC = CIS%NOW(1,ISYUC,MVX)
        if (NUPC == 0) cycle
        NDWNC = CIS%NOW(2,ISYDC,MVX)
        if (NDWNC == 0) cycle

        INDEO = MERGE(merge(IQ, SGS%nLev+IQ, MV==1),merge(IP, SGS%nLev+IP, MV==1),Swap)
        NCP1 = MERGE(EXS%NOCP(INDEO,ISYUSG,MVSGM),EXS%NOCP(INDEO,ISYUC,MVX),Swap)
        if (NCP1 == 0) cycle

        ! CASE IS: UPPER HALF, EXCITE:
        LICP1 = MERGE(EXS%IOCP(INDEO,ISYUSG,MVSGM),EXS%IOCP(INDEO,ISYUC,MVX),Swap)
        IOC = CIS%IOCSF(ISYUC,MVX,ISYCI)

        ! If we have opted in for reuse and there is no vectors stored yet, then construct it.
        ! Store it at position iOff.
        if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) then
           NTMP = NUPSG*NDWNC
           EXS%SGTMP(iOff:iOff+NTMP-1) = Zero
           call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(iOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
        end if

        ! If we are to reuse we need to increment the offset for the next set, jOff is set to iOff
        jOff=iOff
        if (EXS%Reuse_SGTMP) then
          iOff = iOff + NUPSG*NDWNC
          if (iOff-1 > size(EXS%SGTMP)) then
            write(u6,*) 'iOff-1>SIZE(EXS%SGTMP) -- 1'
            write(u6,*) iOff-1,size(EXS%SGTMP)
            call Abend()
          end if
        end if

        INDEO = MERGE(merge(IP, SGS%nLev+IP, MV==1),merge(IQ, SGS%nLev+IQ, MV==1),Swap)
        NCP2 = MERGE(EXS%NOCP(INDEO,ISYDSG,MVSGM),EXS%NOCP(INDEO,ISYDC,MVX),Swap)
        if (NCP2 == 0) cycle

        if (.not. EXS%Reuse_SGTMP) then
           NTMP = NUPSG*NDWNC
           EXS%SGTMP(jOff:jOff+NTMP-1) = Zero
           call Apply_row(One,NDWNC,NUPC,CI(IOC+1),NUPSG,EXS%SGTMP(jOff),NCP1,EXS%ICOUP(1,LICP1+1),swap=Swap)
        end if

        ! CASE IS: LOWER HALF, EXCITE:
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        LICP2 = MERGE(EXS%IOCP(INDEO,ISYDSG,MVSGM),EXS%IOCP(INDEO,ISYDC,MVX),Swap)
        nEntry = nEntry + 1
        call Apply_col(CPQ,NUPSG,NDWNC,EXS%SGTMP(jOff),NDWNSG,SGM(ISGSTA+1),NCP2,EXS%ICOUP(1,LICP2+1),swap=Swap)

      end do
    end do
  end do

  If (SWAP) Then
    i_save_q=IQ
    i_save_p_sym=ISYP
    i_save_p=0
    i_save_q_sym=-1
  else
    i_save_p=IP
    i_save_q_sym=ISYQ
    i_save_q=0
    i_save_p_sym=-1
  end if
  if (EXS%Reuse_SGTMP .and. (.not. Reuse_Sigma)) CI_ID = sum(CI(1:nCSFs))+real(nCSFs,kind=wp)

end if

else

  ! THEN THIS IS A SPECIAL CASE.
  if (IP > SGS%MidLev) then

    ! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
    ! IP=IQ>MIDLEV
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOUW = CIS%IOW(1,ISYUSG,MVSGM)
        IPSHFT = 2*(IP-1-SGS%MidLev)
        LUW = 1+IOUW-CIS%nIpWlk+IPSHFT/(2*nPack)
        IPSHFT = mod(IPSHFT,2*nPack)
        IPPOW = 2**IPSHFT
        do I=1,NUPSG
          IC = CIS%ICase(LUW+I*CIS%nIpWlk)
          ICS = mod(IC/IPPOW,4)
          if (ICS == 0) cycle
          X = CPQ*real((1+ICS)/2,kind=wp)
          ISTA = ISGSTA+I
          call DAXPY_(NDWNSG,X,CI(ISTA),NUPSG,SGM(ISTA),NUPSG)
        end do
      end do
    end do

  else

    ! SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
    ! IP=IQ < MIDLEV.
    do MVSGM=1,CIS%nMidV
      do ISYUSG=1,SGS%nSym
        NS1 = CIS%NOCSF(ISYUSG,MVSGM,ISYSGM)
        if (NS1 == 0) cycle
        ISGSTA = CIS%IOCSF(ISYUSG,MVSGM,ISYSGM)
        NUPSG = CIS%NOW(1,ISYUSG,MVSGM)
        ISYDSG = Mul(ISYUSG,ISYSGM)
        NDWNSG = CIS%NOW(2,ISYDSG,MVSGM)
        IOLW = CIS%IOW(2,ISYDSG,MVSGM)
        IPSHFT = 2*(IP-1)
        LLW = 1+IOLW-CIS%nIpWlk+IPSHFT/(2*nPack)
        IPSHFT = mod(IPSHFT,2*nPack)
        IPPOW = 2**IPSHFT
        do J=1,NDWNSG
          JC = CIS%ICase(LLW+J*CIS%nIpWlk)
          ICS = mod(JC/IPPOW,4)
          if (ICS == 0) cycle
          X = CPQ*real((1+ICS)/2,kind=wp)
          JSTA = ISGSTA + NUPSG*(J-1) + 1
          SGM(JSTA:JSTA+NUPSG-1) = SGM(JSTA:JSTA+NUPSG-1)+X*CI(JSTA:JSTA+NUPSG-1)
        end do
      end do
    end do

  end if
  i_save_q=0
  i_save_p_sym=-1
  i_save_p=0
  i_save_q_sym=-1
  CI_ID=Zero
end if

subroutine SG_Epq_Psi(SGS,CIS,EXS,IP,IQ,CPQ,ISYCI,CI,SGM)
