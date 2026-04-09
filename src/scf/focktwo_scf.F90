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
! Copyright (C) Markus P. Fuelscher                                    *
!               1992, Per Ake Malmqvist                                *
!               2002,2023, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine FOCKTWO_scf(NSYM,NBAS,NFRO,KEEP,DLT,DSQ,FLT,nFlt,FSQ,X1,nX1,X2,nX2,ExFac,nD,nBSQT)

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use RICD_Info, only: Do_DCCD
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, NBAS(nSym), NFRO(nSym), KEEP(nSym), nFlt, nX1, nX2, nD, nBSQT
real(kind=wp), intent(in) :: DLT(nFlt,nD), DSQ(nBSQT,nD), ExFac
real(kind=wp), intent(inout) :: FLT(nFlt,nD), FSQ(nBSQT,nD)
real(kind=wp), intent(out) :: X1(nX1), X2(nX2)
integer(kind=iwp) :: IB, IJB, IJS, IK, iOpt, IRC, ISD, ISF, ISTLT(8), ISTSQ(8), ISX, ISYM, JB, JK, K1, K2, KB, KK, KLB, LB, LK, &
                     LPQ, LSMAX, NB, NB2, NB3, NFI, NFJ, NFK, NFL, NPQ
real(kind=wp) :: Factor, temp, temp_ab
real(kind=wp), external :: DDot_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ISTLTT
#endif

! This routine has been nicked from the MOTRA package. It was
! originally written by Marcus Fuelscher, and has been slightly
! modified by P-A Malmqvist 1992-12-04.
! Further modifications by RL 2002-8-30.
! Purpose: Return FLT, which is the effective one-electron
! Hamiltonian when frozen orbitals are removed. It is the
! same as a closed-shell Fock matrix computed using the
! two-electron contribution only from frozen orbitals.
! FLT is returned as symmetry-blocked lower triangles. FSQ
! contains the same matrix, as symmetry-blocked square matrices.
! DSQ and DLT are computed in the calling routine.
! It is assumed that the SEWARD integral file was opened before
! call to this routine.
! ISTSQ, ISTLT: Offsets to symmetry blocks of FSQ, FLT etc.
!                                                                      *
!***********************************************************************
!                                                                      *
if (Do_DCCD .and. (NSYM /= 1)) then
  write(u6,*) 'DCCD not implemented for nSym/=1'
  call Abend()
end if

Factor = real(nD,kind=wp)*Half
ISTSQ(:) = 0
ISTLT(:) = 0

if (Do_DCCD) then
  if (NSYM /= 1) then
    write(u6,*) 'DCCD not implemented for nSym/=1'
    call Abend()
  end if
  call FOCKTWO_scf_DCCD()
else
  if (NSYM == 1) then
    call FOCKTWO_scf_NoSym()
  else
    call FOCKTWO_scf_Sym()
  end if
end if

if (IRC /= 0) then
  write(u6,*) ' Error return code IRC=',IRC
  write(u6,*) ' from RDORD call, in FTWOI.'
  call Abend()
end if

! Accumulate the contributions
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  K1 = ISTLT(ISYM)
  K2 = ISTSQ(ISYM)
  do IB=1,NB
    do JB=1,IB

      FLT(K1+JB,1) = FLT(K1+JB,1)+FSQ(K2+JB,1)
      if (nD == 2) FLT(K1+JB,2) = FLT(K1+JB,2)+FSQ(K2+JB,2)
#     ifdef _DEBUGPRINT_
      if (nD == 1) then
        write(u6,'(a,i5,a,f12.6)') 'Flt(',K1+JB,',1)=',FLT(K1+JB,1)
      else
        write(u6,'(a,i5,a,2f12.6)') 'Flt(',K1+JB,',:)=',FLT(K1+JB,1),FLT(K1+JB,2)
      end if
#     endif

    end do
    K1 = K1+IB
    K2 = K2+NB
  end do
end do

call GADSum(Flt,nFlt*nD)

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
write(u6,'(6X,A)') 'TEST PRINT FROM FTWOI.'
write(u6,'(6X,A)') 'FROZEN FOCK MATRIX IN AO BASIS:'
ISTLTT = 1
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB > 0) then
    write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FLT(ISTLTT,1),NB)
    if (nD == 2) call TRIPRT(' ',' ',FLT(ISTLTT,2),NB)
    ISTLTT = ISTLTT+nTri_Elem(NB)
  end if
end do
write(u6,'(6X,A)') '----------------------------'
#endif

contains

subroutine FOCKTWO_scf_Sym()

  integer(kind=iwp) :: IP, IPQ, IS, ISYM, JQ, JS, KS, LS
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: ivv
# endif

  do ISYM=2,NSYM
    NB = NBAS(ISYM-1)
    NB2 = NB*NB
    NB3 = (NB2+NB)/2
    ISTSQ(ISYM) = ISTSQ(ISYM-1)+NB2
    ISTLT(ISYM) = ISTLT(ISYM-1)+NB3
  end do

  ! Loop over the symmetry blocks (IS,JS|KS,LS)

  do IS=1,NSYM
    IB = NBAS(IS)
    IK = KEEP(IS)
    NFI = NFRO(IS)
    do JS=1,IS
      JB = NBAS(JS)
      JK = KEEP(JS)
      NFJ = NFRO(JS)
      IJS = Mul(IS,JS)
      IJB = IB*JB
      if (IS == JS) IJB = nTri_Elem(IB)
      do KS=1,IS
        KB = NBAS(KS)
        KK = KEEP(KS)
        NFK = NFRO(KS)
        LSMAX = KS
        if (KS == IS) LSMAX = JS
        LS = Mul(IJS,KS)
        if (LS > LSMAX) cycle
        LB = NBAS(LS)
        LK = KEEP(LS)
        NFL = NFRO(LS)
        KLB = KB*LB
        if (KS == LS) KLB = nTri_Elem(KB)
        ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
        if ((IK+JK+KK+LK) /= 0) cycle
        ! NO FROZEN ORBITALS?
        if ((NFI+NFJ+NFK+NFL) == 0) cycle
        ! NO BASIS FUNCTIONS?
        if ((IJB*KLB) == 0) cycle

        ! Process the different symmetry cases

        if ((IS == JS) .and. (IS == KS)) then
          ! CASE 1: Integrals are of symmetry type (II/II)
          ! Coulomb and exchange terms need to be accumulated
          ! Option code 1: Begin reading at first integral.
          ! NPQ: Nr of submatrices in buffer X1.
          IOPT = 1
          LPQ = 0
          IPQ = 0
          NPQ = 0
          do IP=1,IB
            do JQ=1,IP
              IPQ = IPQ+1
              LPQ = LPQ+1
              if (IPQ > NPQ) then
                call RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
                if (IRC > 1) return
                ! Option code 2: Continue reading at next integral.
                IOPT = 2
                IPQ = 1
              end if
              ISX = (IPQ-1)*KLB+1
              ISF = ISTLT(IS)+LPQ
              ISD = ISTLT(IS)+1
              TEMP = DDOT_(KLB,X1(ISX),1,DLT(ISD,1),1)
              FLT(ISF,1) = FLT(ISF,1)+TEMP
              if (nD == 2) then
                TEMP_ab = DDOT_(KLB,X1(ISX),1,DLT(ISD,2),1)
                FLT(ISF,1) = FLT(ISF,1)+TEMP_ab
                FLT(ISF,2) = FLT(ISF,1)
              end if
#             ifdef _DEBUGPRINT_
              write(u6,'(a,i5,a,f12.6)') '00 Flt(',isf,',1)=',FLT(ISF,1)
              if (nD == 2) write(u6,'(a,i5,a,f12.6)') '00 Flt(',isf,',2)=',FLT(ISF,2)
#             endif
              call SQUARE(X1(ISX),X2,1,KB,LB)
              ISF = ISTSQ(IS)+(JQ-1)*JB+1
              ISD = ISTSQ(IS)+(IP-1)*IB+1

              call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
              if (nD == 2) call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
              if (IP /= JQ) then
                ISF = ISTSQ(IS)+(IP-1)*IB+1
                ISD = ISTSQ(IS)+(JQ-1)*JB+1

                call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                if (nD == 2) call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
              end if
#             ifdef _DEBUGPRINT_
              write(u6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
              if (nD == 2) write(u6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
#             endif

            end do  ! JQ
          end do    ! IP

        else if ((IS == JS) .and. (IS /= KS)) then
          ! CASE 2: Integrals are of symmetry type (II/JJ)
          ! Coulomb terms need to be accumulated only
          IOPT = 1
          LPQ = 0
          IPQ = 0
          NPQ = 0
          do IP=1,IB
            do JQ=1,IP
              IPQ = IPQ+1
              LPQ = LPQ+1
              if (IPQ > NPQ) then
                call RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
                if (IRC > 1) return
                IOPT = 2
                IPQ = 1
              end if
              ISX = (IPQ-1)*KLB+1
              if (NFI /= 0) then
                ISF = ISTLT(KS)+1
                ISD = ISTLT(IS)+LPQ
                TEMP = DLT(ISD,1)
                if (nD == 2) TEMP = DLT(ISD,1)+DLT(ISD,2)
                FLT(ISF:ISF+KLB-1,1) = FLT(ISF:ISF+KLB-1,1)+TEMP*X1(ISX:ISX+KLB-1)
                if (nD == 2) FLT(ISF:ISF+KLB-1,2) = FLT(ISF:ISF+KLB-1,2)+TEMP*X1(ISX:ISX+KLB-1)
              end if
              if (NFK /= 0) then
                ISF = ISTLT(IS)+LPQ
                ISD = ISTLT(KS)+1
                TEMP = DDOT_(KLB,X1(ISX),1,DLT(ISD,1),1)
                FLT(ISF,1) = FLT(ISF,1)+TEMP
                if (nD == 2) then
                  TEMP_ab = DDOT_(KLB,X1(ISX),1,DLT(ISD,2),1)
                  FLT(ISF,1) = FLT(ISF,1)+TEMP_ab
                  FLT(ISF,2) = FLT(ISF,1)
                end if
#               ifdef _DEBUGPRINT_
                write(u6,'(a,i5,a,f12.6)') '02 Flt(',isf,',1)=',FLT(ISF,1)
#               endif

              end if
            end do! JQ
          end do  ! IP
        else if ((IS == KS) .and. (JS == LS)) then
          ! CASE 3: Integrals are of symmetry type (IJ/IJ)
          ! Exchange terms need to be accumulated only
          IOPT = 1
          LPQ = 0
          IPQ = 0
          NPQ = 0
          do IP=1,IB
            do JQ=1,JB
              IPQ = IPQ+1
              LPQ = LPQ+1
              if (IPQ > NPQ) then
                call RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
                if (IRC > 1) return
                IOPT = 2
                IPQ = 1
              end if
              ISX = (IPQ-1)*KLB+1
              if (NFI /= 0) then
                ISD = ISTSQ(IS)+(IP-1)*IB+1
                ISF = ISTSQ(JS)+(JQ-1)*JB+1
                call DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                if (nD == 2) call DGEMV_('N',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
              end if
              if (NFJ /= 0) then
                ISD = ISTSQ(JS)+(JQ-1)*JB+1
                ISF = ISTSQ(IS)+(IP-1)*IB+1
                call DGEMV_('T',LB,KB,-Factor*ExFac,X1(ISX),LB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
                if (nD == 2) call DGEMV_('T',LB,KB,-factor*ExFac,X1(ISX),LB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
              end if
#             ifdef _DEBUGPRINT_
              write(u6,'(a,i5,a,f20.6)') ('03 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
              if (nD == 2) write(u6,'(a,i5,a,f20.6)') ('03 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
#             endif

            end do ! JQ
          end do   ! IP

        end if

      end do ! KS
    end do   ! JS
  end do     ! IS

end subroutine FOCKTWO_scf_Sym

subroutine FOCKTWO_scf_NoSym()

  integer(kind=iwp) :: IP, IPQ, IS, JQ, JS, KS, LS
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: ivv
# endif

  IS = 1
  IB = NBAS(IS)
  IK = KEEP(IS)
  NFI = NFRO(IS)

  JS = 1
  JB = NBAS(JS)
  JK = KEEP(JS)
  NFJ = NFRO(JS)
  IJS = Mul(IS,JS)
  IJB = nTri_Elem(IB)

  KS = 1
  KB = NBAS(KS)
  KK = KEEP(KS)
  NFK = NFRO(KS)
  LSMAX = JS

  LS = 1
  LB = NBAS(LS)
  LK = KEEP(LS)
  NFL = NFRO(LS)
  KLB = nTri_Elem(KB)

  ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?

  if ((IK+JK+KK+LK) /= 0) return
  ! NO FROZEN ORBITALS?
  if ((NFI+NFJ+NFK+NFL) == 0) return
  ! NO BASIS FUNCTIONS?
  if ((IJB*KLB) == 0) return

  ! Process the different symmetry cases

  ! CASE 1: Integrals are of symmetry type (II/II)
  ! Coulomb and exchange terms need to be accumulated
  ! Option code 1: Begin reading at first integral.
  ! NPQ: Nr of submatrices in buffer X1.
  IOPT = 1
  LPQ = 0
  IPQ = 0
  NPQ = 0
  do IP=1,IB
    do JQ=1,IP
      IPQ = IPQ+1
      LPQ = LPQ+1
      if (IPQ > NPQ) then
        call RDORD(IRC,IOPT,IS,JS,KS,LS,X1,nX1,NPQ)
        if (IRC > 1) return
        ! Option code 2: Continue reading at next integral.
        IOPT = 2
        IPQ = 1
      end if
      ISX = (IPQ-1)*KLB+1
      ISF = LPQ
      ISD = 1
      TEMP = DDOT_(KLB,X1(ISX),1,DLT(ISD,1),1)
      FLT(ISF,1) = FLT(ISF,1)+TEMP
      if (nD == 2) then
        TEMP_ab = DDOT_(KLB,X1(ISX),1,DLT(ISD,2),1)
        FLT(ISF,1) = FLT(ISF,1)+TEMP_ab
        FLT(ISF,2) = FLT(ISF,1)
      end if
#     ifdef _DEBUGPRINT_
      write(u6,'(a,i5,a,f12.6)') '00 Flt(',isf,',1)=',FLT(ISF,1)
      if (nD == 2) write(u6,'(a,i5,a,f12.6)') '00 Flt(',isf,',2)=',FLT(ISF,2)
#     endif
      call SQUARE(X1(ISX),X2,1,KB,LB)
      ISF = (JQ-1)*JB+1
      ISD = (IP-1)*IB+1

      call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
      if (nD == 2) call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
      if (IP /= JQ) then
        ISF = (IP-1)*IB+1
        ISD = (JQ-1)*JB+1

        call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,1),1,One,FSQ(ISF,1),1)
        if (nD == 2) call DGEMV_('N',KB,LB,-Factor*ExFac,X2,KB,DSQ(ISD,2),1,One,FSQ(ISF,2),1)
      end if
#     ifdef _DEBUGPRINT_
      write(u6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',1)=',FSQ(ISF+ivv-1,1),ivv=1,kb)
      if (nD == 2) write(u6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv-1,',2)=',FSQ(ISF+ivv-1,2),ivv=1,kb)
#     endif

    end do  ! JQ
  end do    ! IP

end subroutine FOCKTWO_scf_NoSym

subroutine FOCKTWO_scf_DCCD()

  use GetInt_mod, only: Basis_IDs, hash_table, I, ID_IP, lists, LuCVec, nPQ, NumCho, NumV, nVec, Vec2
  use Index_Functions, only: iTri
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: IB, IP, IP_, IPQ, IRP, IRQ, IRS, IS, ISP, ISQ, ISR, iVec1, J, JQ, JQ_, KR, KR_, LS, LS_, nData
# ifdef _DEBUGPRINT_
  integer(kind=iwp) :: ivv
# endif
  logical(kind=iwp) :: Found

  IS = 1
  IB = NBAS(IS)
  IK = KEEP(IS)
  NFI = NFRO(IS)

  call Init_GetInt(IRC)
  call Qpg_iArray('Basis IDs',Found,nData)
  call mma_allocate(Basis_IDs,4,nData/4,Label='Basis_IDs')
  call Get_iArray('Basis IDs',Basis_IDs,nData)

  call mma_allocate(hash_table,size(Basis_IDs,2),Label='hash_table')
  do IP=1,IB
    hash_table(IP) = IP
  end do
  do IP=1,IB-1
    KR = hash_table(IP)
    do JQ=IP+1,IB
      LS = hash_table(JQ)
      if (Basis_IDs(1,KR) > Basis_IDs(1,LS)) then
        IPQ = KR
        KR = LS
        LS = IPQ
        hash_table(IP) = KR
        hash_table(JQ) = LS
      end if
    end do
  end do
  IP = 1
  ID_IP = Basis_IDs(1,hash_table(IP))
  LS = 1
  do IP=2,IB
    if (ID_IP /= Basis_IDs(1,hash_table(IP))) then
      ID_IP = Basis_IDs(1,hash_table(IP))
      LS = LS+1
    end if
  end do
  call mma_allocate(lists,4,LS)
  IP = 1
  ID_IP = Basis_IDs(1,hash_table(IP))
  LS = 1
  lists(1,LS) = 1
  lists(2,LS) = ID_IP
  lists(3,LS) = IP
  lists(4,LS) = IP
  do IP=2,IB
    if (ID_IP /= Basis_IDs(1,hash_table(IP))) then
      ID_IP = Basis_IDs(1,hash_table(IP))
      LS = LS+1
      lists(1,LS) = 1
      lists(2,LS) = ID_IP
      lists(3,LS) = IP
      lists(4,LS) = IP
    else
      lists(1,LS) = lists(1,LS)+1
      lists(4,LS) = IP
    end if
  end do

  IJB = nTri_Elem(IB)
  if (nX1 < IJB) then
    write(u6,*) 'FOCKTWO_SCF_DCCD: nX1<IJB'
    call Abend()
  end if
  call Get_Int_Open(IS,IS,IS,IS)

  ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?

  if (IK /= 0) return
  ! NO FROZEN ORBITALS?
  if (NFI == 0) return
  ! NO BASIS FUNCTIONS?
  !if (IJB == 0) return

  ! Process the different symmetry cases
  do iVec1=1,NumCho(1),nVec
    NumV = min(nVec,NumCho(1)-iVec1+1)
    call RdChoVec(Vec2,nPQ,NumV,iVec1,LuCVec(1))

    ! CASE 1: Integrals are of symmetry type (II/II)
    ! Coulomb and exchange terms need to be accumulated
    ! Option code 1: Begin reading at first integral.
    do J=1,size(lists,2)
      I = J
      ID_IP = lists(2,I)

      do IP_=lists(3,I),lists(4,I)
        IP = hash_table(IP_)
        do JQ_=lists(3,I),IP_
          JQ = hash_table(JQ_)
          ! Skip processing (P,Q|... if they do not share the same center
          IPQ = iTri(IP,JQ)
          ! do batches of integrals for a single fixed pair of pq
          call Get_Int_DCCD(IRC,X1,IPQ,IJB+1)
          if (IRC > 1) return
          ! Do the Coulomb contribution
          if (nD == 1) then
            TEMP = Zero
            do KR_=lists(3,I),lists(4,I)
              KR = hash_table(KR_)
              do LS_=lists(3,I),KR_
                LS = hash_table(LS_)
                IRS = iTri(KR,LS)
                TEMP = TEMP+X1(IRS)*DLT(IRS,1)
              end do
            end do
            FLT(IPQ,1) = FLT(IPQ,1)+TEMP
          else
            TEMP = Zero
            TEMP_ab = Zero
            do KR_=lists(3,I),lists(4,I)
              KR = hash_table(KR_)
              do LS_=lists(3,I),KR_
                LS = hash_table(LS_)
                IRS = iTri(KR,LS)
                TEMP = TEMP+X1(IRS)*DLT(IRS,1)
                TEMP_ab = TEMP_ab+X1(IRS)*DLT(IRS,2)
              end do
            end do
            FLT(IPQ,1) = FLT(IPQ,1)+TEMP+TEMP_ab
            FLT(IPQ,2) = FLT(IPQ,1)
          end if
#         ifdef _DEBUGPRINT_
          write(u6,'(a,i5,a,f12.6)') '00 Flt(',IPQ,',1)=',FLT(IPQ,1)
          if (nD == 2) write(u6,'(a,i5,a,f12.6)') '00 Flt(',IPQ,',2)=',FLT(IPQ,2)
#         endif
          ! Do the exchange contribution
          call SQUARE(X1,X2,1,IB,IB)

          if (nD == 1) then
            do KR_=lists(3,I),lists(4,I)
              KR = hash_table(KR_)
              IRQ = (JQ-1)*IB+KR
              TEMP = Zero
              do LS_=lists(3,I),lists(4,I)
                LS = hash_table(LS_)
                ISR = (KR-1)*IB+LS
                ISP = (IP-1)*IB+LS
                TEMP = TEMP-Factor*ExFac*X2(ISR)*DSQ(ISP,1)
              end do
              FSQ(IRQ,1) = FSQ(IRQ,1)+TEMP
            end do
          else
            do KR_=lists(3,I),lists(4,I)
              KR = hash_table(KR_)
              IRQ = (JQ-1)*IB+KR
              TEMP = Zero
              TEMP_ab = Zero
              do LS_=lists(3,I),lists(4,I)
                LS = hash_table(LS_)
                ISR = (KR-1)*IB+LS
                ISP = (IP-1)*IB+LS
                TEMP = TEMP-Factor*ExFac*X2(ISR)*DSQ(ISP,1)
                TEMP_ab = TEMP_ab-Factor*ExFac*X2(ISR)*DSQ(ISP,2)
              end do
              FSQ(IRQ,1) = FSQ(IRQ,1)+TEMP
              FSQ(IRQ,2) = FSQ(IRQ,2)+TEMP_ab
            end do
          end if
          if (IP /= JQ) then
            ISF = (IP-1)*IB+1
            ISD = (JQ-1)*IB+1

            if (nD == 1) then
              do KR_=lists(3,I),lists(4,I)
                KR = hash_table(KR_)
                IRP = (IP-1)*IB+KR
                TEMP = Zero
                do LS_=lists(3,I),lists(4,I)
                  LS = hash_table(LS_)
                  ISR = (KR-1)*IB+LS
                  ISQ = (JQ-1)*IB+LS
                  TEMP = TEMP-Factor*ExFac*X2(ISR)*DSQ(ISQ,1)
                end do
                FSQ(IRP,1) = FSQ(IRP,1)+TEMP
              end do
            else
              do KR_=lists(3,I),lists(4,I)
                KR = hash_table(KR_)
                IRP = (IP-1)*IB+KR
                TEMP = Zero
                TEMP_ab = Zero
                do LS_=lists(3,I),lists(4,I)
                  LS = hash_table(LS_)
                  ISR = (KR-1)*IB+LS
                  ISQ = (JQ-1)*IB+LS
                  TEMP = TEMP-Factor*ExFac*X2(ISR)*DSQ(ISQ,1)
                  TEMP_ab = TEMP_ab-Factor*ExFac*X2(ISR)*DSQ(ISQ,2)
                end do
                FSQ(IRP,1) = FSQ(IRP,1)+TEMP
                FSQ(IRP,2) = FSQ(IRP,2)+TEMP_ab
              end do
            end if

          end if
#         ifdef _DEBUGPRINT_
          ISF = (JQ-1)*IB
          write(u6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv,',1)=',FSQ(ISF+ivv,1),ivv=1,Ib)
          if (nD == 2) write(u6,'(a,i5,a,f12.6)') ('01 Fsq(',isf+ivv,',2)=',FSQ(ISF+ivv,2),ivv=1,Ib)
#         endif

        end do  ! JQ
      end do    ! IP
    end do      ! I
  end do  ! end of the batch procedure

  call mma_deallocate(lists)
  call Get_Int_Close()
  call mma_deallocate(hash_table)
  call mma_deallocate(Basis_IDs)

end subroutine FOCKTWO_scf_DCCD

end subroutine FOCKTWO_scf
