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
! Copyright (C) 1998, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,LADVICE)
! Advice Sigma routine about best route to take
!
! LADVICE : ADVICE given (short, an integer !!)
!
! For ITERM = 1 :
!   LADVICE = 1 : Business as usual, no transpose of matrix
!                 (resolution on alpha strings, direct exc on beta)
!   LADVICE = 2 = Transpose matrices
!                 (resolution on beta strings, direct exc on alpha)
! (SVC: one call to this routine and ITERM is one, so I removed that
! argument and skipped the checking. Also, the arguments are all scalar,
! so that has been hard-coded now too)
!
! Jeppe Olsen, Tirstrup Airport, Jan 12, 98

use lucia_data, only: IADVICE, IPHGAS, MNHL, NGAS, NOBPT
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IAOCC(*), IBOCC(*), JAOCC(*), JBOCC(*)
integer(kind=iwp), intent(out) :: LADVICE
integer(kind=iwp) :: IGAS, IPHMODI, ITP(16), JTP(16), KHOLEA, KHOLEB, KTP(16), LHOLEA, LHOLEB, LLADVICE, LTP(16), NIJTYP, NKLTYP
real(kind=wp) :: XCJKAJB, XCLJAKB, XFLOPA, XFLOPB, XNIJSX, XNIOB, XNJEL, XNJOB, XNKLSX, XNKOB, XNLEL, XNLOB

! sigma(i,Ka,Ib) = sum(i,kl)<Ib!Eb_kl!Jb>(ij!kl)C(j,Ka,Jb)
!
! Number of ops : Number of sx(kl) N_i*N_j_dimension of C(j,Ka,Jb)
! No absolute calc of flops is made, only a relative measure
!
! Single excitations connecting the two types

!    SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
call SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOCC,JAOCC,IPHGAS)
call SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOCC,JBOCC,IPHGAS)
!write(u6,*) 'NIJTYP, NKLTYP',NIJTYP,NKLTYP
! P-h modifications (I cannot predict these at the moment
if ((NIJTYP >= 1) .and. (NKLTYP >= 1)) then

  if (((IPHGAS(ITP(1)) == 2) .and. (IPHGAS(JTP(1)) == 2)) .or. ((IPHGAS(KTP(1)) == 2) .and. (IPHGAS(LTP(1)) == 2))) then
    IPHMODI = 1
  else
    IPHMODI = 0
  end if
else
  IPHMODI = 0
end if

if ((IPHMODI == 1) .or. (NIJTYP /= 1) .or. (NKLTYP /= 1) .or. (IADVICE == 0)) then
  ! Several connections, i.e. the alpha or the beta blocks are identical,
  ! or ph modifications
  ! just continue
  LADVICE = 1
else
  ! ======================================
  ! Index for flops along C(j,Ka,Jb) route
  ! ======================================
  ! Dim of C(j,Ka,Jb) relative to C(Ja,Jb)
  ! going from Ja to  Ka reduces occ by one elec, changes dim by n/(N-n+1)
  XNJOB = real(NOBPT(JTP(1)),kind=wp)
  XNJEL = real(JAOCC(JTP(1)),kind=wp)
  XCJKAJB = XNJOB*XNJEL/(XNJOB-XNJEL+One)
  ! Number of kl excitations per beta string :
  XNKLSX = real((NOBPT(KTP(1))-JBOCC(KTP(1)))*JBOCC(LTP(1)),kind=wp)
  ! Number of ops (relative to dim of C)
  XNIOB = real(NOBPT(ITP(1)),kind=wp)
  XFLOPA = XCJKAJB*XNKLSX*XNIOB
  ! ======================================
  ! Index for flops along C(l,Ja,Kb) route
  ! ======================================
  ! Dim of C(l,Ja,Kb) relative to C(Ja,Jb)
  XNLOB = real(NOBPT(LTP(1)),kind=wp)
  XNLEL = real(JBOCC(LTP(1)),kind=wp)
  XCLJAKB = XNLOB*XNLEL/(XNLOB-XNLEL+One)
  ! Number of ij excitations per alpha string :
  XNIJSX = real((NOBPT(ITP(1))-JAOCC(ITP(1)))*JAOCC(JTP(1)),kind=wp)
  ! Number of ops (relative to dim of C)
  XNKOB = real(NOBPT(KTP(1)),kind=wp)
  XFLOPB = XCLJAKB*XNIJSX*XNKOB
  ! Switch to second route if atleast 20 percent less work
  if (XFLOPB <= 0.8_wp*XFLOPA) then
    LADVICE = 2
  else
    LADVICE = 1
  end if
  ! Well, an additional consideration :
  ! If the C block involes the smallest allowed number of elecs in hole space,
  ! and the annihilation is in hole space
  ! then we do the annihilation in the space with the smallest number of
  ! hole electrons.
  LHOLEA = 0
  LHOLEB = 0
  do IGAS=1,NGAS
    if (IPHGAS(IGAS) == 2) then
      LHOLEA = LHOLEA+JAOCC(IGAS)
      LHOLEB = LHOLEB+JBOCC(IGAS)
    end if
  end do

  if ((LHOLEA+LHOLEB == MNHL) .and. ((IPHGAS(JTP(1)) == 2) .or. (IPHGAS(LTP(1)) == 2))) then

    if (IPHGAS(JTP(1)) == 2) then
      KHOLEA = LHOLEA-1
      KHOLEB = LHOLEB
    else
      KHOLEA = LHOLEA
      KHOLEB = LHOLEB-1
    end if

    if (KHOLEA == KHOLEB) then
      LLADVICE = LADVICE
    else if (KHOLEA < KHOLEB) then
      LLADVICE = 1
    else
      LLADVICE = 2
    end if
#   ifdef _DEBUGPRINT_
    if (LADVICE /= LLADVICE) then
      write(u6,*) ' Advice changed by hole considetions'
      write(u6,*) ' LADVICE, LLADVICE',LADVICE,LLADVICE
    end if
#   endif
    LADVICE = LLADVICE
  end if

# ifdef _DEBUGPRINT_
  write(u6,*) ' ADVICE active'
  write(u6,*) ' IAOCC JAOCC IBOCC JBOCC'
  call IWRTMA(IAOCC,1,NGAS,1,NGAS)
  call IWRTMA(JAOCC,1,NGAS,1,NGAS)
  call IWRTMA(IBOCC,1,NGAS,1,NGAS)
  call IWRTMA(JBOCC,1,NGAS,1,NGAS)
  write(u6,*) ' ITP JTP KTP LTP ',ITP(1),JTP(1),KTP(1),LTP(1)
  write(u6,*) ' XFLOPA,XFLOPB',XFLOPA,XFLOPB
  write(u6,*) ' ADVICE given : ',LADVICE
# endif
end if
! End if several types/ph modi

end subroutine ADVICE_SIGMA
