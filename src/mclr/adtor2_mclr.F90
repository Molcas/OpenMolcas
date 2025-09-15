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
! Copyright (C) 1996, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine AdToR2_MCLR(RHO2,RHO2T,ITYPE,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
! Add contributions to two electron density matrix RHO2
! output density matrix is in the form Rho2(ij,kl),(ij) >= (kl)
!
! Jeppe Olsen, Fall of 96
!
! Itype = 1 => alpha-alpha or beta-beta loop
!              input is in form Rho2t(ik,jl)
! Itype = 2 => alpha-beta loop
!              input is in form Rho2t(ij,kl)

use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(inout) :: RHO2(*)
real(kind=wp), intent(in) :: RHO2T(*)
integer(kind=iwp), intent(in) :: ITYPE, NI, IOFF, NJ, JOFF, NK, KOFF, NL, LOFF, NORB
integer(kind=iwp) :: I, IACTIVE, II, IIOFF, IJ, IJKL, IJKLT, IKIND, IKJLT, IPERM, J, JJ, JJOFF, JLIND, K, KK, KKOFF, KL, L, LL, &
                     LLOFF, NII, NIK, NJJ, NKK, NLL
real(kind=wp) :: FACTOR, SGN, SIGNIK, SIGNJL

! dummy initialize
I = 0
J = 0
K = 0
L = 0
NII = 0
NJJ = 0
NKK = 0
NLL = 0
IIOFF = 0
JJOFF = 0
KKOFF = 0
LLOFF = 0
SGN = Zero
IACTIVE = 0

#ifdef _DEBUGPRINT_
write(u6,*) ' Welcome to ADTOR2'
write(u6,*) ' ================='
write(u6,*) ' NI NJ NK NL = ',NI,NJ,NK,NL
write(u6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
write(u6,*) ' ITYPE = ',ITYPE
write(u6,*) ' Initial two body density matrix'
call PRSYM(RHO2,NORB**2)
!write(u6,*) ' RHO2T :'
!if (ITYPE == 1) then
!  if (IOFF == KOFF) then
!    NROW = nTri_Elem(NI)
!  else
!    NROW = NI*NK
!  end if
!  if (JOFF == LOFF) then
!    NCOL = nTri_Elem(NJ)
!  else
!    NCOL = NJ*NL
!  end if
!else if (ITYPE == 2) then
!  NROW = NI*NJ
!  NCOL = NK*NL
!end if
!call WRTMAT(RHO2T,NROW,NCOL,NROW,NCOL)
#endif

if (ITYPE == 1) then

  ! =======================================
  !     Alpha-alpha or beta-beta term
  ! =======================================

  ! Four permutations
  do IPERM=1,4
    if (IPERM == 1) then
      NII = NI
      IIOFF = IOFF
      NJJ = NJ
      JJOFF = JOFF
      NKK = NK
      KKOFF = KOFF
      NLL = NL
      LLOFF = LOFF
      SGN = One
      IACTIVE = 1
    else if (IPERM == 2) then
      if (IOFF /= KOFF) then
        NII = NK
        IIOFF = KOFF
        NKK = NI
        KKOFF = IOFF
        NJJ = NJ
        JJOFF = JOFF
        NLL = NL
        LLOFF = LOFF
        IACTIVE = 1
      else
        IACTIVE = 0
      end if
      SGN = -One
    else if (IPERM == 3) then
      if (JOFF /= LOFF) then
        NII = NI
        IIOFF = IOFF
        NKK = NK
        KKOFF = KOFF
        NJJ = NL
        JJOFF = LOFF
        NLL = NJ
        LLOFF = JOFF
        SGN = -One
        IACTIVE = 1
      else
        IACTIVE = 0
      end if
    else if (IPERM == 4) then
      if ((IOFF /= KOFF) .and. (JOFF /= LOFF)) then
        NKK = NI
        KKOFF = IOFF
        NII = NK
        IIOFF = KOFF
        NJJ = NL
        JJOFF = LOFF
        NLL = NJ
        LLOFF = JOFF
        SGN = One
        IACTIVE = 1
      else
        IACTIVE = 0
      end if
    end if

    !IJOFF = (JJOFF-1)*NORB+IIOFF
    !KLOFF = (LLOFF-1)*NORB+KKOFF
    !if ((IACTIVE == 1) .and. (IJOFF >= KLOFF)) then
    if (IACTIVE == 1) then
      do II=1,NII
        do JJ=1,NJJ
          do KK=1,NKK
            do LL=1,NLL
              IJ = (JJ+JJOFF-2)*NORB+II+IIOFF-1
              KL = (LL+LLOFF-2)*NORB+KK+KKOFF-1
              if (IJ >= KL) then
                IJKL = nTri_Elem(IJ-1)+KL
                if (IPERM == 1) then
                  I = II
                  K = KK
                  J = JJ
                  L = LL
                else if (IPERM == 2) then
                  I = KK
                  K = II
                  J = JJ
                  L = LL
                else if (IPERM == 3) then
                  I = II
                  K = KK
                  J = LL
                  L = JJ
                else if (IPERM == 4) then
                  I = KK
                  K = II
                  J = LL
                  L = JJ
                end if
                if (IOFF /= KOFF) then
                  IKIND = (K-1)*NI+I
                  NIK = NI*NK
                  SIGNIK = One
                else
                  IKIND = iTri(I,K)
                  NIK = nTri_Elem(NI)
                  if (I >= K) then
                    SIGNIK = One
                  else
                    SIGNIK = -One
                  end if
                end if
                if (JOFF /= LOFF) then
                  JLIND = (L-1)*NJ+J
                  SIGNJL = One
                else
                  JLIND = iTri(J,L)
                  if (J >= L) then
                    SIGNJL = One
                  else
                    SIGNJL = -One
                  end if
                end if
                IKJLT = (JLIND-1)*NIK+IKIND
                RHO2(IJKL) = RHO2(IJKL)-SGN*SIGNJL*SIGNIK*RHO2T(IKJLT)
                ! The minus: Rho2t comes as <a+i a+k aj al>, but we want
                ! <a+ia+k al aj>
              end if
            end do
          end do
        end do
      end do
      ! End of active/inactive if
    end if
    ! End of loop over permutations
  end do
else if (ITYPE == 2) then

  ! =======================================
  !     Alpha-alpha or beta-beta term
  ! =======================================

  do I=1,NI
    do J=1,NJ
      do K=1,NK
        do L=1,NL
          IJ = (J+JOFF-2)*NORB+I+IOFF-1
          KL = (L+LOFF-2)*NORB+K+KOFF-1
          if (IJ == KL) then
            FACTOR = Two
          else
            FACTOR = One
          end if
          IJKL = iTri(IJ,KL)
          IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI+(J-1)*NI+I
          RHO2(IJKL) = RHO2(IJKL)+FACTOR*RHO2T(IJKLT)
        end do
      end do
    end do
  end do

else if (itype == 3) then
  do I=1,NI
    do J=1,NJ
      do K=1,NK
        do L=1,NL
          IJ = (J+JOFF-2)*NORB+I+IOFF-1
          KL = (L+LOFF-2)*NORB+K+KOFF-1
          !IJKL = iTri(IJ,KL)
          IJKL = (IJ-1)*NORB**2+KL
          IJKLT = (L-1)*NJ*NK*NI+(K-1)*NJ*NI+(J-1)*NI+I
          RHO2(IJKL) = RHO2(IJKL)+RHO2T(IJKLT)
        end do
      end do
    end do
  end do

end if

!i = 1
!j = 4
!k = 5
!l = 6
!iA1 = iTri((i-1)*nOrb+j,(k-1)*nOrb+l)
!iA2 = iTri((j-1)*nOrb+i,(k-1)*nOrb+l)
!iA3 = iTri((i-1)*nOrb+j,(l-1)*nOrb+k)
!iA4 = iTri((j-1)*nOrb+i,(l-1)*nOrb+k)
!rfel = RHO2(ia1)+RHO2(ia2)+RHO2(ia3)+RHO2(ia4)
!if (rfel /= Zero) write(u6,*) 'STOP'
!#ifdef _DEBUGPRINT_
!write(u6,*) ' Updated two-body density matrix'
!call PRSYM(RHO2,NORB**2)
!#endif

return

end subroutine AdToR2_MCLR
