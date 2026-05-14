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
!               2003, Jesper Wisborg Krogh                             *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ADTOR2(RHO2,RHO2S,RHO2A,RHO2T,ITYPE,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB,IPACK)
! Add contributions to two electron density matrix RHO2
! output density matrix is in the form Rho2(ij,kl),(ij) >= (kl)
!
! Jeppe Olsen, Fall of 96
! Jesper Wisborg Krogh, September 03: Can now symmetry pack on the fly
!
! Itype = 1 => alpha-alpha or beta-beta loop
!              input is in form Rho2t(ik,jl)
! Itype = 2 => alpha-beta loop
!              input is in form Rho2t(ij,kl)

use Index_Functions, only: iTri, nTri_Elem
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: RHO2(*), RHO2S(*), RHO2A(*)
real(kind=wp), intent(in) :: RHO2T(*)
integer(kind=iwp), intent(in) :: ITYPE, NI, IOFF, NJ, JOFF, NK, KOFF, NL, LOFF, NORB
logical(kind=iwp), intent(in) :: IPACK
integer(kind=iwp) :: I, I_PACK, IACTIVE, II, IIOFF, IJ, IJ_PACK, IJKL, IJKL_PACK, IJKLT, IKIND, IKJLT, IPERM, J, J_PACK, JI_PACK, &
                     JIKL_PACK, JJ, JJOFF, JLIND, K, K_PACK, KK, KKOFF, KL, KL_PACK, L, L_PACK, LL, LLOFF, NELMNT, NII, NIK, NJJ, &
                     NKK, NLL
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NCOL, NROW
#endif
real(kind=wp) :: FACTOR, FACTOR_PACK, SGN, SIGNIK, SIGNJL, TERM
logical(kind=iwp) :: DO_IJKL, DO_IJKL_PACK, DO_JIKL_PACK, IPROCEED

! Some dummy initializations
NII = 0 ! jwk-cleanup
NJJ = 0 ! jwk-cleanup
NKK = 0 ! jwk-cleanup
NLL = 0 ! jwk-cleanup
KKOFF = 0 ! jwk-cleanup
LLOFF = 0 ! jwk-cleanup
SGN = One ! jwk-cleanup
IACTIVE = 1 ! jwk-cleanup
I = 0 ! jwk-cleanup
K = 0 ! jwk-cleanup
J = 0 ! jwk-cleanup
L = 0 ! jwk-cleanup
I_PACK = 0
J_PACK = 0
K_PACK = 0
L_PACK = 0
IJ_PACK = 0
JI_PACK = 0
KL_PACK = 0
FACTOR = Zero

#ifdef _DEBUGPRINT_
write(u6,*) ' Welcome to ADTOR2'
write(u6,*) ' ================='
write(u6,*) ' NI NJ NK NL = ',NI,NJ,NK,NL
write(u6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
write(u6,*) ' ITYPE = ',ITYPE
write(u6,*) ' NORB ',NORB
if (.not. IPACK) then
  write(u6,*) ' Initial two body density matrix'
  call PRSYM(RHO2,NORB**2)
end if

write(u6,*) ' RHO2T :'
if (ITYPE == 1) then
  if (IOFF == KOFF) then
    NROW = nTri_Elem(NI)
  else
    NROW = NI*NK
  end if
  if (JOFF == LOFF) then
    NCOL = nTri_Elem(NJ)
  else
    NCOL = NJ*NL
  end if
else if (ITYPE == 2) then
  NROW = NI*NJ
  NCOL = NK*NL
end if
call WRTMAT(RHO2T,NROW,NCOL,NROW,NCOL)
#endif

!return

NELMNT = nTri_Elem(NORB**2)

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
      !IJOFF = (JJOFF-1)*NORB+IIOFF
      !KLOFF = (LLOFF-1)*NORB+LLOFF
      do II=1,NII
        do JJ=1,NJJ
          do KK=1,NKK
            do LL=1,NLL
              IJ = (JJ+JJOFF-2)*NORB+II+IIOFF-1
              KL = (LL+LLOFF-2)*NORB+KK+KKOFF-1

              IPROCEED = .false.
              DO_IJKL = .false.
              DO_IJKL_PACK = .false.
              DO_JIKL_PACK = .false.
              if (IPACK) then
                I_PACK = II+IIOFF-1
                J_PACK = JJ+JJOFF-1
                K_PACK = KK+KKOFF-1
                L_PACK = LL+LLOFF-1

                IJ_PACK = nTri_Elem(I_PACK-1)+J_PACK
                JI_PACK = nTri_Elem(J_PACK-1)+I_PACK
                KL_PACK = nTri_Elem(K_PACK-1)+L_PACK

                if (K_PACK == L_PACK) then
                  FACTOR = Quart
                else
                  FACTOR = Half
                end if

                if (K_PACK >= L_PACK) then
                  if ((I_PACK >= J_PACK) .and. (IJ_PACK >= KL_PACK)) then
                    DO_IJKL_PACK = .true.
                    IPROCEED = .true.
                  end if
                  if ((J_PACK >= I_PACK) .and. (JI_PACK >= KL_PACK)) then
                    DO_JIKL_PACK = .true.
                    IPROCEED = .true.
                  end if
                end if
              else if (IJ >= KL) then
                DO_IJKL = .true.
                IPROCEED = .true.
                FACTOR = One
              end if

              if (IPROCEED) then
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
                TERM = FACTOR*SGN*SIGNJL*SIGNIK*RHO2T(IKJLT)

                if (DO_IJKL_PACK) then
                  IJKL_PACK = nTri_Elem(IJ_PACK-1)+KL_PACK
                  RHO2S(IJKL_PACK) = RHO2S(IJKL_PACK)-TERM
                  RHO2A(IJKL_PACK) = RHO2A(IJKL_PACK)-TERM
                end if
                if (DO_JIKL_PACK) then
                  JIKL_PACK = nTri_Elem(JI_PACK-1)+KL_PACK
                  RHO2S(JIKL_PACK) = RHO2S(JIKL_PACK)-TERM
                  RHO2A(JIKL_PACK) = RHO2A(JIKL_PACK)+TERM
                end if
                if (DO_IJKL) then
                  IJKL = nTri_Elem(IJ-1)+KL
                  if (IJKL > NELMNT) then
                    write(u6,*) ' Problemo 1 : IJKL > NELMNT'
                    write(u6,*) ' IJKL, NELMNT',IJKL,NELMNT
                    write(u6,*) ' IJ, KL',IJ,KL
                    write(u6,*) ' JJ JJOFF ',JJ,JJOFF
                    write(u6,*) ' II IIOFF ',II,IIOFF
                    write(u6,*) ' IPERM = ',IPERM
                  end if
                  RHO2(IJKL) = RHO2(IJKL)-TERM
                end if
                ! The minus : Rho2t comes as <a+i a+k aj al>, but we want
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
          if (IJKL > NELMNT) then
            write(u6,*) ' Problemo 2 : IJKL > NELMNT'
            write(u6,*) ' IJKL, NELMNT',IJKL,NELMNT
          end if
          TERM = FACTOR*RHO2T(IJKLT)
          if (IPACK) then
            do IPERM=1,2
              if (IPERM == 1) then
                I_PACK = I+IOFF-1
                J_PACK = J+JOFF-1
                K_PACK = K+KOFF-1
                L_PACK = L+LOFF-1
                IJ_PACK = nTri_Elem(I_PACK-1)+J_PACK
                JI_PACK = nTri_Elem(J_PACK-1)+I_PACK
                KL_PACK = nTri_Elem(K_PACK-1)+L_PACK
              else
                I_PACK = K+KOFF-1
                J_PACK = L+LOFF-1
                K_PACK = I+IOFF-1
                L_PACK = J+JOFF-1
                IJ_PACK = nTri_Elem(I_PACK-1)+J_PACK
                JI_PACK = nTri_Elem(J_PACK-1)+I_PACK
                KL_PACK = nTri_Elem(K_PACK-1)+L_PACK
              end if

              if (K_PACK == L_PACK) then
                FACTOR_PACK = Quart
              else
                FACTOR_PACK = Half
              end if
              if ((I_PACK == K_PACK) .and. (J_PACK == L_PACK)) FACTOR_PACK = FACTOR_PACK*Half

              if ((I_PACK >= J_PACK) .and. (K_PACK >= L_PACK) .and. (IJ_PACK >= KL_PACK)) then
                IJKL_PACK = nTri_Elem(IJ_PACK-1)+KL_PACK
                RHO2S(IJKL_PACK) = RHO2S(IJKL_PACK)+FACTOR_PACK*TERM
                RHO2A(IJKL_PACK) = RHO2A(IJKL_PACK)+FACTOR_PACK*TERM
              end if
              if ((J_PACK >= I_PACK) .and. (K_PACK >= L_PACK) .and. (JI_PACK >= KL_PACK)) then
                JIKL_PACK = nTri_Elem(JI_PACK-1)+KL_PACK
                IJKL_PACK = nTri_Elem(IJ_PACK-1)+KL_PACK
                RHO2S(JIKL_PACK) = RHO2S(JIKL_PACK)+FACTOR_PACK*TERM
                RHO2A(JIKL_PACK) = RHO2A(JIKL_PACK)-FACTOR_PACK*TERM
              end if
            end do
          else
            RHO2(IJKL) = RHO2(IJKL)+TERM
          end if
        end do
      end do
    end do
  end do

end if

end subroutine ADTOR2
