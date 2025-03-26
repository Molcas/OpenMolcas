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

implicit real*8(A-H,O-Z)
! Input
dimension RHO2T(*)
! Input and output
dimension RHO2(*)

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
SIGN = 0.0d0
IACTIVE = 0

#ifdef _DEBUGPRINT_
write(6,*) ' Welcome to ADTOR2'
write(6,*) ' ================='
write(6,*) ' NI NJ NK NL = ',NI,NJ,NK,NL
write(6,*) ' IOFF JOFF KOFF LOFF =',IOFF,JOFF,KOFF,LOFF
write(6,*) ' ITYPE = ',ITYPE
write(6,*) ' Initial two body density matrix'
call PRSYM(RHO2,NORB**2)
!write(6,*) ' RHO2T :'
!if (ITYPE == 1) then
!  if (IOFF == KOFF) then
!    NROW = NI*(NI+1)/2
!  else
!    NROW = NI*NK
!  end if
!  if (JOFF == LOFF) then
!    NCOL = NJ*(NJ+1)/2
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
      SIGN = 1.0d0
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
      SIGN = -1.0d0
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
        SIGN = -1.0d0
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
        SIGN = 1.0d0
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
                IJKL = IJ*(IJ-1)/2+KL
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
                  SIGNIK = 1.0d0
                else
                  IKIND = max(I,K)*(max(I,K)-1)/2+min(I,K)
                  NIK = NI*(NI+1)/2
                  if (I == max(I,K)) then
                    SIGNIK = 1.0d0
                  else
                    SIGNIK = -1.0d0
                  end if
                end if
                if (JOFF /= LOFF) then
                  JLIND = (L-1)*NJ+J
                  SIGNJL = 1.0d0
                else
                  JLIND = max(J,L)*(max(J,L)-1)/2+min(J,L)
                  if (J == max(J,L)) then
                    SIGNJL = 1.0d0
                  else
                    SIGNJL = -1.0d0
                  end if
                end if
                IKJLT = (JLIND-1)*NIK+IKIND
                RHO2(IJKL) = RHO2(IJKL)-SIGN*SIGNJL*SIGNIK*RHO2T(IKJLT)
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
            FACTOR = 2.0d0
          else
            FACTOR = 1.0d0
          end if
          IJKL = max(IJ,KL)*(max(IJ,KL)-1)/2+min(IJ,KL)
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
          !IJKL = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
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
!iA1 = itri((i-1)*nOrb+j,(k-1)*nOrb+l)
!iA2 = itri((j-1)*nOrb+i,(k-1)*nOrb+l)
!iA3 = itri((i-1)*nOrb+j,(l-1)*nOrb+k)
!iA4 = itri((j-1)*nOrb+i,(l-1)*nOrb+k)
!rfel = RHO2(ia1)+RHO2(ia2)+RHO2(ia3)+RHO2(ia4)
!if (rfel /= 0.0D0) write(6,*) 'STOP'
!#ifdef _DEBUGPRINT_
!write(6,*) ' Updated two-body density matrix'
!call PRSYM(RHO2,NORB**2)
!#endif

return

end subroutine AdToR2_MCLR
