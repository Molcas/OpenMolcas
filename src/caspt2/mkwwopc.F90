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

subroutine MKWWOPC(IVEC,JVEC,OP1,NOP2,OP2,NOP3,OP3)
! Presently symmetry blocking is disregarded, but index pair
! permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
! NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)

! Given the coefficients for two excitation operators of the
! type ATVX = Case C, available in vectors numbered IVEC and
! JVEC on file, construct the zero-, one-, two-, and three-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.
! Formula used:
!  W1(tuv,a)(conj)*W2(xyz,b) = dab * ( Evutxyz +dyu Evztx
!                       + dyx Evutz + dtu Evxyz + dtu dyx Evz )

use SUPERINDEX, only: MTUV
use EQSOLV, only: MODVEC
use caspt2_module, only: NASHT, NASUP, NINDEP, NISUP, NSYM, NTUVES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, NOP2, NOP3
real(kind=wp), intent(inout) :: OP1(NASHT,NASHT), OP2(NOP2), OP3(NOP3)
integer(kind=iwp) :: ICASE, IIEND, IISTA, ISCT, ISYM, ITABS, ITUV, ITUVABS, ITUVEND, ITUVSTA, ITX, ITZ, IUABS, IVABS, IVU, IVX, &
                     IVZ, IW1, IW2, IWPROD, IXABS, IXYZ, IXYZABS, IXYZEND, IXYZSTA, IYABS, IYZ, IZABS, JTX, JVU, JVUTXYZ, JVUTZ, &
                     JVXYZ, JVZTX, JYZ, LW1A, LW2A, MDVEC, MWS1, MWS2, NAS, NCOL, NIS, NWPROD, NWSCT
real(kind=wp) :: W_PROD
real(kind=wp), allocatable :: W1(:), W2(:), WPROD(:)

ICASE = 4
! Loop over symmetry ISYM
do ISYM=1,NSYM
  NAS = NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE)
  if (NINDEP(ISYM,ICASE) == 0) cycle
  ! Allocate space for one section of excitation amplitudes:
  MDVEC = MODVEC(ISYM,ICASE)
  call mma_allocate(W1,NAS*MDVEC,Label='W1')
  call mma_allocate(W2,NAS*MDVEC,Label='W2')
  NWSCT = min(NAS,1000)
  NWPROD = NWSCT**2
  ! Allocate space for the contraction:
  call mma_allocate(WPROD,NWPROD,Label='WPROD')
  ! Sectioning loop added:
  ISCT = 0
  do IISTA=1,NIS,MDVEC
    ISCT = ISCT+1
    IIEND = min(IISTA-1+MDVEC,NIS)
    NCOL = 1+IIEND-IISTA
    call RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
    call RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
    ! Loop over sections of WW1 and WW2:
    do ITUVSTA=1,NAS,NWSCT
      LW1A = ITUVSTA
      ITUVEND = min(ITUVSTA-1+NWSCT,NAS)
      MWS1 = ITUVEND+1-ITUVSTA
      do IXYZSTA=1,NAS,NWSCT
        IXYZEND = min(IXYZSTA-1+NWSCT,NAS)
        LW2A = IXYZSTA
        MWS2 = IXYZEND+1-IXYZSTA
        ! Multiply WProd = (W1 sect )*(W2 sect transpose)
        WPROD(:) = Zero
        call DGEMM_('N','T',MWS1,MWS2,NCOL,One,W1(LW1A),NAS,W2(LW2A),NAS,One,WPROD,NWSCT)

        ! Loop over (TUV) in its section
        do ITUV=ITUVSTA,ITUVEND
          IW1 = ITUV+1-ITUVSTA
          ITUVABS = ITUV+NTUVES(ISYM)
          ITABS = MTUV(1,ITUVABS)
          IUABS = MTUV(2,ITUVABS)
          IVABS = MTUV(3,ITUVABS)
          IVU = IVABS+NASHT*(IUABS-1)
          ! Loop over (XYZ) in its section
          do IXYZ=IXYZSTA,IXYZEND
            IW2 = IXYZ+1-IXYZSTA
            IXYZABS = IXYZ+NTUVES(ISYM)
            IXABS = MTUV(1,IXYZABS)
            IYABS = MTUV(2,IXYZABS)
            IZABS = MTUV(3,IXYZABS)
            ITX = ITABS+NASHT*(IXABS-1)
            IYZ = IYABS+NASHT*(IZABS-1)
            IWPROD = IW1+NWSCT*(IW2-1)
            W_PROD = WPROD(IWPROD)
            ! Remember:
            !  W1(tuv,a)(conj)*W2(xyz,b) = dab * ( Evutxyz + dyu Evztx + dyx Evutz + dtu Evxyz + dtu dyx Evz )
            ! Contrib to 3-particle operator:
            if (IVU < ITX) then
              if (IVU >= IYZ) then
                JVU = ITX
                JTX = IVU
                JYZ = IYZ
              else if (ITX < IYZ) then
                JVU = IYZ
                JTX = ITX
                JYZ = IVU
              else
                JVU = ITX
                JTX = IYZ
                JYZ = IVU
              end if
            else
              if (IVU < IYZ) then
                JVU = IYZ
                JTX = IVU
                JYZ = ITX
              else if (ITX >= IYZ) then
                JVU = IVU
                JTX = ITX
                JYZ = IYZ
              else
                JVU = IVU
                JTX = IYZ
                JYZ = ITX
              end if
            end if
            JVUTXYZ = ((JVU+1)*JVU*(JVU-1))/6+(JTX*(JTX-1))/2+JYZ
            OP3(JVUTXYZ) = OP3(JVUTXYZ)+WPROD(IWPROD)
            ! Contrib to 2-particle operator, from  dyu Evztx:
            if (IYABS == IUABS) then
              IVZ = IVABS+NASHT*(IZABS-1)
              ITX = ITABS+NASHT*(IXABS-1)
              if (IVZ >= ITX) then
                JVZTX = (IVZ*(IVZ-1))/2+ITX
              else
                JVZTX = (ITX*(ITX-1))/2+IVZ
              end if
              OP2(JVZTX) = OP2(JVZTX)+W_PROD
            end if
            ! Contrib to 2-particle operator, from  dyx Evutz:
            if (IYABS == IXABS) then
              IVU = IVABS+NASHT*(IUABS-1)
              ITZ = ITABS+NASHT*(IZABS-1)
              if (IVU >= ITZ) then
                JVUTZ = (IVU*(IVU-1))/2+ITZ
              else
                JVUTZ = (ITZ*(ITZ-1))/2+IVU
              end if
              OP2(JVUTZ) = OP2(JVUTZ)+W_PROD
            end if
            ! Contrib to 2-particle operator, from  dtu Evxyz:
            if (ITABS == IUABS) then
              IVX = IVABS+NASHT*(IXABS-1)
              IYZ = IYABS+NASHT*(IZABS-1)
              if (IVX >= IYZ) then
                JVXYZ = (IVX*(IVX-1))/2+IYZ
              else
                JVXYZ = (IYZ*(IYZ-1))/2+IVX
              end if
              OP2(JVXYZ) = OP2(JVXYZ)+W_PROD
              ! Contrib to 1-particle operator, from  dtu dyx Evz:
              if (IYABS == IXABS) OP1(IVABS,IZABS) = OP1(IVABS,IZABS)+W_PROD
            end if
          end do
        end do
      end do
    end do
    ! Extra sectioning loop added...
  end do
  ! Deallocate temporary space:
  call mma_deallocate(W1)
  call mma_deallocate(W2)
  call mma_deallocate(WPROD)
end do

end subroutine MKWWOPC
