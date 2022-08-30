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
! Copyright (C) 2000, Per Ake Malmqvist                                *
!***********************************************************************

subroutine CMBN2DC(RNXYZ,NZETA,LA,LB,ZETA,RKAPPA,RFINAL,ALPHA,BETA,IFGRAD)
!***********************************************************************
!
! OBJECT: COMPUTE THE SECOND DERIVATIVE NON-ADIABATIC COUPLING
! MATRIX ELEMENTS, OF TYPE < D/DX CHI_1 | D/DX CHI_2 >
! WITH DIFFERENTIATION WRT NUCLEAR COORDINATES
!
!     AUTHOR: PER AKE MALMQVIST NOV 2000
!             DEPT. OF THEORETICAL CHEMISTRY,
!             UNIVERSITY OF LUND, SWEDEN
!             FOLLOWING THE PATTERN OF R. LINDH,
!             SAME PLACE.
!
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two, Four, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NZETA, LA, LB
real(kind=wp) :: RNXYZ(NZETA,3,0:LA+1,0:LB+1), ZETA(NZETA), RKAPPA(NZETA), RFINAL(NZETA,nTri_Elem1(LA),nTri_Elem1(LB),1), &
                 ALPHA(NZETA), BETA(NZETA)
logical(kind=iwp) :: IFGRAD(3)
integer(kind=iwp) :: IPA, IPB, IXA, IXB, IYA, IYAMAX, IYB, IYBMAX, IZA, IZB, IZETA
real(kind=wp) :: DIFFX, DIFFY, DIFFZ, OVLX, OVLY, OVLZ

! PREFACTOR FOR THE PRIMITIVE OVERLAP MATRIX
RKAPPA(:) = RKAPPA*(ZETA**(-OneHalf))

! LOOP STRUCTURE FOR THE CARTESIAN ANGULAR PARTS
do IXA=0,LA
  IYAMAX = LA-IXA
  do IXB=0,LB
    IYBMAX = LB-IXB
    do IYA=0,IYAMAX
      IZA = LA-IXA-IYA
      IPA = C_Ind(LA,IXA,IZA)
      do IYB=0,IYBMAX
        IZB = LB-IXB-IYB
        IPB = C_Ind(LB,IXB,IZB)

        ! COMBINE 1-DIM PRIMITIVE OVERLAP INTEGRALS
        if (IFGRAD(1)) then
          ! COMPUTE INTEGRALS TYPE <D/DX,D/DX>
          do IZETA=1,NZETA
            DIFFX = Four*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,1,IXA+1,IXB+1)
            if (IXB > 0) then
              DIFFX = DIFFX-Two*ALPHA(IZETA)*real(IXB,kind=wp)*RNXYZ(IZETA,1,IXA+1,IXB-1)
              if (IXA > 0) then
                DIFFX = DIFFX+real(IXA*IXB,kind=wp)*RNXYZ(IZETA,1,IXA-1,IXB-1)
              end if
            end if
            if (IXA > 0) then
              DIFFX = DIFFX-real(2*IXA,kind=wp)*BETA(IZETA)*RNXYZ(IZETA,1,IXA-1,IXB+1)
            end if
            OVLY = RNXYZ(IZETA,2,IYA,IYB)
            OVLZ = RNXYZ(IZETA,3,IZA,IZB)
            RFINAL(IZETA,IPA,IPB,1) = RKAPPA(IZETA)*DIFFX*OVLY*OVLZ
          end do
        end if
        if (IFGRAD(2)) then
          ! COMPUTE INTEGRALS TYPE <D/DY,D/DY>
          do IZETA=1,NZETA
            DIFFY = Four*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,2,IYA+1,IYB+1)
            if (IYB > 0) then
              DIFFY = DIFFY-Two*ALPHA(IZETA)*real(IYB,kind=wp)*RNXYZ(IZETA,2,IYA+1,IYB-1)
              if (IYA > 0) then
                DIFFY = DIFFY+real(IYA*IYB,kind=wp)*RNXYZ(IZETA,2,IYA-1,IYB-1)
              end if
            end if
            if (IYA > 0) then
              DIFFY = DIFFY-real(2*IYA,kind=wp)*BETA(IZETA)*RNXYZ(IZETA,1,IYA-1,IYB+1)
            end if
            OVLX = RNXYZ(IZETA,1,IXA,IXB)
            OVLZ = RNXYZ(IZETA,3,IZA,IZB)
            RFINAL(IZETA,IPA,IPB,1) = RKAPPA(IZETA)*OVLX*DIFFY*OVLZ
          end do
        end if
        if (IFGRAD(1)) then
          ! COMPUTE INTEGRALS TYPE <D/DZ,D/DZ>
          do IZETA=1,NZETA
            DIFFZ = Four*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,1,IZA+1,IZB+1)
            if (IZB > 0) then
              DIFFZ = DIFFZ-Two*ALPHA(IZETA)*real(IZB,kind=wp)*RNXYZ(IZETA,1,IZA+1,IZB-1)
              if (IZA > 0) then
                DIFFZ = DIFFZ+real(IZA*IZB,kind=wp)*RNXYZ(IZETA,1,IZA-1,IZB-1)
              end if
            end if
            if (IZA > 0) then
              DIFFZ = DIFFZ-real(2*IZA,kind=wp)*BETA(IZETA)*RNXYZ(IZETA,1,IZA-1,IZB+1)
            end if
            OVLX = RNXYZ(IZETA,1,IXA,IXB)
            OVLY = RNXYZ(IZETA,2,IYA,IYB)
            RFINAL(IZETA,IPA,IPB,1) = RKAPPA(IZETA)*OVLX*OVLY*DIFFZ
          end do
        end if

        ! END OF LOOP NEST OVER CARTESIAN ANGULAR COMPONENT
      end do
    end do
  end do
end do

return

end subroutine CMBN2DC
