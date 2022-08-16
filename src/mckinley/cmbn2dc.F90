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

subroutine CMBN2DC(RNXYZ,NZETA,LA,LB,ZETA,RKAPPA,final,ALPHA,BETA,IFGRAD)
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

implicit real*8(A-H,O-Z)
real*8 final(NZETA,(LA+1)*(LA+2)/2,(LB+1)*(LB+2)/2,1), ZETA(NZETA), RKAPPA(NZETA), BETA(NZETA), RNXYZ(NZETA,3,0:LA+1,0:LB+1), &
       ALPHA(NZETA)
logical IFGRAD(3)
! STATEMENT FUNCTION FOR CARTESIAN INDEX
IND(IXYZ,IX,IZ) = (IXYZ-IX)*(IXYZ-IX+1)/2+IZ+1

! PREFACTOR FOR THE PRIMITIVE OVERLAP MATRIX
do IZETA=1,NZETA
  RKAPPA(IZETA) = RKAPPA(IZETA)*(ZETA(IZETA)**(-1.5d0))
end do

! LOOP STRUCTURE FOR THE CARTESIAN ANGULAR PARTS
do IXA=0,LA
  IYAMAX = LA-IXA
  do IXB=0,LB
    IYBMAX = LB-IXB
    do IYA=0,IYAMAX
      IZA = LA-IXA-IYA
      IPA = IND(LA,IXA,IZA)
      do IYB=0,IYBMAX
        IZB = LB-IXB-IYB
        IPB = IND(LB,IXB,IZB)

        ! COMBINE 1-DIM PRIMITIVE OVERLAP INTEGRALS
        if (IFGRAD(1)) then
          ! COMPUTE INTEGRALS TYPE <D/DX,D/DX>
          do IZETA=1,NZETA
            DIFFX = 4d0*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,1,IXA+1,IXB+1)
            if (IXB > 0) then
              DIFFX = DIFFX-2d0*ALPHA(IZETA)*dble(IXB)*RNXYZ(IZETA,1,IXA+1,IXB-1)
              if (IXA > 0) then
                DIFFX = DIFFX+dble(IXA*IXB)*RNXYZ(IZETA,1,IXA-1,IXB-1)
              end if
            end if
            if (IXA > 0) then
              DIFFX = DIFFX-dble(2*IXA)*BETA(IZETA)*RNXYZ(IZETA,1,IXA-1,IXB+1)
            end if
            OVLY = RNXYZ(IZETA,2,IYA,IYB)
            OVLZ = RNXYZ(IZETA,3,IZA,IZB)
            final(IZETA,IPA,IPB,1) = RKAPPA(IZETA)*DIFFX*OVLY*OVLZ
          end do
        end if
        if (IFGRAD(2)) then
          ! COMPUTE INTEGRALS TYPE <D/DY,D/DY>
          do IZETA=1,NZETA
            DIFFY = 4d0*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,2,IYA+1,IYB+1)
            if (IYB > 0) then
              DIFFY = DIFFY-2d0*ALPHA(IZETA)*dble(IYB)*RNXYZ(IZETA,2,IYA+1,IYB-1)
              if (IYA > 0) then
                DIFFY = DIFFY+dble(IYA*IYB)*RNXYZ(IZETA,2,IYA-1,IYB-1)
              end if
            end if
            if (IYA > 0) then
              DIFFY = DIFFY-dble(2*IYA)*BETA(IZETA)*RNXYZ(IZETA,1,IYA-1,IYB+1)
            end if
            OVLX = RNXYZ(IZETA,1,IXA,IXB)
            OVLZ = RNXYZ(IZETA,3,IZA,IZB)
            final(IZETA,IPA,IPB,1) = RKAPPA(IZETA)*OVLX*DIFFY*OVLZ
          end do
        end if
        if (IFGRAD(1)) then
          ! COMPUTE INTEGRALS TYPE <D/DZ,D/DZ>
          do IZETA=1,NZETA
            DIFFZ = 4d0*ALPHA(IZETA)*BETA(IZETA)*RNXYZ(IZETA,1,IZA+1,IZB+1)
            if (IZB > 0) then
              DIFFZ = DIFFZ-2d0*ALPHA(IZETA)*dble(IZB)*RNXYZ(IZETA,1,IZA+1,IZB-1)
              if (IZA > 0) then
                DIFFZ = DIFFZ+dble(IZA*IZB)*RNXYZ(IZETA,1,IZA-1,IZB-1)
              end if
            end if
            if (IZA > 0) then
              DIFFZ = DIFFZ-dble(2*IZA)*BETA(IZETA)*RNXYZ(IZETA,1,IZA-1,IZB+1)
            end if
            OVLX = RNXYZ(IZETA,1,IXA,IXB)
            OVLY = RNXYZ(IZETA,2,IYA,IYB)
            final(IZETA,IPA,IPB,1) = RKAPPA(IZETA)*OVLX*OVLY*DIFFZ
          end do
        end if

        ! END OF LOOP NEST OVER CARTESIAN ANGULAR COMPONENT
      end do
    end do
  end do
end do

return

end subroutine CMBN2DC
