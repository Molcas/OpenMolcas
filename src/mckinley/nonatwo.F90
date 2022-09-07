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

subroutine NONATWO( &
#                  define _CALLING_
#                  include "grd_mck_interface.fh"
                  )
!***********************************************************************
! OBJECT: TO COMPUTE THE 2ND DERIVATIVE NONADIABATIC COUPLING
! INTEGRALS, OF TYPE
!     < D/DX CHI_1 | D/DX CHI_2 >
!
!     AUTHOR: PER AKE MALMQVIST, MAX PLANCK INSTITUT F ASTROPHYSIK
!             GARCHING, MUENCHEN NOV 2000
!     AFTER PROGRAMMING PATTERN ESTABLISHED BY ROLAND LINDH
!
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Her_RW, only: HerR, HerW, iHerR, iHerW
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_mck_interface.fh"
integer(kind=iwp) :: IBETA, IP, IPALPH, IPAXYZ, IPBETA, IPBXYZ, IPRNXYZ, IPRXYZ, IPSCRT, NIP
logical(kind=iwp) :: ABEQ(3)

#include "macros.fh"
unused_var(ZInv)
unused_var(lOper)
unused_var(iDCnt)
unused_var(iStabM)
unused_var(nStabM)

ABEQ(1) = A(1) == RB(1)
ABEQ(2) = A(2) == RB(2)
ABEQ(3) = A(3) == RB(3)

NIP = 1
IPAXYZ = NIP
NIP = NIP+NZETA*3*NHER*(LA+2)
IPBXYZ = NIP
NIP = NIP+NZETA*3*NHER*(LB+2)
IPRXYZ = NIP
NIP = NIP+NZETA*3*NHER*(NORDOP+1)
IPRNXYZ = NIP
NIP = NIP+NZETA*3*(LA+2)*(LB+2)*(NORDOP+1)
IPALPH = NIP
NIP = NIP+NZETA
IPBETA = NIP
NIP = NIP+NZETA
IPSCRT = NIP
NIP = NIP+nTri_Elem1(LA)*nTri_Elem1(LB)*NZETA*2

if (NIP-1 > NARR) then
  write(u6,*) ' NONATWO: Too small array.'
  write(u6,*) ' Submitted array size NARR=',NARR
  write(u6,*) ' Needed size at least NIP =',NIP
  call Abend()
end if

! COMPUTE THE CARTESIAN VALUES OF THE BASIS FUNCTIONS ANGULAR PART
call CRTCMP(ZETA,P,NZETA,A,ARRAY(IPAXYZ),LA+1,HerR(iHerR(NHER)),NHER,ABEQ)
call CRTCMP(ZETA,P,NZETA,RB,ARRAY(IPBXYZ),LB+1,HerR(iHerR(NHER)),NHER,ABEQ)

!PAM: WILL WE NEED THIS??
! COMPUTE THE CONTRIBUTION FROM THE MULTIPOLE MOMENT OPERATOR
ABEQ(1) = .false.
ABEQ(2) = .false.
ABEQ(3) = .false.
call CRTCMP(ZETA,P,NZETA,CCOOR,ARRAY(IPRXYZ),NORDOP,HerR(iHerR(NHER)),NHER,ABEQ)

! COMPUTE THE PRIMITIVE 1-DIMENSIONAL OVERLAP INTEGRALS.
call ASSMBL(ARRAY(IPRNXYZ),ARRAY(IPAXYZ),LA+1,ARRAY(IPRXYZ),NORDOP,ARRAY(IPBXYZ),LB+1,NZETA,HerW(iHerW(NHER)),NHER)

! COMBINE THE CARTESIAN COMPONENTS OF THE 2DC MATRIX ELEMENTS
IP = IPALPH
do IBETA=1,NBETA
  ARRAY(IP:IP+NALPHA-1) = ALPHA
  IP = IP+NALPHA
end do
IP = IPBETA
do IBETA=1,NBETA
  ARRAY(IP:IP+NALPHA-1) = BETA(IBETA)
  IP = IP+NALPHA
end do
call CMBN2DC(ARRAY(IPRNXYZ),NZETA,LA,LB,ZETA,RKAPPA,ARRAY(IPSCRT),ARRAY(IPALPH),ARRAY(IPBETA),IFGRAD)

! SYMMETRY ADAPT THE 2ND DERIVATIVE COUPLING INTEGRALS
call SYMADO_MCK(ARRAY(IPSCRT),NZETA*nTri_Elem1(LA)*nTri_Elem1(LB),rFinal,NROP,nOP,INDGRD,IU,IV,IFGRAD,IDCAR,TRANS)

return

end subroutine NONATWO
