************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2016, Sebastian Wouters                                *
*               2016, Quan Phung                                       *
************************************************************************

#include "compiler_features.h"

#ifdef _ENABLE_CHEMPS2_DMRG_
      Subroutine mkfg3chemps2(mkF,NLEV,G1,F1,G2,F2,G3,F3,idxG3,NG3)
      use Symmetry_Info, only: Mul
      use gugx, only: SGS
      use caspt2_module, only: jState, nActel, EPSA, mState
      use caspt2_module, only: NG3
      use definitions, only: iwp, wp, Byte, u6
      IMPLICIT NONE

      LOGICAL(KIND=IWP), INTENT(IN) :: mkF
      INTEGER(KIND=IWP), INTENT(IN) :: NLEV, NG3
      REAL(KIND=WP), INTENT(OUT) ::G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)
      REAL(KIND=WP), INTENT(OUT) ::F1(NLEV,NLEV),F2(NLEV,NLEV,NLEV,NLEV)
      REAL(KIND=WP), INTENT(OUT) :: G3(NG3), F3(nG3)
      INTEGER(KIND=BYTE), INTENT(IN) :: idxG3(6,nG3)

      INTEGER(KIND=IWP) IY,IZ,IW
      INTEGER(KIND=IWP) IYSYM,IXYSYM
      INTEGER(KIND=IWP) NAC4

      If(NACTEL.GT.1) Then
        NAC4 = NLEV * NLEV * NLEV * NLEV
        Call chemps2_load2pdm( nlev, G2, MSTATE(JSTATE) )
        Call two2onerdm( nlev, NACTEL, G2, G1 )
      Else
        write(u6,*) "FATAL ERROR: DMRG-CASPT2 with
     & CHEMPS2 does not work with NACTEL=1"
      End If

! Double checked with CheMPS2::CASPT2::create_f_dots()
      Do iz=1,nlev
        iySym=SGS%ism(iz)
        Do iy=1,nlev
          ixySym=Mul(SGS%ism(iy),iySym)
          If(mkF.AND.ixySym.EQ.1) Then
            F1(iy,iz) = 0.0
            Do iw=1,nlev
              F1(iy,iz)=F1(iy,iz)+G2(iw,iw,iy,iz)*EPSA(iw)
            End Do
          End If
        End Do
      End Do

      If(NACTEL.GE.3) THEN

        If (mkF) call chemps2_load3pdm( nlev, idxG3, NG3, F3, .false.,
     &                                  EPSA, F2, MSTATE(JSTATE) )

        call chemps2_load3pdm( nlev, idxG3, NG3, G3, mkF , EPSA,
     &                         F2, MSTATE(JSTATE) )

      End If

      End Subroutine mkfg3chemps2

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#     include "macros.fh"
      subroutine empty_mkfg3chemps2()
      end subroutine empty_mkfg3chemps2

#endif
