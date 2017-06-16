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
#ifdef _ENABLE_CHEMPS2_DMRG_
      Subroutine mkfg3chemps2(IFF,G1,F1,G2,F2,G3,F3,idxG3)
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"

      INTEGER, INTENT(IN) :: IFF
      REAL*8, INTENT(OUT) :: G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: F1(NLEV,NLEV),F2(NLEV,NLEV,NLEV,NLEV)
      REAL*8, INTENT(OUT) :: G3(*), F3(*)
      INTEGER*1, INTENT(IN) :: idxG3(6,*)

      INTEGER IY,IZ,IW
      INTEGER IYSYM,IXYSYM
      INTEGER NAC4

      If(NACTEL.GT.1) Then
        NAC4 = NLEV * NLEV * NLEV * NLEV
        Call chemps2_load2pdm( nlev, G2, MSTATE(JSTATE) )
        Call two2onerdm_bis( nlev, NACTEL, G2, G1 )
      Else
        write(6,*) "FATAL ERROR: DMRG-CASPT2 with
     & CHEMPS2 does not work with NACTEL=1"
      End If

! Double checked with CheMPS2::CASPT2::create_f_dots()
      Do iz=1,nlev
        iySym=ism(iz)
        Do iy=1,nlev
          ixySym=Mul(ism(iy),iySym)
          If(IFF.NE.0.AND.ixySym.EQ.1) Then
            F1(iy,iz) = 0.0
            Do iw=1,nlev
              F1(iy,iz)=F1(iy,iz)+G2(iw,iw,iy,iz)*EPSA(iw)
            End Do
          End If
        End Do
      End Do

      If(NACTEL.GE.3) THEN

        call chemps2_load3pdm( nlev, idxG3, NG3, F3, .false., EPSA,
     &                         F2, MSTATE(JSTATE) )
        call chemps2_load3pdm( nlev, idxG3, NG3, G3, .true. , EPSA,
     &                         F2, MSTATE(JSTATE) )

      End If

      Return
      End
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_mkfg3chemps2()
      End
#endif
