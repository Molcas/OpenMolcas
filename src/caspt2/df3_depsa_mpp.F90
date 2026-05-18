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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

#ifdef _MOLCAS_MPP_
      Subroutine DF3_DEPSA_MPP(NG3,NASHT,DF3,DEPSA,lg_S,idxG3)

      use Symmetry_Info, only: Mul
      USE SUPERINDEX, only: KTUV
      USE Para_Info, only: nProcs
      use definitions, only: wp, iwp, byte
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: IASYM, NTUVES

      implicit none

#include "global.fh"
#include "mafdecls.fh"

      integer(kind=iwp), intent(in) :: NG3, NASHT, lg_S
      real(kind=wp), intent(in) :: DF3(NG3)
      real(kind=wp), intent(inout) :: DEPSA(NASHT,NASHT)
      integer(kind=byte), intent(in) :: idxG3(6,NG3)

      real(kind=wp), allocatable :: WRK(:)
      integer(kind=iwp) :: isym, iRank, ILO, IHI, JLO, JHI, NROW, NCOL, &
     &  iG3, iT, iU, iV, iX, iY, iZ, iST, iSU, iSV, iSX, iSY, iSZ,      &
     &  ituvs, ixyzs, ISUP1, ISUP2, JSUP1, JSUP2, IW, NSEQ
      real(kind=wp) :: F3VAL

      !! do depsa
      isym=1
      Do iRank = 0, NPROCS-1
        CALL GA_Distribution(lg_S,iRank,ILO,IHI,JLO,JHI)
        NROW=iHi-iLo+1
        NCOL=jHi-jLo+1
        Call mma_allocate(WRK,NROW*NCOL,Label='WRK')
        Call GA_Get(lg_S,iLo,iHi,jLo,jHi,WRK,NROW)
        DO iG3=1,NG3
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=Mul(IST,Mul(ISU,ISV))
          ixyzs=Mul(ISX,Mul(ISY,ISZ))
          if(ituvs /= ixyzs) cycle

          F3VAL = DF3(iG3)

          ISUP1=KTUV(iV,iU,iT)-nTUVES(iSYM)
          If (ISUP1 >= iLo .and. ISUP1 <= iHi) Then
            Do iW = 1, nAshT
              JSUP1=KTUV(iX,iW,iZ)-nTUVES(iSYM)
              NSEQ=1+ISUP1-iLo+NROW*(JSUP1-jLo)
              DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*WRK(NSEQ)
            End Do
          End If
          JSUP2=KTUV(iX,iY,iZ)-nTUVES(iSYM)
          If (JSUP2 >= iLo .and. JSUP2 <= iHi) Then
            Do iW = 1, nAshT
              ISUP2=KTUV(iV,iW,iT)-nTUVES(iSYM)
              NSEQ=1+JSUP2-iLo+NROW*(ISUP2-jLo)
              DEPSA(iW,iU) = DEPSA(iW,iU) - F3VAL*WRK(NSEQ)
            End Do
          End If
        END DO
        call mma_deallocate(WRK)
      End Do

      Return

      End Subroutine DF3_DEPSA_MPP

#elif defined (NAGFOR)
      ! Some compilers do not like empty files
      subroutine empty_DF3_DEPSA_MPP()
      end subroutine empty_DF3_DEPSA_MPP
#endif
