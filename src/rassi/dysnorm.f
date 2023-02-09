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
* Copyright (C) 2020, Bruno Tenorio                                    *
************************************************************************
*  SUBROUTINE DYSNORM
*  PURPOSE: CALCULATE CORRECTED DYSON NORMS FOR CI EXPANSIONS IN A
*  BIORTH. BASIS
************************************************************************
      SUBROUTINE DYSNORM(CMOA,DYSCMO,NORM)

      Implicit Integer (A-Z)

      integer :: isym
      real*8 NORM,NORMSCR,DDOT_
      integer :: nb, nbast, nbast1, nbast2
      real*8, allocatable :: SAO(:),IAO(:),Scr(:),Scr2(:)
      real*8, allocatable :: Dysab(:),Dysab2(:)
      character(len=8) :: Label
      integer :: iOff1, iOff2
      integer :: iOpt,iSyLbl,iRc
      integer :: IC,istca,istcb,ist,ista,istcc,istc,ndys
      REAL*8 DYSCMO(*),CMOA(*)
      integer :: istcmo(8), istao(8), istacc(8)
      Integer no1,nb1,nscr,isy1
      EXTERNAL DDOT_
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "stdalloc.fh"

C============================================================
      nbast=0
      nbast1=0
      nbast2=0
      do isym=1,nsym
        nb=NBASF(isym)
        nbast=nbast+nb
        nbast1=nbast1+(nb*(nb+1))/2
        nbast2=nbast2+nb**2
      end do

      call mma_allocate(DYSAB,NOSHT)
      CALL DCOPY_(NOSHT,DYSCMO,1,DYSAB,1)

      call mma_allocate(SAO,NBAST1)
      call mma_allocate(IAO,NBAST2)
      IAO=0.0D0

      iRc=-1
      iOpt=6
      IC=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,IC,SAO,iSyLbl)
      iOff1 = 0
      iOff2 = 0
      Do iSym = 1,nSym
        nb = nBasf(iSym)
        If ( nb.gt.0 ) then
          call mma_allocate(Scr,nb*nb)
          scr=0.0D0
          Call Square(SAO(1+iOff1),Scr,1,nb,nb)
          CALL DCOPY_(nb*nb,Scr,1,IAO(iOff2+1),1)
          call mma_deallocate(Scr)
        end if
        iOff1 = iOff1 + (nb*nb+nb)/2
        iOff2 = iOff2 + (nb*nb)
      end do
      call mma_deallocate(SAO)

C============================================================

      IST=1
      ISTA=1
      ISTCC=1
      DO ISY1=1,NSYM
        ISTCMO(ISY1)=IST
        ISTAO(ISY1)=ISTA
        ISTACC(ISY1)=ISTCC
        NO1=NOSH(ISY1)
        IST=IST+NO1*NBASF(ISY1)
        ISTA=ISTA+(NBASF(ISY1)*NBASF(ISY1) )
        ISTCC=ISTCC+NO1*NO1
      END DO

      NSCR=NCMO
      NDYS=1
      DO ISY1=1,NSYM
        ISTCB=ISTCMO(ISY1)
        ISTCA=ISTAO(ISY1)
        ISTC=ISTACC(ISY1)
        NO1=NOSH(ISY1)
        NB1=NBASF(ISY1)
        if (NB1*NO1 == 0) cycle

        call mma_allocate(scr,nscr)
        call mma_allocate(scr2,nscr)
        Scr=0.0D0
        Scr2=0.0D0

        CALL DGEMM_('N','N', NB1, NO1, NB1, 1.0D0,
     &                 IAO(ISTCA),NB1, CMOA(ISTCB), NB1,
     &         0.0D0, Scr(ISTCB), NB1)

        CALL DGEMM_('T','N', NO1, NO1, NB1, 1.0D0,
     &                 CMOA(ISTCB),NB1, Scr(ISTCB), NB1,
     &         0.0D0, Scr2(ISTC), NO1)

! Src2 is the M matrix
! Proceed to compute norm=Dab*(Dab*M)
! Where Dab*M=Dysab2

        call mma_allocate(DYSAB2,NO1)
        DYSAB2=0.0D0
        CALL DGEMV_('N', NO1, NO1, 1.0D0, Scr2(ISTC), NO1,
     &              DYSAB(NDYS),1, 0.0D0, DYSAB2, 1)

        normscr= DDOT_(NO1, DYSAB(NDYS), 1, DYSAB2, 1)
        norm=norm+normscr
        call mma_deallocate(DYSAB2)

        NDYS=NDYS+NO1
        Call mma_deallocate(Scr)
        Call mma_deallocate(Scr2)
      END DO

      Call mma_deallocate(IAO)
      call mma_deallocate(DYSAB)

      RETURN

      END
