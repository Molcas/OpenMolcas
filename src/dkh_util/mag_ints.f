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
* Copyright (C) 2019, Thomas J. Duignan                                *
*               2021, Rulin Feng                                       *
************************************************************************
C
C----------------------------------------------------------------------|
C
***************************************************************
*
*     An interface for handling magnetic integrals from Gen1Int
*
***************************************************************
      subroutine copy_mag_ints(natoms)
      implicit none
#include "WrkSpc.fh"
      integer, dimension(1) :: IDUM
      integer nmag, irc, iopt, icomp, toper, iat, jtri
      integer iscrt, ncomp, natoms, lu_one
      character*8 Label
*
* Copy magnetic integrals from ONEREL to ONEINT
*
      lu_one=2
      iopt = 0
      irc = -1
      call OpnOne(irc,iopt,'ONEREL',lu_one)
      if (irc.ne.0) goto 9999

      ! Primitive integrals stored on ONEREL
      Label='MAGXP  1'
      iOpt=1
      iComp=1
      tOper=255
      ! Integral dimensions
      call iRdOne(irc,iOpt,Label,iComp,idum,tOper)
      if (irc.ne.0) goto 9999
      nmag=idum(1)
      ! Some scratch space
      call GetMem('scratch ','ALLO','REAL',iscrt,nmag+4)
      iOpt=0
      nComp=9
      do iat=1,nAtoms
        do jtri=1,2
          if (jtri.eq.1) then
            write(Label,'(A,I3)') 'MAGXP',iat
          else
            write(Label,'(A,I3)') 'MAGPX',iat
          endif
          do iComp=1,nComp
            ! Read the primitives from ONEREL
            call RdOne(iRC,iOpt,Label,iComp,Work(iscrt),tOper)
            If (iRC.ne.0) goto 9999
            ! Close ONEREL
            call ClsOne(iRC,iOpt)
            ! Open ONEINT
            call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
            If (iRC.ne.0) goto 9999
            ! Write the primitives to ONEINT ?
            call WrOne(iRC,iOpt,Label,iComp,Work(iscrt),tOper)
            call ClsOne(iRC,iOpt)
            call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
            If (iRC.ne.0) goto 9999
          enddo
        enddo
      enddo
      call GetMem('scratch ','FREE','REAL',iscrt,nmag+4)
      call ClsOne(iRC,iOpt)
      return

 9999 Continue
      Write (6,*) ' *** Error in subroutine Copy_Mag_ints ***'
      Write (6,'(A,A)') '     Label = ', Label
      Call Abend
      End


*
***************************************************************
*


      subroutine merge_mag_ints(nb, jz, lt, ut, dotran)
*
*     Splice together square matrices stored artificially as
*     lower triangular matrices.  Result is that upper triangular
*     portion of lt is set to ut and ut is then transposed so
*     the "upper triangular" portion is in the "lower triangular"
*     elements to circumvent Molcas' internal integral handling.
*
      implicit none
      integer nb, jz, im, jm, km, lm
      logical dotran
      real*8 lt(jz), ut(jz)

      do im = 1, nb
        do jm = 1, nb
          km = (im - 1) * nb + jm
          if (im.le.jm) lt(km) = ut(km)
        enddo
      enddo

      if (dotran) then
        do im = 1, nb
          do jm = 1, nb
            km = (im - 1) * nb + jm
            lm = (jm - 1) * nb + im
            ut(km) = lt(lm)
          enddo
        enddo
      else
        do im = 1, jz
          ut(im) = lt(im)
        enddo
      endif

      end subroutine merge_mag_ints
