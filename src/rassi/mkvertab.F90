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

subroutine mkVERTAB(NACTEL,M2SPIN,LSYM,NPART,NGASORB,NGASLIM,ISSTAB,NFSB0,NRDETS0,NFSB,NRDETS)

use stdalloc, only: mma_allocate, mma_deallocate
use rassi_global_arrays, only: FSBARR
use cntrl, only: MORSBITS
use Symmetry_Info, only: nSym => nIrrep, MUL
use Definitions, only: u6

implicit none
integer ISSTAB(*)
integer NASPRT, NACTEL, M2SPIN, LSYM
integer NPART
integer NGASORB(0:NSYM,0:NPART)
integer NGASLIM(2,NPART)
integer NFSB, NFSB0, NRDETS, NRDETS0
integer NVERT, IVER1, IVER2
integer IPOS, ISPART, NPC1, M2C1, ISC1
integer ISST, KPOP, KMS2, KSYM
integer MX1, MN1, MX2, MN2, IPART, NOREM, MINEL, MAXEL
integer ISPEND, NGO, NSP, ISPSTA, MNL, MXL, MNR, MXR, NO
integer NPC2, M2C2, ISC2
integer NSSTP, KSSTP, KPOPMAX, KMS2MAX, NPOP, MS2, NSSTPTR
integer N, IEND, ISTA, IADDR, ISAVE, NEXT
integer MS2MIN, MS2MAX, MSFACT, NVIDX, IVIDX, IVERT
integer LEVUP, LEVDWN, IADDR1, IADDR2, NARC
integer NARCVRT, NVERTAB, NDWNTAB, NARCLEV
integer NFSBARR, IVUP, LOWEST, ISW, NSW
integer KPOPLIM
integer IARC, IA, ISUM, KEEP
integer NSBS
integer NSDBLK
integer LIMARR(2,0:50)
integer ITRY(50)
integer, allocatable :: SSTPTR(:), SSTARR(:), VIDX(:)
integer, allocatable :: VERTAB(:), DWNTAB(:), LTDA(:)
integer, allocatable :: WEIGHT(:), SWITCH(:)

!----------------------------------------------------------------
! Unbutton the substring table:
NSYM = ISSTAB(4)
NASPRT = ISSTAB(5)
NSSTP = ISSTAB(7)
KSSTP = 15

!TEST write(u6,*) ' Test print in VERTAB.'
!TEST write(u6,*) ' NSYM  :',NSYM
!TEST write(u6,*) ' NASPRT:',NASPRT
!TEST write(u6,*) ' NSSTP :',NSSTP
!TEST write(u6,*) ' Substrings:'
!TEST do isst=1,nsstp
!TEST   NSBS = ISSTAB(KSSTP+0+5*(ISST-1))
!TEST   NPOP = ISSTAB(KSSTP+1+5*(ISST-1))
!TEST   KSYM = ISSTAB(KSSTP+2+5*(ISST-1))
!TEST   MS2 = ISSTAB(KSSTP+3+5*(ISST-1))
!TEST   ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
!TEST   write(u6,'(1x,6i8)') isst,nsbs,npop,ksym,ms2,ispart
!TEST end do
! Now the properties of the substring types can be picked up as
!       NSBS   = ISSTAB(KSSTP+0+5*(ISST-1))
!       NPOP   = ISSTAB(KSSTP+1+5*(ISST-1))
!       ISYM   = ISSTAB(KSSTP+2+5*(ISST-1))
!       MS2    = ISSTAB(KSSTP+3+5*(ISST-1))
!       ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
!---------------------------------------------------------------
! Determine KPOPMAX = Max population of any subpartition
! Determine KMS2MAX = Max spin proj in any subpartition
KPOPMAX = 0
KMS2MAX = 0
do ISST=1,NSSTP
  NPOP = ISSTAB(KSSTP+1+5*(ISST-1))
  MS2 = ISSTAB(KSSTP+3+5*(ISST-1))
  KPOPMAX = max(NPOP,KPOPMAX)
  KMS2MAX = max(MS2,KMS2MAX)
end do
!TEST write(u6,'(1x,a,8i8)') 'KPOPMAX:',KPOPMAX
!TEST write(u6,'(1x,a,8i8)') 'KMS2MAX:',KMS2MAX
!----------------------------------------------------------------
! Using the GAS restrictions, set up an array with min and max nr of
! electrons (cumulative) up to and including subpartition ISPART:
MX1 = 0
MN1 = 0
MX2 = NACTEL
MN2 = NACTEL
do IPART=1,NPART
  NOREM = NGASORB(0,IPART)
  MINEL = max(0,NGASLIM(1,IPART))
  MAXEL = min(NOREM,NGASLIM(2,IPART))
  MX2 = MX2-MINEL
  MN2 = MN2-MAXEL
end do
ISPEND = 0

LIMARR(1,0) = 0
LIMARR(2,0) = 0
do IPART=1,NPART
  NGO = NGASORB(0,IPART)
  MINEL = max(0,NGASLIM(1,IPART))
  MAXEL = min(NGO,NGASLIM(2,IPART))
  NSP = (NGO+MORSBITS-1)/MORSBITS
  ISPSTA = ISPEND+1
  ISPEND = ISPEND+NSP
  MNL = max(MN1,MN2)
  MXL = min(MX1,MX2)
  MN1 = MN1+MINEL
  MX1 = MX1+MAXEL
  MN2 = MN2+MAXEL
  MX2 = MX2+MINEL
  MNR = max(MN1,MN2)
  MXR = min(MX1,MX2)
  NOREM = NGO
  do ISPART=ISPSTA,ISPEND
    NO = min(NOREM,MORSBITS)
    NOREM = NOREM-NO
    LIMARR(1,ISPART) = max(MNL,MNR-NOREM)
    LIMARR(2,ISPART) = min(MXL+NGO-NOREM,MXR)
  end do
end do
!TEST write(u6,*) ' The LIMARR array:'
!TEST write(u6,'(1x,8i5)') (limarr(1,ispart),ispart=0,nasprt)
!TEST write(u6,'(1x,8i5)') (limarr(2,ispart),ispart=0,nasprt)

!----------------------------------------------------------------
! First, set up an array containing substring types.
! It would be nice to be able to quickly loop through just those
! substrings that have a limited population, and belong to a
! specific subpartition. To achieve this, we need a temporary
! array with index limits (0:KPOPMAX+1,1:NASPRT), and one with
! index limits (1:NSSTP).
NSSTPTR = (KPOPMAX+2)*NASPRT
call mma_allocate(SSTPTR,NSSTPTR,Label='SSTPTR')
do ISPART=1,NASPRT
  do KPOP=1,KPOPMAX+2
    SSTPTR(KPOP+(KPOPMAX+2)*(ISPART-1)) = 0
  end do
end do
call mma_allocate(SSTARR,NSSTP,Label='SSTARR')
do ISST=1,NSSTP
  NPOP = ISSTAB(KSSTP+1+5*(ISST-1))
  ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
  N = 1+SSTPTR(1+NPOP+(KPOPMAX+2)*(ISPART-1))
  SSTPTR(1+NPOP+(KPOPMAX+2)*(ISPART-1)) = N
end do
IEND = 0
do ISPART=1,NASPRT
  do KPOP=0,KPOPMAX
    N = SSTPTR(1+KPOP+(KPOPMAX+2)*(ISPART-1))
    ISTA = IEND+1
    IEND = IEND+N
    SSTPTR(1+KPOP+(KPOPMAX+2)*(ISPART-1)) = ISTA
  end do
  SSTPTR(1+KPOPMAX+1+(KPOPMAX+2)*(ISPART-1)) = IEND+1
end do
do ISST=1,NSSTP
  NPOP = ISSTAB(KSSTP+1+5*(ISST-1))
  ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
  IADDR = SSTPTR(1+NPOP+(KPOPMAX+2)*(ISPART-1))
  SSTPTR(1+NPOP+(KPOPMAX+2)*(ISPART-1)) = IADDR+1
  SSTARR(IADDR) = ISST
end do
ISAVE = 1
do ISPART=1,NASPRT
  do KPOP=0,KPOPMAX
    NEXT = SSTPTR(1+KPOP+(KPOPMAX+2)*(ISPART-1))
    SSTPTR(1+KPOP+(KPOPMAX+2)*(ISPART-1)) = ISAVE
    ISAVE = NEXT
  end do
end do
!TEST test print
!TEST IADDR1 = sstptr(1+(kpopmax+2)*(1-1))
!TEST IADDR2 = sstptr(1+kpopmax+1+(kpopmax+2)*(nasprt-1))
!TEST write(u6,'(1x,a,8i8)') 'iaddr1,iaddr2:',iaddr1,iaddr2
!TEST do IADDR=IADDR1,IADDR2-1
!TEST   isst = sstarr(IADDR)
!TEST   kpop = ISSTAB(KSSTP+1+5*(ISST-1))
!TEST   ksym = ISSTAB(KSSTP+2+5*(ISST-1))
!TEST   kms2 = ISSTAB(KSSTP+3+5*(ISST-1))
!TEST write(u6,'(1x,a,8i8)') 'iaddr,isst,kpop,ksym,kms2:',iaddr,isst,kpop,ksym,kms2
!TEST end do
!TEST End of test print

!----------------------------------------------------------------
! Here follows the construction of a graph, to be used for
! producing the FSB table.
! Use a single unique compound index for the vertices
! This way, we can use a temporary counter array to
! identify vertices.
! Subpartition is in 1..NASPRT, so we need levels 0..NASPRT
! Cumulative population is in 0..NACTEL
! Intermediate spin projection is in MS2MIN..MS2MAX:
MS2MIN = (M2SPIN-KMS2MAX*NASPRT)/2
MS2MAX = (M2SPIN+KMS2MAX*NASPRT)/2
! and it runs in steps of two units.
MSFACT = 1+(MS2MAX-MS2MIN)/2
! Intermediate symmetry label is in 1..NSYM
! All in all, the index
!  IND=ISCUM+NSYM*((M2CUM-MS2MIN)/2 + MSFACT*(NPCUM + (NACTEL+1)*ISPART))
! should do the trick.
! Run through the generation procedure twice. First time, we will
! determine the sizes of the vertex and downchain tables, and
! allocate them. Second time, we will fill in the values in these
! tables, and then the temporary counter table can be deallocated.
!----------------------------------
! Allocate temporary counter array:
NVIDX = NSYM*MSFACT*(NACTEL+1)*(NASPRT+1)
call mma_allocate(VIDX,NVIDX,Label='VIDX')
do IVIDX=1,NVIDX
  VIDX(IVIDX) = 0
end do
!----------------------------------
! Initialize top vertex:
IPOS = LSYM+NSYM*((M2SPIN-MS2MIN)/2+MSFACT*(NACTEL+(NACTEL+1)*NASPRT))
VIDX(IPOS) = 1
IVERT = 1
!----------------------------------
! Recursively, find all other reachable vertices:
! Upper and lower level:
do LEVUP=NASPRT,1,-1
  LEVDWN = LEVUP-1
  ! A wasteful loop to pick out vertices on the upper level:
  do NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
    do M2C1=-NPC1,NPC1,2
      if (M2C1 < MS2MIN) goto 130
      if (M2C1 > MS2MAX) goto 130
      do ISC1=1,NSYM
        IPOS = ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(NPC1+(NACTEL+1)*LEVUP))
        IVER1 = VIDX(+IPOS)
        if (IVER1 == 0) goto 120
        ! This is a valid upper vertex on the upper level.
        ! Loop over valid arcs:
        KPOPLIM = min(NPC1,KPOPMAX)
        IADDR1 = SSTPTR(1+(KPOPMAX+2)*(LEVUP-1))
        IADDR2 = SSTPTR(1+KPOPLIM+1+(KPOPMAX+2)*(LEVUP-1))
        do IADDR=IADDR1,IADDR2-1
          ISST = SSTARR(IADDR)
          KPOP = ISSTAB(KSSTP+1+5*(ISST-1))
          KSYM = ISSTAB(KSSTP+2+5*(ISST-1))
          KMS2 = ISSTAB(KSSTP+3+5*(ISST-1))
          ! Just temporary sanity test:
          ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
          if (ISPART /= LEVUP) then
            write(u6,*) ' THIS IS INSANE!'
            call ABEND()
          end if
          NPC2 = NPC1-KPOP
          if (NPC2 < LIMARR(1,LEVDWN)) goto 110
          if (NPC2 > LIMARR(2,LEVDWN)) goto 110
          M2C2 = M2C1-KMS2
          if (M2C2 < MS2MIN) goto 110
          if (M2C2 > MS2MAX) goto 110
          ISC2 = MUL(KSYM,ISC1)
          if (abs(M2C2) > NPC2) goto 110
          if ((NPC2 == 0) .and. (ISC2 /= 1)) goto 110
          ! Inspect the lower vertex: Is it a new one?
          IPOS = ISC2+NSYM*((M2C2-MS2MIN)/2+MSFACT*(NPC2+(NACTEL+1)*LEVDWN))
          IVER2 = VIDX(+IPOS)
          if (IVER2 > 0) goto 110
          ! This is a new vertex. Register its number
          IVERT = IVERT+1
          IVER2 = IVERT
          VIDX(+IPOS) = IVER2
110       continue
        end do
        ! Finished looping over possible arcs.
120     continue
      end do
130   continue
    end do
  end do
  ! Finished looping over possible upper vertices.
end do
! Finished loop over upper level, LEVUP.
!----------------------------------
! Now the index array holds a positive integer IV if the combination
! of properties is such that this vertex could be reached from above.
! We now need to keep just those that can also be reached from below.
! This is accomplished by a similar loop, now in the other direction:
!----------------------------------
NARC = 0
! Upper and lower level:
do LEVUP=1,NASPRT
  LEVDWN = LEVUP-1
  ! A wasteful loop to pick out vertices on the upper level:
  do NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
    do M2C1=-NPC1,NPC1,2
      if (M2C1 < MS2MIN) goto 230
      if (M2C1 > MS2MAX) goto 230
      do ISC1=1,NSYM
        IPOS = ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(NPC1+(NACTEL+1)*LEVUP))
        IVER1 = VIDX(+IPOS)
        if (IVER1 == 0) goto 220
        ! This is a valid upper vertex on the upper level.
        ! Is it reachable from below?
        NARCVRT = 0
        ! Loop over valid arcs:
        KPOPLIM = min(NPC1,KPOPMAX)
        IADDR1 = SSTPTR(1+(KPOPMAX+2)*(LEVUP-1))
        IADDR2 = SSTPTR(1+KPOPLIM+1+(KPOPMAX+2)*(LEVUP-1))
        do IADDR=IADDR1,IADDR2-1
          ISST = SSTARR(IADDR)
          KPOP = ISSTAB(KSSTP+1+5*(ISST-1))
          KSYM = ISSTAB(KSSTP+2+5*(ISST-1))
          KMS2 = ISSTAB(KSSTP+3+5*(ISST-1))
          ! Just temporary sanity test:
          ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
          if (ISPART /= LEVUP) then
            write(u6,*) ' THIS IS INSANE!'
            call ABEND()
          end if
          NPC2 = NPC1-KPOP
          if ((LEVDWN == 0) .and. (NPC2 > 0)) goto 210
          M2C2 = M2C1-KMS2
          if (M2C2 < MS2MIN) goto 210
          if (M2C2 > MS2MAX) goto 210
          ISC2 = MUL(KSYM,ISC1)
          if (abs(M2C2) > NPC2) goto 210
          if ((NPC2 == 0) .and. (ISC2 /= 1)) goto 210
          ! Inspect the lower vertex: Is it reachable?
          IPOS = ISC2+NSYM*((M2C2-MS2MIN)/2+MSFACT*(NPC2+(NACTEL+1)*LEVDWN))
          IVER2 = VIDX(+IPOS)
          if (IVER2 > 0) NARCVRT = NARCVRT+1
210       continue
        end do
        ! Finished looping over possible arcs.
        ! Remove this vertex if it was unreachable.
        if (NARCVRT == 0) then
          IPOS = ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(NPC1+(NACTEL+1)*LEVUP))
          VIDX(+IPOS) = 0
        end if
        NARC = NARC+NARCVRT
220     continue
      end do
230   continue
    end do
  end do
  ! Finished looping over possible upper vertices.
end do
! Finished loop over upper level, LEVUP.
!----------------------------------
! Renumber the vertices:
NVERT = 0
! Upper and lower level:
do LEVUP=NASPRT,0,-1
  ! A wasteful loop to pick out vertices on the upper level:
  do NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
    do M2C1=-NPC1,NPC1,2
      if (M2C1 < MS2MIN) goto 330
      if (M2C1 > MS2MAX) goto 330
      do ISC1=1,NSYM
        IPOS = ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(NPC1+(NACTEL+1)*LEVUP))
        IVER1 = VIDX(+IPOS)
        if (IVER1 == 0) goto 320
        NVERT = NVERT+1
        VIDX(+IPOS) = NVERT
320     continue
      end do
330   continue
    end do
  end do
end do
!----------------------------------
! Now allocate the vertex and downchain tables, and loop again:
NVERTAB = 6*NVERT
call mma_allocate(VERTAB,NVERTAB,Label='VERTAB')
NDWNTAB = 3*NARC
call mma_allocate(DWNTAB,NDWNTAB,Label='DWNTAB')
IARC = 0
!TEST write(u6,'(1x,a,8i8)') 'Nr of vertices NVERT=',NVERT
!TEST write(u6,'(1x,a,8i8)') 'Nr of arcs     NARC =',NARC
! --------------------------------------------------------------
! Allocate LTDA(1:2,1:NASPRT), Level-to-Downarc array:
call mma_allocate(LTDA,2*(NASPRT+1),Label='LTDA')
!----------------------------------
! Upper and lower level:
do LEVUP=NASPRT,1,-1
  LEVDWN = LEVUP-1
  ! Level-to-Downarc array:
  LTDA(1+2*LEVUP) = IARC+1
  LTDA(1+2*LEVUP+1) = 0
  ! Initialize counter of arcs from this upper level:
  NARCLEV = 0
  ! A wasteful loop to pick out vertices on the upper level:
  do NPC1=LIMARR(1,LEVUP),LIMARR(2,LEVUP)
    do M2C1=-NPC1,NPC1,2
      if (M2C1 < MS2MIN) goto 430
      if (M2C1 > MS2MAX) goto 430
      do ISC1=1,NSYM
        IPOS = ISC1+NSYM*((M2C1-MS2MIN)/2+MSFACT*(NPC1+(NACTEL+1)*LEVUP))
        IVER1 = VIDX(+IPOS)
        if (IVER1 == 0) goto 420
        ! This is a valid upper vertex on the upper level.
        VERTAB(+1+6*(IVER1-1)) = LEVUP
        VERTAB(+2+6*(IVER1-1)) = NPC1
        VERTAB(+3+6*(IVER1-1)) = M2C1
        VERTAB(+4+6*(IVER1-1)) = ISC1
        VERTAB(+5+6*(IVER1-1)) = IARC+1
        VERTAB(+6+6*(IVER1-1)) = 0
        ! Initialize counter of arcs from this upper vertex:
        NARCVRT = 0
        ! Loop over valid arcs:
        KPOPLIM = min(NPC1,KPOPMAX)
        IADDR1 = SSTPTR(1+(KPOPMAX+2)*(LEVUP-1))
        IADDR2 = SSTPTR(1+KPOPLIM+1+(KPOPMAX+2)*(LEVUP-1))
        do IADDR=IADDR1,IADDR2-1
          ISST = SSTARR(IADDR)
          KPOP = ISSTAB(KSSTP+1+5*(ISST-1))
          KSYM = ISSTAB(KSSTP+2+5*(ISST-1))
          KMS2 = ISSTAB(KSSTP+3+5*(ISST-1))
          ! Just temporary sanity test:
          ISPART = ISSTAB(KSSTP+4+5*(ISST-1))
          if (ISPART /= LEVUP) then
            write(u6,*) ' THIS IS INSANE!'
            call ABEND()
          end if
          NPC2 = NPC1-KPOP
          if ((LEVDWN == 0) .and. (NPC2 > 0)) goto 410
          M2C2 = M2C1-KMS2
          if (M2C2 < MS2MIN) goto 410
          if (M2C2 > MS2MAX) goto 410
          ISC2 = MUL(KSYM,ISC1)
          if (abs(M2C2) > NPC2) goto 410
          if ((NPC2 == 0) .and. (ISC2 /= 1)) goto 410
          ! Is it a valid arc? See if a lower vertex exists.
          IPOS = ISC2+NSYM*((M2C2-MS2MIN)/2+MSFACT*(NPC2+(NACTEL+1)*LEVDWN))
          IVER2 = VIDX(+IPOS)
          if (IVER2 == 0) goto 410
          ! A valid arc has been found. The upper vertex IVER1 is
          ! joined to the lower vertex IVER2 by an arc associated
          ! with the substring type ISST.
          NARCVRT = NARCVRT+1
          NARCLEV = NARCLEV+1
          IARC = IARC+1
          DWNTAB(+1+3*(IARC-1)) = ISST
          DWNTAB(+2+3*(IARC-1)) = IVER1
          DWNTAB(+3+3*(IARC-1)) = IVER2
410       continue
        end do
        ! Finished looping over possible arcs.
        VERTAB(+6+6*(IVER1-1)) = NARCVRT
420     continue
      end do
430   continue
    end do
  end do
  ! Finished looping over vertices on this level.
  LTDA(1+2*LEVUP+1) = NARCLEV
end do
! Finished loop over upper level, LEVUP.
VERTAB(+1+6*(NVERT-1)) = 0
VERTAB(+2+6*(NVERT-1)) = 0
VERTAB(+3+6*(NVERT-1)) = 0
VERTAB(+4+6*(NVERT-1)) = 1
VERTAB(+5+6*(NVERT-1)) = IARC
VERTAB(+6+6*(NVERT-1)) = 0
LTDA(1) = IARC
LTDA(1+1) = 0
!---------------------------------------------------------
! The graph is finished. we do no longer need these arrays:
call mma_deallocate(SSTPTR)
call mma_deallocate(SSTARR)
call mma_deallocate(VIDX)
!---------------------------------------------------------
! A graph has been constructed. We may construct a weight array:
call mma_allocate(WEIGHT,NVERT,Label='WEIGHT')
do IVERT=1,NVERT-1
  WEIGHT(+IVERT) = 0
end do
WEIGHT(+NVERT) = 1
do LEVUP=1,NASPRT
  LEVDWN = LEVUP-1
  ! Loop over arcs leading down from upper level:
  NARCLEV = LTDA(1+2*LEVUP+1)
  do IA=1,NARCLEV
    IARC = LTDA(1+2*LEVUP)-1+IA
    ISST = DWNTAB(+1+3*(IARC-1))
    IVER1 = DWNTAB(+2+3*(IARC-1))
    IVER2 = DWNTAB(+3+3*(IARC-1))
    ISUM = WEIGHT(+IVER1)+WEIGHT(+IVER2)
    WEIGHT(+IVER1) = ISUM
  end do
end do
! No more use for Level-to-Downarc array:
call mma_deallocate(LTDA)
!---------------------------------------------------------
! Now, for any vertex iv, WEIGHT(+iv) is the total number of
! downwalks that can reach vertex nr nvert (The bottom vertex).
! Traverse the graph. Use a swich array. The graph is traversed
! from the top. At each vertex, one of the available arcs downwards
! is selected. Which arc to select is given by ISWITCH(LEVUP), where
! LEVUP is the level of the upper vertex.
call mma_allocate(SWITCH,NASPRT,Label='SWITCH')
! Initial values: Lowest walk
do LEVUP=1,NASPRT
  SWITCH(+LEVUP) = 1
end do
! The total number of walks in the graph:
NFSB0 = WEIGHT(1)
NFSBARR = (NASPRT+2)*NFSB0
call mma_allocate(FSBARR,NFSBARR,Label='FSBARR')
NRDETS0 = 0
NFSB = 0
NRDETS = 0
!TEST write(u6,*) ' A list of all generated FS blocks:'
500 continue
! Construct this walk. While constructing it, also find out if
! it has a successor.
IVUP = 1
LOWEST = NASPRT+1
do LEVUP=NASPRT,1,-1
  ISW = SWITCH(+LEVUP)
  ! Select arc nr isw from those available to vertex ivup.
  IARC = VERTAB(+5+6*(IVUP-1))-1+ISW
  ! Could it be incremented?
  NSW = VERTAB(+6+6*(IVUP-1))
  !TEST write(u6,'(1x,a,8i8)') 'IVUP,ISW,NSW:',IVUP,ISW,NSW
  if (NSW > ISW) LOWEST = LEVUP
  ! Consult the downarc table:
  ISST = DWNTAB(+1+3*(IARC-1))
  IVER1 = DWNTAB(+2+3*(IARC-1))
  IVER2 = DWNTAB(+3+3*(IARC-1))
  !TEST write(u6,'(1x,a,8i8)') 'ISST,IVER1,IVER2:',ISST,IVER1,IVER2
  if (IVER1 /= IVUP) then
    write(u6,*) ' OOOOPS! THIS SHOULD NOT HAPPEN.'
    call ABEND()
  end if
  ITRY(LEVUP) = ISST
  IVUP = IVER2
end do
! Check GAS restrictions:
NSDBLK = 1
ISPEND = 0
KEEP = 1
do IPART=1,NPART
  NGO = NGASORB(0,IPART)
  NSP = (NGO+MORSBITS-1)/MORSBITS
  ISPSTA = ISPEND+1
  ISPEND = ISPEND+NSP
  NOREM = NGO
  ISUM = 0
  do ISPART=ISPSTA,ISPEND
    ISST = ITRY(ISPART)
    NSBS = ISSTAB(KSSTP+0+5*(ISST-1))
    NPOP = ISSTAB(KSSTP+1+5*(ISST-1))
    ISUM = ISUM+NPOP
    NSDBLK = NSDBLK*NSBS
  end do
  if (ISUM < NGASLIM(1,IPART)) KEEP = 0
  if (ISUM > NGASLIM(2,IPART)) KEEP = 0
end do
NRDETS0 = NRDETS0+NSDBLK
if (KEEP == 1) then
  NFSB = NFSB+1
  do ISPART=1,NASPRT
    FSBARR(ISPART+(NASPRT+2)*(NFSB-1)) = ITRY(ISPART)
  end do
  FSBARR(1+NASPRT+(NASPRT+2)*(NFSB-1)) = NSDBLK
  FSBARR(1+NASPRT+1+(NASPRT+2)*(NFSB-1)) = NRDETS+1
  NRDETS = NRDETS+NSDBLK

end if
! Next walk: Lowest increasable switch value is at level lowest.
if (LOWEST > NASPRT) goto 600
do LEVUP=1,LOWEST-1
  SWITCH(+LEVUP) = 1
end do
SWITCH(+LOWEST) = 1+SWITCH(+LOWEST)
goto 500

600 continue
call mma_deallocate(SWITCH)
call mma_deallocate(VERTAB)
call mma_deallocate(DWNTAB)
call mma_deallocate(WEIGHT)
! Here we are. No more FS blocks can be constructed.
!TEST write(u6,*) ' Total nr of Slater determinants, NRDETS=',NRDETS
!TEST write(u6,*) ' Returning from VERTAB.'

end subroutine mkVERTAB
