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

subroutine Tessera(IPRINT,MaxT,Nesf,NS,NV,XE,YE,ZE,RE,IntSph,PTS,CCC,PP,AREA)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), parameter :: MxVert = 20
integer(kind=iwp), intent(in) :: IPRINT, MaxT, Nesf, NS
integer(kind=iwp), intent(inout) :: NV
real(kind=wp), intent(in) :: XE(*), YE(*), ZE(*), RE(*)
integer(kind=iwp), intent(_OUT_) :: IntSph(MxVert,*)
real(kind=wp), intent(inout) :: PTS(3,MxVert)
real(kind=wp), intent(out) :: CCC(3,MxVert), PP(3), AREA
integer(kind=iwp) :: I, ICOP, ICUT, IND(MxVert), INTSCR(MxVerT), II, IV1, IV2, J, JJ, L, LTYP(MxVert), N, NSFE1
real(kind=wp) :: CCCP(3,MxVert), DE2, DELR, DELR2, DIST, DNORM, P1(3), P2(3), P3(3), P4(3), POINT(3), POINTL(3,MxVert), &
                 PSCR(3,MxVert), RC, RC2
real(kind=wp), parameter :: Small = 1.0e-12_wp, TOL = -1.0e-10_wp

! Coord. del centro che sottende l'arco tra i vertici
! n e n+1 (per i primi tre vertici e' sicuramente il centro della
! sfera) e sfera alla cui intersezione con NS appartiene l'arco (se
! appartiene alla sfera originaria INTSPH(N,MxTs)=NS)

AREA = Zero
do J=1,3
  CCC(1,J) = XE(NS)
  CCC(2,J) = YE(NS)
  CCC(3,J) = ZE(NS)
end do
! INTSPH viene riferito alla tessera MxTs, e in seguito riceve il
! numero corretto.
do N=1,3
  INTSPH(N,Maxt) = NS
end do
! Loop sulle altre sfere
do NSFE1=1,NESF
  if (NSFE1 == NS) cycle
  ! Memorizza i vertici e i centri che sottendono gli archi
  do J=1,NV
    INTSCR(J) = INTSPH(J,MaxT)
    do I=1,3
      PSCR(I,J) = PTS(I,J)
      CCCP(I,J) = CCC(I,J)
    end do
  end do
  ICOP = 0
  do J=1,MxVert
    IND(J) = 0
    LTYP(J) = 0
  end do
  ! Loop sui vertici della tessera considerata
  do I=1,NV
    DELR2 = (PTS(1,I)-XE(NSFE1))**2+(PTS(2,I)-YE(NSFE1))**2+(PTS(3,I)-ZE(NSFE1))**2
    DELR = sqrt(DELR2)
    if (DELR < (RE(NSFE1)-Small)) then
      IND(I) = 1
      ICOP = ICOP+1
    end if
  end do
  ! Se la tessera e' completamente coperta, la trascura
  if (ICOP == NV) return
  ! Controlla e classifica i lati della tessera: LTYP = 0 (coperto),
  ! 1 (tagliato con il II vertice coperto), 2 (tagliato con il I
  ! vertice coperto), 3 (bitagliato), 4 (libero)
  ! Loop sui lati
  do L=1,NV
    IV1 = L
    IV2 = L+1
    if (L == NV) IV2 = 1
    if ((IND(IV1) == 1) .and. (IND(IV2) == 1)) then
      LTYP(L) = 0
    elseif ((IND(IV1) == 0) .and. (IND(IV2) == 1)) then
      LTYP(L) = 1
    elseif ((IND(IV1) == 1) .and. (IND(IV2) == 0)) then
      LTYP(L) = 2
    elseif ((IND(IV1) == 0) .and. (IND(IV2) == 0)) then
      LTYP(L) = 4
      RC2 = (CCC(1,L)-PTS(1,L))**2+(CCC(2,L)-PTS(2,L))**2+(CCC(3,L)-PTS(3,L))**2
      RC = sqrt(RC2)
      ! Su ogni lato si definiscono 11 punti equispaziati, che vengono
      ! controllati
      do II=1,11
        POINT(1) = PTS(1,IV1)+II*(PTS(1,IV2)-PTS(1,IV1))/11
        POINT(2) = PTS(2,IV1)+II*(PTS(2,IV2)-PTS(2,IV1))/11
        POINT(3) = PTS(3,IV1)+II*(PTS(3,IV2)-PTS(3,IV1))/11
        POINT(1) = POINT(1)-CCC(1,L)
        POINT(2) = POINT(2)-CCC(2,L)
        POINT(3) = POINT(3)-CCC(3,L)
        DNORM = sqrt(POINT(1)**2+POINT(2)**2+POINT(3)**2)
        POINT(1) = POINT(1)*RC/DNORM+CCC(1,L)
        POINT(2) = POINT(2)*RC/DNORM+CCC(2,L)
        POINT(3) = POINT(3)*RC/DNORM+CCC(3,L)
        DIST = sqrt((POINT(1)-XE(NSFE1))**2+(POINT(2)-YE(NSFE1))**2+(POINT(3)-ZE(NSFE1))**2)
        if ((DIST-RE(NSFE1)) < TOL) then
          !if (DIST < RE(NSFE1)) then
          LTYP(L) = 3
          do JJ=1,3
            POINTL(JJ,L) = POINT(JJ)
          end do
          exit
        end if
      end do
    end if
  end do
  ! Se la tessera e' spezzata in due o piu' tronconi, la trascura
  ICUT = 0
  do L=1,NV
    if ((LTYP(L) == 1) .or. (LTYP(L) == 2)) ICUT = ICUT+1
    if (LTYP(L) == 3) ICUT = ICUT+2
  end do
  ICUT = ICUT/2
  if (ICUT > 1) return
  ! Creazione dei nuovi vertici e lati della tessera
  ! Loop sui lati
  N = 1
  do L=1,NV
    ! Se il lato L e' coperto:
    if (LTYP(L) == 0) cycle
    IV1 = L
    IV2 = L+1
    if (L == NV) IV2 = 1
    !*******************************************************************
    ! Se il lato L e' tagliato (con il I vertice scoperto):
    if (LTYP(L) == 1) then
      do JJ=1,3
        PTS(JJ,N) = PSCR(JJ,IV1)
        CCC(JJ,N) = CCCP(JJ,IV1)
      end do
      INTSPH(N,MaxT) = INTSCR(IV1)
      N = N+1
      ! Trova l'intersezione tra i due vertici del lato L
      !
      ! P1 = coord. del primo vertice
      ! P2 = coord. del secondo vertice
      ! P3 = coord. del centro dell'arco sotteso
      ! P4 = coord. dell'intersezione

      do JJ=1,3
        P1(JJ) = PSCR(JJ,IV1)
        P2(JJ) = PSCR(JJ,IV2)
        P3(JJ) = CCCP(JJ,IV1)
      end do

      call INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,0,IPRINT)
      ! Aggiorna i vertici della tessera e il centro dell'arco
      do JJ=1,3
        PTS(JJ,N) = P4(JJ)
      end do

      ! Il nuovo arco sara' sotteso tra questo e il prossimo punto
      ! di intersezione: il centro che lo sottende
      ! sara' il centro del cerchio di intersezione tra la sfera NS
      ! e la sfera NSFE1.

      DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+(ZE(NSFE1)-ZE(NS))**2
      CCC(1,N) = XE(NS)+(XE(NSFE1)-XE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(2,N) = YE(NS)+(YE(NSFE1)-YE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(3,N) = ZE(NS)+(ZE(NSFE1)-ZE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      INTSPH(N,MaxT) = NSFE1
      N = N+1
    end if
    !*******************************************************************
    ! Se il lato L e' tagliato (con il II vertice scoperto):
    if (LTYP(L) == 2) then
      ! Trova l'intersezione tra i due vertici del lato L
      !
      ! P1 = coord. del primo vertice
      ! P2 = coord. del secondo vertice
      ! P3 = coord. del centro dell'arco sotteso
      ! P4 = coord. dell'intersezione

      do JJ=1,3
        P1(JJ) = PSCR(JJ,IV1)
        P2(JJ) = PSCR(JJ,IV2)
        P3(JJ) = CCCP(JJ,IV1)
      end do

      call INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,1,IPRINT)
      ! Aggiorna i vertici della tessera e il centro dell'arco
      do JJ=1,3
        PTS(JJ,N) = P4(JJ)
        CCC(JJ,N) = CCCP(JJ,IV1)
      end do
      INTSPH(N,MaxT) = INTSCR(IV1)
      N = N+1
    end if
    !*******************************************************************
    ! Se il lato e' intersecato due volte:
    if (LTYP(L) == 3) then
      do JJ=1,3
        PTS(JJ,N) = PSCR(JJ,IV1)
        CCC(JJ,N) = CCCP(JJ,IV1)
      end do
      INTSPH(N,MaxT) = INTSCR(IV1)
      N = N+1
      ! Trova l'intersezione tra il primo vertice e un punto intermedio
      ! coperto
      !
      ! P1 = coord. del primo vertice
      ! P2 = coord. del secondo vertice
      ! P3 = coord. del centro dell'arco sotteso
      ! P4 = coord. dell'intersezione

      do JJ=1,3
        P1(JJ) = PSCR(JJ,IV1)
        P2(JJ) = POINTL(JJ,L)
        P3(JJ) = CCCP(JJ,IV1)
      end do
      call INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,0,IPRINT)
      ! Aggiorna i vertici della tessera e il centro dell'arco
      do JJ=1,3
        PTS(JJ,N) = P4(JJ)
      end do

      ! Il nuovo arco sara' sotteso tra questo e il prossimo punto
      ! di intersezione: il centro che lo sottende
      ! sara' il centro del cerchio di intersezione tra la sfera NS
      ! e la sfera NSFE1.

      DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+(ZE(NSFE1)-ZE(NS))**2
      CCC(1,N) = XE(NS)+(XE(NSFE1)-XE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(2,N) = YE(NS)+(YE(NSFE1)-YE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(3,N) = ZE(NS)+(ZE(NSFE1)-ZE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      INTSPH(N,MaxT) = NSFE1
      N = N+1

      ! Trova l'intersezione tra un punto intermedio coperto e il
      ! secondo vertice
      !
      ! P1 = coord. del primo vertice
      ! P2 = coord. del secondo vertice
      ! P3 = coord. del centro dell'arco sotteso
      ! P4 = coord. dell'intersezione

      do JJ=1,3
        P1(JJ) = POINTL(JJ,L)
        P2(JJ) = PSCR(JJ,IV2)
        P3(JJ) = CCCP(JJ,IV1)
      end do

      call INTER_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,NSFE1,1,IPRINT)
      ! Aggiorna il vertice e il centro dell'arco
      do JJ=1,3
        PTS(JJ,N) = P4(JJ)
        CCC(JJ,N) = CCCP(JJ,IV1)
      end do
      INTSPH(N,MaxT) = INTSCR(IV1)
      N = N+1
    end if

    !*******************************************************************
    ! Se il lato e' scoperto:
    if (LTYP(L) == 4) then
      do JJ=1,3
        PTS(JJ,N) = PSCR(JJ,IV1)
        CCC(JJ,N) = CCCP(JJ,IV1)
      end do
      INTSPH(N,MaxT) = INTSCR(IV1)
      N = N+1
    end if
    ! Controlla che il numero di vertici creati non sia eccessivo
    if (N > 11) then
      write(u6,'(/,a)') ' TESSERA: too many vertices in a tessera'
      call Abend()
    end if
  end do

  NV = N-1
end do

! Se la tessera non e' stata scartata, a questo punto ne troviamo
! l'area e il punto rappresentativo

call GAUBON(MaxT,XE,YE,ZE,RE,IntSph,NV,NS,PTS,CCC,PP,AREA,IPRINT)

return

end subroutine Tessera
