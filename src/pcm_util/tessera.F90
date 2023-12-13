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

subroutine Tessera(IPRINT,Nesf,NS,NV,XE,YE,ZE,RE,IntSph,PTS,CCC,PP,AREA)

use PCM_Arrays, only: MxVert
use Constants, only: Zero, Two, Eleven
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IPRINT, Nesf, NS
integer(kind=iwp), intent(inout) :: NV
real(kind=wp), intent(in) :: XE(*), YE(*), ZE(*), RE(*)
integer(kind=iwp), intent(out) :: IntSph(MxVert)
real(kind=wp), intent(inout) :: PTS(3,MxVert)
real(kind=wp), intent(out) :: CCC(3,MxVert), PP(3), AREA
integer(kind=iwp) :: I, ICOP, ICUT, IND(MxVert), INTSCR(MxVert), II, IV1, IV2, L, LTYP(MxVert), N, NSFE1
real(kind=wp) :: CCCP(3,MxVert), DE2, DELR, DELR2, DIST, DNORM, P1(3), P2(3), P3(3), P4(3), POINT(3), POINTL(3,MxVert), &
                 PSCR(3,MxVert), RC, RC2
real(kind=wp), parameter :: Small = 1.0e-12_wp, TOL = -1.0e-10_wp

! Coord. del centro che sottende l'arco tra i vertici
! n e n+1 (per i primi tre vertici e' sicuramente il centro della
! sfera) e sfera alla cui intersezione con NS appartiene l'arco (se
! appartiene alla sfera originaria INTSPH(N)=NS)

AREA = Zero
CCC(1,1:3) = XE(NS)
CCC(2,1:3) = YE(NS)
CCC(3,1:3) = ZE(NS)
! INTSPH viene riferito alla tessera MxTs, e in seguito riceve il
! numero corretto.
INTSPH(1:3) = NS
! Loop sulle altre sfere
do NSFE1=1,NESF
  if (NSFE1 == NS) cycle
  ! Memorizza i vertici e i centri che sottendono gli archi
  INTSCR(1:NV) = INTSPH(1:NV)
  PSCR(:,1:NV) = PTS(:,1:NV)
  CCCP(:,1:NV) = CCC(:,1:NV)
  ICOP = 0
  IND(:) = 0
  LTYP(:) = 0
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
    else if ((IND(IV1) == 0) .and. (IND(IV2) == 1)) then
      LTYP(L) = 1
    else if ((IND(IV1) == 1) .and. (IND(IV2) == 0)) then
      LTYP(L) = 2
    else if ((IND(IV1) == 0) .and. (IND(IV2) == 0)) then
      LTYP(L) = 4
      RC2 = (CCC(1,L)-PTS(1,L))**2+(CCC(2,L)-PTS(2,L))**2+(CCC(3,L)-PTS(3,L))**2
      RC = sqrt(RC2)
      ! Su ogni lato si definiscono 11 punti equispaziati, che vengono
      ! controllati
      do II=1,11
        POINT(:) = PTS(:,IV1)+II*(PTS(:,IV2)-PTS(:,IV1))/Eleven
        POINT(:) = POINT(:)-CCC(:,L)
        DNORM = sqrt(POINT(1)**2+POINT(2)**2+POINT(3)**2)
        POINT(:) = POINT(:)*RC/DNORM+CCC(:,L)
        DIST = sqrt((POINT(1)-XE(NSFE1))**2+(POINT(2)-YE(NSFE1))**2+(POINT(3)-ZE(NSFE1))**2)
        if ((DIST-RE(NSFE1)) < TOL) then
          !if (DIST < RE(NSFE1)) then
          LTYP(L) = 3
          POINTL(:,L) = POINT(:)
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
      PTS(:,N) = PSCR(:,IV1)
      CCC(:,N) = CCCP(:,IV1)
      INTSPH(N) = INTSCR(IV1)
      N = N+1
      ! Trova l'intersezione tra i due vertici del lato L
      !
      ! P1 = coord. del primo vertice
      ! P2 = coord. del secondo vertice
      ! P3 = coord. del centro dell'arco sotteso
      ! P4 = coord. dell'intersezione

      P1(:) = PSCR(:,IV1)
      P2(:) = PSCR(:,IV2)
      P3(:) = CCCP(:,IV1)

      call INTER_PCM(XE(NSFE1),YE(NSFE1),ZE(NSFE1),RE(NSFE1),P1,P2,P3,P4,0,IPRINT)
      ! Aggiorna i vertici della tessera e il centro dell'arco
      PTS(:,N) = P4(:)

      ! Il nuovo arco sara' sotteso tra questo e il prossimo punto
      ! di intersezione: il centro che lo sottende
      ! sara' il centro del cerchio di intersezione tra la sfera NS
      ! e la sfera NSFE1.

      DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+(ZE(NSFE1)-ZE(NS))**2
      CCC(1,N) = XE(NS)+(XE(NSFE1)-XE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(2,N) = YE(NS)+(YE(NSFE1)-YE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(3,N) = ZE(NS)+(ZE(NSFE1)-ZE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      INTSPH(N) = NSFE1
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

      P1(:) = PSCR(:,IV1)
      P2(:) = PSCR(:,IV2)
      P3(:) = CCCP(:,IV1)

      call INTER_PCM(XE(NSFE1),YE(NSFE1),ZE(NSFE1),RE(NSFE1),P1,P2,P3,P4,1,IPRINT)
      ! Aggiorna i vertici della tessera e il centro dell'arco
      PTS(:,N) = P4(:)
      CCC(:,N) = CCCP(:,IV1)
      INTSPH(N) = INTSCR(IV1)
      N = N+1
    end if
    !*******************************************************************
    ! Se il lato e' intersecato due volte:
    if (LTYP(L) == 3) then
      PTS(:,N) = PSCR(:,IV1)
      CCC(:,N) = CCCP(:,IV1)
      INTSPH(N) = INTSCR(IV1)
      N = N+1
      ! Trova l'intersezione tra il primo vertice e un punto intermedio
      ! coperto
      !
      ! P1 = coord. del primo vertice
      ! P2 = coord. del secondo vertice
      ! P3 = coord. del centro dell'arco sotteso
      ! P4 = coord. dell'intersezione

      P1(:) = PSCR(:,IV1)
      P2(:) = POINTL(:,L)
      P3(:) = CCCP(:,IV1)
      call INTER_PCM(XE(NSFE1),YE(NSFE1),ZE(NSFE1),RE(NSFE1),P1,P2,P3,P4,0,IPRINT)
      ! Aggiorna i vertici della tessera e il centro dell'arco
      PTS(:,N) = P4(:)

      ! Il nuovo arco sara' sotteso tra questo e il prossimo punto
      ! di intersezione: il centro che lo sottende
      ! sara' il centro del cerchio di intersezione tra la sfera NS
      ! e la sfera NSFE1.

      DE2 = (XE(NSFE1)-XE(NS))**2+(YE(NSFE1)-YE(NS))**2+(ZE(NSFE1)-ZE(NS))**2
      CCC(1,N) = XE(NS)+(XE(NSFE1)-XE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(2,N) = YE(NS)+(YE(NSFE1)-YE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      CCC(3,N) = ZE(NS)+(ZE(NSFE1)-ZE(NS))*(RE(NS)**2-RE(NSFE1)**2+DE2)/(Two*DE2)
      INTSPH(N) = NSFE1
      N = N+1

      ! Trova l'intersezione tra un punto intermedio coperto e il
      ! secondo vertice
      !
      ! P1 = coord. del primo vertice
      ! P2 = coord. del secondo vertice
      ! P3 = coord. del centro dell'arco sotteso
      ! P4 = coord. dell'intersezione

      P1(:) = POINTL(:,L)
      P2(:) = PSCR(:,IV2)
      P3(:) = CCCP(:,IV1)

      call INTER_PCM(XE(NSFE1),YE(NSFE1),ZE(NSFE1),RE(NSFE1),P1,P2,P3,P4,1,IPRINT)
      ! Aggiorna il vertice e il centro dell'arco
      PTS(:,N) = P4(:)
      CCC(:,N) = CCCP(:,IV1)
      INTSPH(N) = INTSCR(IV1)
      N = N+1
    end if

    !*******************************************************************
    ! Se il lato e' scoperto:
    if (LTYP(L) == 4) then
      PTS(:,N) = PSCR(:,IV1)
      CCC(:,N) = CCCP(:,IV1)
      INTSPH(N) = INTSCR(IV1)
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

call GAUBON(XE,YE,ZE,RE,IntSph,NV,NS,PTS,CCC,PP,AREA,IPRINT)

return

end subroutine Tessera
