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

subroutine DVer(IOpt,IC,L0,L,L2,DX,DY,DZ,Vert,Centr,Sphere,IntSph)

use PCM_Arrays, only: MxVert
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOpt, IC, L0, L, L2, IntSph(MxVert)
real(kind=wp), intent(out) :: DX, DY, DZ
real(kind=wp), intent(in) :: Vert(3,MxVert), Centr(3,MxVert), Sphere(4,*)
integer(kind=iwp) :: L1, NSJ
real(kind=wp) :: DNORM3, FACT, P(3), P1(3), P2(3), P3(3), PROD

! Trova la derivata della posizione del vertice L della tessera ITS
!
! IOpt = 0 : rispetto alla coordinata IC della sfera che, intersecando
!            ISPHE, forma il lato in esame. (ex derver)
! IOpt = 1 : rispetto al raggio della sfera aggiunta NSJ che,
!            intersecando la tessera ITS, crea il lato considerato.
!            (ex derver1, IC not referenced)
!
! Se stiamo considerando il primo vertice del lato : L > 0,
! il lato che si muove e' L-L2, e il vettore derivata e' diretto
! lungo la tangente al lato precedente (L0-L).
! Se stiamo considerando il secondo vertice del lato : L < 0,
! il lato che si muove e' L0-L, e il vettore derivata e' tangente
! al lato seguente (L-L2).

if (L > 0) then
  L1 = L
  NSJ = INTSPH(L1)
else
  L1 = -L
  NSJ = INTSPH(L0)
end if

! Il vettore P indica la posizione del vertice rispetto al centro
! della sfera NSJ
P(:) = VERT(:,L1)-Sphere(1:3,NSJ)

! Trova il versore tangente al lato L0-L1 se stiamo considerando il
! primo vertice, L2-L1 se consideriamo il secondo
if (L > 0) then
  P1(:) = VERT(:,L1)-CENTR(:,L0)
  P2(:) = VERT(:,L0)-CENTR(:,L0)
else
  P1(:) = VERT(:,L1)-CENTR(:,L1)
  P2(:) = VERT(:,L2)-CENTR(:,L1)
end if
call crprod(P1,P2,P3)
P2(:) = P3(:)
call crprod(P1,P2,P3)
DNORM3 = sqrt(P3(1)*P3(1)+P3(2)*P3(2)+P3(3)*P3(3))
P3(:) = P3(:)/DNORM3
! Trova la derivata
PROD = P(1)*P3(1)+P(2)*P3(2)+P(3)*P3(3)
if (IOpt == 0) then
  if ((PROD == ZERO) .and. (P(IC) /= ZERO)) then
    write(u6,'(a)') 'Stop in DVer.'
    call Abend()
  end if
  if (PROD == ZERO) PROD = One
  FACT = P(IC)/PROD
else if (IOpt == 1) then
  if (PROD == ZERO) then
    write(u6,'(a)') 'Stop in DVer.'
    call Abend()
  end if
  FACT = Sphere(4,NSJ)/PROD
else
  FACT = Zero
  write(u6,'(a)') 'Illegal IOpt in DVer.'
  call Abend()
end if
DX = FACT*P3(1)
DY = FACT*P3(2)
DZ = FACT*P3(3)

return

end subroutine DVer
