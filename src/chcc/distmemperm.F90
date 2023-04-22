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

subroutine DistMemPerm(PosT)
! This routine does:
! define initial positions of Permanent
!
! PosT - initial and last position (I/O)

use Index_Functions, only: nTri_Elem
use chcc_global, only: intkey, no, nv, PosA, PosAex, PosFoo, PosFree, PosFvo, PosFvv, PosGoo, PosGvv, PosHoo, PosHvo, PosHvv, &
                       PosOE, PosT1n, PosT1o, printkey
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: PosT
integer(kind=iwp) :: length, nbas(1)

!1.1 Foo file
length = no*no
PosFoo = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Foo ',PosFoo,length

!1.2 Fvo file
length = no*nv
PosFvo = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Fvo ',PosFvo,length

!1.3 Fvv file
length = nv*nv
PosFvv = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Fvv ',PosFvv,length

!2 OE file

call Get_iArray('nBas',nBas,1)
length = nbas(1)

PosOE = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM OE ',PosOE,length

!3.1 T1o file
length = no*nv
PosT1o = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM T1o ',PosT1o,length

!3.2 T1n file
length = no*nv
PosT1n = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM T1n ',PosT1n,length

!4.1 Hoo file
length = no*no
PosHoo = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Hoo ',PosHoo,length

!4.2 Hvo file
length = no*nv
PosHvo = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Hvo ',PosHvo,length

!4.3 Hvv file
length = nv*nv
PosHvv = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Hvv ',PosHvv,length

!5.1 Goo file
length = no*no
PosGoo = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Goo ',PosGoo,length

!5.2 Hvv file
length = nv*nv
PosGvv = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Gvv ',PosGvv,length

!6.1 A files @@ A file medzi fixnymi je na zamyslenie (lebo je o4)
length = no*no*nTri_Elem(no)
PosA = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM A   ',PosA,length
if (intkey == 0) then
  PosAex = PosT
  PosT = PosT+length
  if (printkey >= 10) write(u6,*) 'DM Aex ',PosAex,length
else
  PosAex = PosT
end if

!omega PosFree - position of the space, where work arrays started
PosFree = PosT

return

end
