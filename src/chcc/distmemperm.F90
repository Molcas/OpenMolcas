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
        subroutine DistMemPerm (PossT)
!
!       This routine do:
!       define initial possitions of Permanent
!
!       PossT    - initial and last possition (I/O)
!
        implicit none
#include "chcc1.fh"
!
        integer PossT
!
!       help variables
        integer length,nbas(1)
!
!
!1.1    Foo file
        length=no*no
        PossFoo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Foo ',PossFoo,length
        end if
!
!1.2    Fvo file
        length=no*nv
        PossFvo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Fvo ',PossFvo,length
        end if
!
!1.3    Fvv file
        length=nv*nv
        PossFvv=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Fvv ',PossFvv,length
        end if
!
!
!2      OE file
!
        call Get_iArray('nBas',nBas,1)
        length=nbas(1)
!
        PossOE=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM OE ',PossOE,length
        end if
!
!
!3.1    T1o file
        length=no*nv
        PossT1o=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM T1o ',PossT1o,length
        end if
!
!3.2    T1n file
        length=no*nv
        PossT1n=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM T1n ',PossT1n,length
        end if
!
!
!4.1    Hoo file
        length=no*no
        PossHoo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Hoo ',PossHoo,length
        end if
!
!4.2    Hvo file
        length=no*nv
        PossHvo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Hvo ',PossHvo,length
        end if
!
!4.3    Hvv file
        length=nv*nv
        PossHvv=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Hvv ',PossHvv,length
        end if
!
!
!5.1    Goo file
        length=no*no
        PossGoo=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Goo ',PossGoo,length
        end if
!
!5.2    Hvv file
        length=nv*nv
        PossGvv=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM Gvv ',PossGvv,length
        end if
!
!6.1    A files @@ A file medzi fixnymi je na zamyslenie (lebo je o4)
        length=no*no*no*(no+1)/2
        PossA=possT
        PossT=PossT+length
        if (printkey.ge.10) then
        write (6,*) 'DM A   ',PossA,length
        end if
        if (intkey.eq.0) then
          PossAex=PossT
          PossT=PossT+length
          if (printkey.ge.10) then
          write (6,*) 'DM Aex ',PossAex,length
          end if
        else
          PossAex=PossT
        end if
!
!
!omega  PossFree - Possition of the space, where work arrays started
        PossFree=PossT
!
        return
        end
