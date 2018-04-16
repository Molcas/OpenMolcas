************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Integer Function MyGetKey_nosupport (InUnit, What, IValue,
     &    RValue, SValue,
     &    N, IArray, RArray)
      Implicit Real*8 (A-H,O-Z)
      character SValue *(*)
      Dimension IArray(*), RArray(*)
      character What
      character KWord*120
      MyGetKey_nosupport=0
      iptr=1
      i=1
 1    Read(InUnit,'(A)',Err=20, end=20) KWord
      If (KWord(1:1).eq.'*'.or.KWord.eq.' ') Go To 1
      Call UpCase(KWord)
        if(What.eq.'I') then
              Read(KWord,*,Err=20) IValue
          else
           if(What.eq.'R') then
               Read(KWord,*,Err=20) RValue
             return
           else
           if(What.eq.'A') then
               Read(KWord,*,Err=20, End=40) (IArray(i),i=iptr,N)
           else

            if(What.eq.'D') then
               Read(KWord,*,Err=20, End=40) (RArray(i),i=iptr,N)
            else

            if(What.eq.'S') then
               Call NoBlanks(SValue, 120, KWord)
               goto 100
            else
             if(What.eq.'U') then
               Read(KWord,*,Err=2) IValue
               What='I'
               Goto 100
2              Read(KWord,*,Err=3) RValue
               What='R'
               Goto 100
3              Call NoBlanks(SValue, 120, KWord)
               What='S'
               Goto 100
             endif
            endif
           endif
          endif
         endif
        endif
100    return
40      iptr=i
        goto 1
20      MyGetKey_nosupport=1
        return
       end
        subroutine NoBlanks(out,n,in)
        character out*(*), in*(n)
        integer flag
        flag=-1
        j=0
        do 1 i=1,len(in)
         if(flag.eq.-1.and.in(i:i).eq.' ') goto 1
         flag=0
         if(i.le.len(in)-1) then
           if(in(i:i+1).eq.'  ') goto 1
         endif
         j=j+1
         out(j:j)=in(i:i)
1       continue
        out(j+1:)=' '
        return
        end
