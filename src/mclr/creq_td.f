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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************
      SubRoutine creq_td(q,rint,G2,idsym)
!
!     Constructs the Q matrix
!
      use Constants, only: Zero
      use MCLR_Data, only: nDens2, ipMatBA, ipMO, nA
      use input_mclr, only: nSym,nAsh,ntAsh,nBas
      Implicit None
      Integer idSym
      Real*8 Q(nDens2),rint(*),G2(ntash,ntash,ntash,ntash)
      Integer iS, jS, kS, lS, ipS, ijS, iAsh, jAsh, kAsh, lAsh, ipQ, ipi
!
!      Q = (pj|kl)d
!       pi         ijkl
!
       Q(:)=Zero
       Do iS=1,nSym
        ipS=iEOr(is-1,idsym-1)+1
         if (nBas(ips).ne.0) Then

        Do jS=1,nsym
         ijS=iEOR(is-1,js-1)+1
         Do kS=1,nSym
          ls=iEOr(ijs-1,ks-1)+1
          Do iAsh=1,nAsh(is)
           Do jAsh=1,nAsh(js)
            Do kAsh=1,nAsh(ks)
             Do lAsh=1,nAsh(ls)
              ipQ=ipMatba(ips,is)+nBas(ips)*(iAsh-1)
              ipi=ipMO(js,ks,ls)+                                       &
     &         (nBas(ips)*(jAsh-1+nAsh(js)*                             &
     &          (kAsh-1+nAsh(ks)*(lAsh-1))))
              call daxpy_(nBas(ips),                                    &
     &                   G2(iAsh+nA(is),jAsh+nA(js),                    &
     &                      kAsh+nA(ks),lAsh+nA(ls)),                   &
     &                   rint(ipI),1,Q(ipQ),1)
             End Do
            End Do
           End Do
          End Do
         End Do
        End Do
        end if
       End Do
       end SubRoutine creq_td
