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
      SUBROUTINE rdinp_surfacehop()
      use Tully_variables
      IMPLICIT none
#include "Molcas.fh"
#include "surfacehop.fh"
      character*180 Get_Ln, Key, Line
      logical Found
      integer LuSpool, isfreeunit, NSTATE,i,j,ndata,maxHop
      real*8  temp(mxroot),AmatrixVR(mxroot*mxroot)
      real*8  AmatrixVI(mxroot*mxroot)
      external Get_Ln, isfreeunit

      LuSpool = isfreeunit(21)
      CALL SpoolInp(LuSpool)

      CALL RdNLst(LuSpool,'surfacehop')

      do
      Key = Get_Ln(LuSpool)
      Line = Key
      CALL UpCase(Line)

      select case(Line(1:4))
        case('TULL')
             tullyL=.true.
        case('DECO')
             Line = Get_Ln(LuSpool)
             CALL Get_F(1,DECO,1)
             decoherence=.true.
        case('SUBS')
             Line = Get_Ln(LuSpool)
             CALL Get_I(1,NSUBSTEPS,1)
        case('ETHR')
             Line = Get_Ln(LuSpool)
             CALL Get_F(1,Ethreshold,1)
        case('RTHR')
             Line = Get_Ln(LuSpool)
             CALL Get_F(1,RandThreshold,1)
        case('PSUB')
             tullySubVerb=.true.
        case('DMTX')
             Line = Get_Ln(LuSpool)
             CALL Get_I(1,NSTATE,1)
             do i=1, NSTATE
                Line = Get_Ln(LuSpool)
                CALL Get_F(1,temp,NSTATE)
                do j=1, NSTATE
                   AmatrixVR((i-1)*NSTATE+j)=temp(j)
                end do
             end do
             do i=1, NSTATE
                Line = Get_Ln(LuSpool)
                CALL Get_F(1,temp,NSTATE)
                do j=1, NSTATE
                   AmatrixVI((i-1)*NSTATE+j)=temp(j)
                end do
             end do
             call Qpg_zArray('AmatrixV', Found, ndata)
             IF (.NOT.Found) THEN
               call Put_dArray('RAmatrixV', AmatrixVR, NSTATE*NSTATE)
               call Put_dArray('IAmatrixV', AmatrixVI, NSTATE*NSTATE)
             else if (ndata.ne.NSTATE*NSTATE) then
               call warningmessage(2,'The A matrix dimension on the '//
     &         'input is different from the one in runfile')
               CALL Abend()
             end if
        case('FRAN')
             Line = Get_Ln(LuSpool)
             CALL Get_F(1,FixedRand,1)
             fixedrandL=.true.
        case('ISEE') ! initial seed number
             Line = Get_Ln(LuSpool)
             CALL Get_I(1,InitSeed,1)
             iseedL=.true.
        case('MAXH')
             Line = Get_Ln(LuSpool)
             CALL Get_I(1,maxHop,1)
             CALL Put_iScalar('MaxHopsTully',maxHop)
C             Write (6,*) 'MaxHops set to ', maxHop
        case('H5RE')
#ifdef _HDF5_
             lH5Restart = .True.
             Line = Get_Ln(LuSpool)
             CALL Get_S(1,FILE_H5RES,1)
#else
             write (6,*) 'The user asks to restart the dynamics '
             write (6,*) 'calculation from a HDF5 file, but this is not'
             write (6,*) ' supported in this installation.'
             call Quit_OnUserError()
#endif
        case('END ')
             exit
        case default
             Write (6,*) 'Unknown keyword: ', trim(Key)
             CALL Abend()
      end select
      end do
      CLOSE(LuSpool)
      RETURN
      END
