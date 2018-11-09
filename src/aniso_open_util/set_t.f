      Subroutine set_T( nT, nTempMagn, TINPUT, TempMagn, Tmin, Tmax,
     &                  chit_exp, Texp,
     &                  T, XTexp )

      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
c input:
      Integer, intent(in)       :: nT, nTempMagn
      Logical, intent(in)       :: TINPUT
      Real(kind=wp), intent(in) :: Tmin, Tmax, TempMagn(nTempMagn),
     &                             Texp(nT), chit_exp(nT)
      Real(kind=wp), intent(out):: T(nT+nTempMagn), XTexp(nT+nTempMagn)
c local variables:
      Integer :: i
      Real(kind=wp) :: dltt

      Call qEnter('set_T')

      ! set nT, T(i) and XTexp(i) arrays:
      T=0.0_wp
      XTexp=0.0_wp
!---------------------------------------------------------------------!
      If ( TINPUT ) Then
      ! case 1:  T(iT) computed from input values: Texp, and nTempMagn
      !          nTempMagn = 0
         If ( nTempMagn > 0  ) Then
            Do i=1, nTempMagn
               T(i)=TempMagn(i)
            End Do
            Do i= 1, nT
                   T(i+nTempMagn) =    Texp(i)
               XTexp(i+nTempMagn) =chit_exp(i)
            End Do
         Else
            Do i = 1, nT
               T(i)     =     Texp(i)
               XTexp(i) = chit_exp(i)
            End Do
         End If
      Else
      ! case 2:  T(iT) computed from input values: Tmin, Tmax, nT
      !          and nTempMagn
         dltt=0.0_wp
         dltt=(tmax-tmin)/(dble(nT-1))

         If ( nTempMagn > 0 ) Then
            Do i=1, nTempMagn
               T(i)=TempMagn(i)
            End Do

            T(1+nTempMagn)=0.0001_wp
            Do i = 2, nT
               T(i+nTempMagn) = Tmin + dltt*dble(i-1)
            End Do
         Else !compute_magnetization

            T(1)=0.0001_wp
            Do i = 2, nT
               T(i) = Tmin + dltt*dble(i-1)
            End Do
         End If
      End If !tinput

!---------------------------------------------------------------------!

      ! check for T=0 values and replace them
      ! with a definite nonzero value:
      Do i=1,nT+nTempMagn
        If( T(i) == 0.0_wp ) Then
            T(i) = 0.0001_wp
        End If
      End Do

      Call qExit('set_T')
      Return
      End subroutine set_T
