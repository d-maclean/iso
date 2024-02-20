module eep

  !MESA modules
  use interp_1d_lib, only: interp_4pt_pm

  !local modules
  use iso_eep_support

  implicit none
integer:: RGB_loc, TAMS_loc
  logical, parameter :: eep_verbose = .false.

contains

  subroutine secondary_eep(t,s)
    ! takes a regular track (read from history) and creates an
    ! EEP track with a fixed number of regularly spaced points
    ! requires that the primary EEPs have already been found.
    type(track), intent(in) :: t
    type(track), intent(out) :: s
    integer :: i, j, k, num_p, num_s

    s% version_string = t% version_string
    s% initial_Y = t% initial_Y
    s% initial_Z = t% initial_Z
    s% Fe_div_H  = t% Fe_div_H
    s% alpha_div_Fe = t% alpha_div_Fe
    s% v_div_vcrit = t% v_div_vcrit
    s% ncol = t% ncol
    s% star_type = t% star_type
    !determine total number of EEPs in this track
    num_p = 0
    num_s = 0
    do j=1,primary
       if(t% eep(j) > 0) then
          num_p = num_p + 1
          if(j > 1) then
             num_s = num_s + eep_interval(j-1)
          endif
       endif
    enddo
    if(eep_verbose) write(*,*) ' track, number of p, s eeps = ', i, num_p, num_s
    ! identify the attributes of the new EEP track
    s% filename = t% filename
    s% initial_mass = t% initial_mass
    s% MESA_revision_number = t% MESA_revision_number
    s% ntrack = num_p + num_s
    s% neep = num_p
    s% ncol = t% ncol
    allocate(s% cols(s% ncol))
    s% cols = t% cols

    ! allocate and fill new track
    allocate(s% tr(s% ncol,s% ntrack))
    allocate(s% dist(s% ntrack))
    allocate(s% eep(s% neep))
    s% tr = 0d0
    s% dist = 0d0

    k = 1
    ! fill with primary and secondary EEPs
    do j=1,num_p
       if(eep_verbose .and. t% eep(j) > 0) write(*,*) t% tr(i_age,t% eep(j))
       s% tr(:,k) = t% tr(:,t% eep(j))
       s% dist(k) = t% dist(t% eep(j))
       s% eep(j) = k
       !interpolation on distance to fill secondary EEPs
       if(j < num_p) call eep_interpolate(t,j,k,s)
       if(j<num_p) k = k + eep_interval(j) + 1
    enddo
  end subroutine secondary_eep

  subroutine eep_interpolate(t,j,k,s)
    !t is the input track from MESA history file
    !s is the EEP-based track under construction
    !j is the primary EEP number 
    !k is the secondary EEP number
    !this routine takes a slice from t and interpolates into the
    !corresponding slice in s
    type(track), intent(in) :: t
    integer, intent(in) :: j, k
    type(track), intent(inout) :: s
    real(dp) :: alpha, beta, dist, delta, new_dist, dx, a(3), x(4), y(4)
    real(dp), pointer :: vec(:)
    integer :: i, jj, j0, j1, n, num_eep
    logical, parameter :: linear = .true.

    if(t% eep(j) == t% eep(j+1)) return

    num_eep = eep_interval(j)

    ! determine distance between primary eeps
    dist = t% dist(t% eep(j+1)) - t% dist(t% eep(j))

    ! interval = (total distance) / (number of secondary eeps)
    delta = dist / real(num_eep+1,kind=dp)

    j0 = t% eep(j)
    n = t% ntrack
    do i=1,num_eep
       new_dist = t% dist(t% eep(j)) + delta*real(i,kind=dp)
       allocate(vec(n))
       vec = t% dist
       j0=binary_search(n, vec, j0, new_dist)
       deallocate(vec); nullify(vec)
       j1=j0+1
       if(linear)then
          alpha = (t% dist(j1) - new_dist)/(t% dist(j1)-t% dist(j0))
          beta = 1d0-alpha
          s% tr(:,k+i) = alpha*t% tr(:,j0) + beta*t% tr(:,j1)
       else
          dx = new_dist - t% dist(j0)

          !for interpolation purposes, reset j0 if near the boundaries
          if(j0 > t% ntrack-2) j0 = t% ntrack-2
          if(j0 < 2) j0=2

          do jj=1,t% ncol
             x = t% dist(j0-1:j0+2)
             y = t% tr(jj,j0-1:j0+2)
             call interp_4pt_pm( x, y, a)
             s% tr(jj,k+i) = t% tr(jj,j0) + dx*(a(1) + dx*(a(2) + dx*a(3)))
          enddo
       endif
       s% dist(k+i) = new_dist
    enddo
  end subroutine eep_interpolate


  subroutine primary_eep(t)
    ! sets the locations of the primary EEPs in a track read from a history data file
    type(track), intent(inout) :: t
    integer :: ieep, inc
    t% EEP = 0 !initialize
    ieep=1
    inc=1
    if(t% he_star) then !only do this section if starting with a He star
        t% EEP(ieep) = ZAHB(t,1)
        print*, 'ZAHB',t% EEP(ieep)
        if(check(t,ieep)) return; ieep=ieep+1
        
        t% EEP(ieep) = TAHB(t,1d-4,t% EEP(ieep-1)+inc)
        print*, 'TAHB',t% EEP(ieep)
        if(check(t,ieep)) return; ieep=ieep+1
        
        if(t% star_type == star_high_mass)then
            t% EEP(ieep) = CarbonBurn(t,t% EEP(ieep-1)+inc)
            if (t% EEP(ieep)==0)t% EEP(ieep) = t% ntrack
            print*, 'cBurn', t% EEP(ieep)  ! end of HG

        elseif(t% star_type <= star_low_mass)then
            t% EEP(ieep) = TPAGB(t,t% EEP(ieep-1)+inc)
            center_gamma_limit = 7

            if (t% EEP(ieep)==0)t% EEP(ieep) = WDCS(t,t% EEP(ieep-1)+inc)
            print*, 'TPAGB',t% EEP(ieep)
            if(t% EEP(ieep)==t% ntrack) return; ieep=ieep+1

!            t% EEP(ieep) = PostAGB(t,t% EEP(ieep-1)+inc)
!            ; if(check(t,ieep)) return; ieep=ieep+1]
            ! above doesn't work for helium stars
            ! use t.mass_conv_envelope/(t.mass-t.co_core_mass) for giants
            
            !center_gamma_limit = 3
            !t% EEP(ieep) = WDCS(t,t% EEP(ieep-1)+inc); if(t% EEP(i)==t% ntrack) return; ieep=ieep+1
            
            center_gamma_limit = 9.0
            t% EEP(ieep) = WDCS(t, t% EEP(ieep-1)+inc)
            print*, 'PostAGB',t% EEP(ieep)
        endif
    else !normal H star 
        !t% EEP(ieep) = PreMS_Tc(t,5.0d0,1); if(check(t,ieep)) return; ieep=ieep+1
        t% EEP(ieep) = PreMS_fudge(1); if(check(t,ieep)) return; ieep=ieep+1
        t% EEP(ieep) = ZAMS(t,t% EEP(ieep-1)+inc)
        print*, trim(t% filename),t% EEP(ieep), 'zams'
        if(check(t,ieep)) return; ieep=ieep+1

        t% EEP(ieep) = TAMS(t,3.5d-1,t% EEP(ieep-1)+inc)
        print*, trim(t% filename),t% EEP(ieep), 'iams'
        if(check(t,ieep)) return; ieep=ieep+1

        t% EEP(ieep) = TAMS(t,1d-12,t% EEP(ieep-1)+inc)
        print*, trim(t% filename),t% EEP(ieep), 'tams'
        if(check(t,ieep)) return; ieep=ieep+1

        t% EEP(ieep) = RGBTip(t,t% EEP(ieep-1)+inc)
        print*, trim(t% filename),t% EEP(ieep), 'rgb_tip'
        if(check(t,ieep)) return; ieep=ieep+1

        t% EEP(ieep) = ZAHB(t,t% EEP(ieep-1)+inc)
        print*, trim(t% filename),t% EEP(ieep), 'zahb'
        if(check(t,ieep)) return; ieep=ieep+1

        t% EEP(ieep) = TAHB(t,1d-4,t% EEP(ieep-1)+inc)
        print*, trim(t% filename),t% EEP(ieep), 'tahb'
        if(check(t,ieep)) return; ieep=ieep+1

        if(t% star_type <= star_low_mass)then
            t% EEP(ieep) = TPAGB(t,t% EEP(ieep-1)+inc)
            print*, trim(t% filename),t% EEP(ieep), 'tpagb'
            if(check(t,ieep)) return; ieep=ieep+1
            
!           if (abs(t% EEP(ieep-1)-t% ntrack)<1000) force_last = .true.
!           if (force_last) then
!                t% EEP(ieep-1) = t% ntrack
!                print*, 'forcing TPAGB to be the last point'
!                force_last = .false.
!                return
!           endif
            t% EEP(ieep) = PostAGB(t,t% EEP(ieep-1)+inc)
            print*, trim(t% filename),t% EEP(ieep), 'postagb'
            if(check(t,ieep)) return; ieep=ieep+1

            ! t% EEP(ieep) = WDCS(t, t% EEP(ieep-1)+inc)
!            print*, trim(t% filename),t% EEP(ieep), 'wdcs'
        elseif(t% star_type == star_high_mass)then
            t% EEP(ieep) = CarbonBurn(t,t% EEP(ieep-1)+inc)
            print*, trim(t% filename),t% EEP(ieep),'cburn'
        endif
    endif
  end subroutine primary_eep
  
  logical function check(t,i)
    type(track), intent(in) :: t
    integer, intent(in) :: i
    check = (t% EEP(i)==0 .or. t% EEP(i)==t% ntrack)
  end function check

  function PreMS_Tc(t,logTc,guess) result(PreMS)
    type(track), intent(in) :: t
    real(dp), intent(in) :: logTc
    integer, intent(in) :: guess
    integer :: i, my_guess, PreMS
    real(dp) :: my_logTc
    PreMS = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    else
       my_guess = guess
    endif

    ! check for high mass stars, which may
    ! have hotter Tc than the input logTc
    ! in that case try the first value + 0.01
    if(t% tr(i_Tc,1) > logTc) then
       my_logTc = t% tr(i_Tc,1) + 0.01d0
    else
       my_logTc = logTc
    endif

    do i=my_guess,t% ntrack
       if(t% tr(i_Tc,i) > my_logTc) then
          PreMS = i
          return
       endif
    enddo

  end function PreMS_Tc

  !an alternative to constant central temperature is constant age
  function PreMS_age(t,logAge,guess) result(PreMS)
    type(track), intent(in) :: t
    real(dp), intent(in) :: logAge
    integer, intent(in) :: guess
    integer :: i, my_guess, PreMS
    PreMS = 0

    if(t% initial_mass > 3d1)then
       PreMS=1
       return
    endif

    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    else
       my_guess = guess
    endif

    do i=my_guess,t% ntrack
       if(log10(t% tr(i_age,i)) > logAge) then
          PreMS = i
          return
       endif
    enddo
  end function PreMS_age

  integer function preMS_fudge(guess)
    integer, intent(in) :: guess
    preMS_fudge = guess
  end function preMS_fudge
  
  integer function ZAMS(t,guess)
    type(track), intent(in) :: t
    real(dp) :: Xmax, Xmin, Xc, LH, Lmin !g_max, 
    real(dp), parameter :: Lfac = 9.99d-1 !Xfac = 9.8d-1
    integer, intent(in) :: guess
    integer :: i, my_guess, ZAMS1, ZAMS2, ZAMS3

    ZAMS = 0; ZAMS1 = 0; ZAMS2 = 0; ZAMS3 = 0

    Xmax = t% tr(i_Xc,1)
    Xmin = Xmax - 1.0d-3
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    else
       my_guess = guess
    endif

    !small amount of central H burned
    i = my_guess
    Xc = t% tr(i_Xc,i)
    do while(Xc > Xmin .and. i < t% ntrack )
       i = i+1
       Xc = t% tr(i_Xc,i)
    enddo

    ZAMS1 = i

!!$    !test of L_H/L_tot > some fraction
!!$    LH = pow10(t% tr(i_logLH,i))
!!$    Lmin = Lfac * pow10(t% tr(i_logL,i))
!!$    do while(LH > Lmin)
!!$       i = i-1
!!$       LH = pow10(t% tr(i_logLH,i))
!!$       Lmin = Lfac * pow10(t% tr(i_logL,i))
!!$    enddo
!!$
!!$    ZAMS2 = i

    ZAMS3=maxloc(t% tr(i_logg,1:ZAMS1),dim=1)

    !currently we choose to use ZAMS3 definition but can also
    !incorporate the ZAMS2 definition, or even average 1,2,&3

    ZAMS=ZAMS3

    !in case no H-burning occurs, take the location of the highest 
    !central temperature
    if(t% star_type == sub_stellar) ZAMS = maxloc(t% tr(i_Tc,:),dim=1)
  end function ZAMS

  integer function TAMS(t,Xmin,guess) !also IAMS with a higher Xmin
    type(track), intent(in) :: t
    real(dp), intent(in) :: Xmin
    integer, intent(in) :: guess
    integer :: i, my_guess
    real(dp), parameter :: max_age = 1.5d10 ! 15 Gyr
    TAMS = 0
    if(guess < 1) then 
       my_guess = 1
    elseif(guess >= t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i=my_guess,t% ntrack
       if(t% tr(i_Xc,i) < Xmin) then
          TAMS = i
          return
       endif
    enddo
    ! Xc test fails so consider instead the age for low mass tracks
    ! if the age of the last point > Max_Age, accept it
    if(t% initial_mass <= very_low_mass_limit .and. t% tr(i_age,t% ntrack) >= max_age) TAMS = t% ntrack
    if(t% merger) TAMS = t% ntrack
  end function TAMS

  integer function RGBTip(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    integer :: i, my_guess, RGBTip1, RGBTip2
    real(dp) :: Ymin, Lmax, Tmin
    RGBTip = 0; RGBTip1 = 0; RGBTip2 = 0
    if(guess < 1 .or. guess >= t% ntrack) then 
       return
    else
       my_guess = guess
    endif
    if(eep_verbose)then
       write(*,*) '    guess = ', guess
       write(*,*) ' my_guess = ', my_guess
       write(*,*) '   ntrack = ', t% ntrack
    endif
    Ymin = t% tr(i_Yc, my_guess) - 3d-2 !was 1d-2
    Lmax = -99d0
    Tmin = 6d0
    do i=my_guess, t% ntrack
       if(t% tr(i_Yc,i) > Ymin .and. t% tr(i_logL,i) > Lmax)then
          Lmax = t% tr(i_logL,i)
          RGBTip1 = i
       endif
       if(t% tr(i_Yc,i) > Ymin .and. t% tr(i_logTe,i) < Tmin)then
          Tmin = t% tr(i_logTe,i)
          RGBTip2 = i
       endif
    enddo
    RGBTip = min(RGBTip1,RGBTip2)
  end function RGBTip

  integer function ZAHB(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    integer :: i, my_guess, my_guess_2
    real(dp) :: Ymin, Tmin, LHemax
    ZAHB = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif
    Ymin = t% tr(i_Yc, my_guess) - 3d-2
    LHemax = -99d0
    my_guess_2 = my_guess

    do i=my_guess, t% ntrack
       if(t% tr(i_Yc,i) > Ymin .and. t% tr(i_logLHe,i) > LHemax)then
          LHemax = t% tr(i_logLHe,i)
          my_guess_2 = i
       endif
    enddo
    Tmin = 99d0
    do i=my_guess_2,t% ntrack
       if(t% tr(i_Yc,i) > Ymin .and. t% tr(i_Tc,i) < Tmin) then
          Tmin = t% tr(i_Tc,i)
          ZAHB = i
       endif
    enddo
!    if(ZAHB==guess) ZAHB=0
  end function ZAHB

  integer function TAHB(t,Ymin,guess)
    type(track), intent(in) :: t
    real(dp), intent(in) :: Ymin
    integer, intent(in) :: guess
    integer :: i, my_guess
    TAHB = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i=my_guess,t% ntrack
       if(t% tr(i_Yc,i) < Ymin)then
          TAHB = i
          return
       endif
    enddo
  end function TAHB

  integer function TPAGB(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    real(dp), parameter :: Ymin = 1d-6
    real(dp), parameter :: HeShellMin = 5d-2 !was 1d-1
    real(dp) :: Yc, HeShell
    integer :: i, my_guess
    TPAGB = 0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i=my_guess, t% ntrack
       Yc = t% tr(i_Yc,i)
       HeShell = t% tr(i_He_Core,i) - t% tr(i_CO_Core,i)
       if(Yc < Ymin .and. HeShell < HeShellMin)then
          TPAGB = i
          return
       endif
    enddo
  end function TPAGB

  integer function PostAGB(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    real(dp), parameter :: core_mass_frac_limit = 8d-1
    !real(dp), parameter :: log_Teff_limit = 3.6d0
    real(dp) :: Tc_now, TC_end, core_mass_frac
    integer :: i, my_guess
    PostAGB=0
    if(guess < 1 .or. guess > t% ntrack) then 
       my_guess = 1
    elseif(guess == t% ntrack)then
       return
    else
       my_guess = guess
    endif

    Tc_now = t% tr(i_Tc,my_guess)
    Tc_end = t% tr(i_Tc,t% ntrack)

    if(Tc_now > Tc_end)then
       !has TP-AGB
       do i=my_guess, t% ntrack
          core_mass_frac = t% tr(i_co_core,i) / t% tr(i_mass,i) 
          if(core_mass_frac > core_mass_frac_limit)then
             PostAGB = i
             return
          endif
       enddo
    else
       !has no TP-AGB
       return
    endif
    
  end function PostAGB

  integer function WDCS(t,guess)
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    integer :: my_guess, i
    WDCS=0
    if(guess < 1) then 
       my_guess = 1
    elseif(guess >= t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i = my_guess, t% ntrack
       if(t% tr(i_gamma,i) >= center_gamma_limit)then
          WDCS=i
          return
       endif
    enddo
  end function WDCS

  integer function CarbonBurn(t,guess) !for high-mass stars
    type(track), intent(in) :: t
    integer, intent(in) :: guess
    integer :: my_guess, i
    real(dp), parameter :: limit_XY=1d-2    !was 1d-8
    CarbonBurn = 0
    if(guess < 1)then
       my_guess = 1
    elseif(guess>=t% ntrack)then
       return
    else
       my_guess = guess
    endif
    do i=my_guess, t% ntrack
       if(t% tr(i_Xc,i) < limit_XY .and. t% tr(i_Yc,i) < limit_XY .and. t% tr(i_Cc,i) < center_carbon_limit)then
          CarbonBurn=i
          exit
       endif
    enddo
  end function CarbonBurn

 
    subroutine check_for_bgb(s)
    type(track) :: s,t
    integer :: j_bgb,j,l,k,temp_eep_interval(primary-1)
    integer, allocatable :: temp_eep(:)
    !check for BGB if it's RGB

          if (s% ntrack>=s% eep(5)) then
            j_bgb = BGB(s)
            if (j_bgb>0) then

                print*,trim(s% filename),' has BGB at ',j_bgb
                s% has_BGB = .true.
                l = s% neep+1
                allocate(temp_eep(l))

                do j = 1, l
                    if (j<=4) then !TAMS
                        temp_eep(j) = s% eep(j)
                    elseif(j==5) then
                        temp_eep(j) = 554
                    else
                        temp_eep(j) = s% eep(j-1)
                    endif
                enddo
                print*, 'new eeps', temp_eep
                temp_eep_interval = eep_interval
                eep_interval(4)= 99
                eep_interval(5)= 50
                k = 454
                t = s
                deallocate(t% eep)
                t% neep = 11
                allocate(t% eep(t% neep))
                t% eep = temp_eep
                t% eep(5) = s% eep(4)+j_bgb
                if (s% neep>5) t% eep(6) = s% eep(5)
                s% tr(:,554) = t% tr(:,t% eep(5))
                s% dist(554) = t% dist(t% eep(5))
                do j = 4,5
!                   s% eep(j) = k
                   !interpolation on distance to fill secondary EEPs
                   call eep_interpolate(t,j,k,s)
                   k = k + eep_interval(j) + 1
                enddo
                deallocate(s% eep)
                s% neep = l
                allocate(s% eep(l))
                s% eep = temp_eep
                deallocate(temp_eep)
                eep_interval = temp_eep_interval
            endif
          endif
    end subroutine check_for_bgb


    integer function BGB(t) result(j_bgb)
        type(track), intent(in) :: t
        real(dp) :: T_bgb_limit = 3.8

        integer :: peak, jfinal, jini, j_diff,k
        real(dp) :: mass,l_calc,diff
        real(dp), allocatable ::Lum(:),Teff(:),core_mass(:)
        real(dp), allocatable ::diff_L(:),diff_Te(:),dLdTe(:)

        j_bgb = -1

        jini = 454
        jfinal = 605
        j_diff = jfinal - jini
        if (j_diff<=0) return

        allocate(Lum(j_diff),Teff(j_diff),core_mass(j_diff))
        allocate (diff_L(j_diff -1),diff_Te(j_diff -1),dLdTe(j_diff -1))


        Teff = t% tr(i_logTe,jini:jfinal)
        if (Teff(j_diff)> T_bgb_limit) return


        mass = t% initial_mass
        Lum =  t% tr(i_logL,jini:jfinal)
        core_mass = t%tr(i_he_core,jini:jfinal)

        diff_L = Lum(2:j_diff)-Lum(1:j_diff-1)
        diff_Te = Teff(2:j_diff)-Teff(1:j_diff-1)
        dLdTe = diff_L/diff_Te

        peak = maxloc(dLdTe,dim = 1)
!        print*,"peak = ",mass,peak, Teff(peak),j_diff

        if(peak>=j_diff-1) return
        if (dLdTe(peak)>0.0 .and. dLdTe(peak+1)<0.0) then
        !checking for oscillations
            peak = maxloc(dLdTe(:peak-1),dim = 1)     !!whatif ends are very close??
        endif
        if (dLdTe(peak)>0.0) then          !convective core
            do k = peak,j_diff-1
                if (dLdTe(k)<1E-32 .and. Teff(k)<3.8) then
                    j_bgb = k
!                    t% has_RGB=.true.
                    return
                endif
            end do
        elseif (mass<2.0) then            !radiative core  !TODO: -- mhook?
            do k = 1,j_diff-1
                l_calc = log10(2.3E+5*(core_mass(k)**6))
                diff =  l_calc-Lum(k)
                if (abs(diff)<0.12) then
                    j_bgb = k
!                    t% has_RGB=.true.
                    return
                endif
            end do
        end if

        deallocate(Lum,Teff,core_mass)
        deallocate(diff_L,diff_Te,dLdTe)
    end function BGB
end module eep
