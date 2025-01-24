program make_eeps

  !local modules
  use iso_eep_support
  use eep
  use phases

  implicit none

  integer :: i, ierr, io, num
  character(len=8) :: version_string='unknown'
  character(len=file_path) :: input_file, history_columns_list
  character(len=file_path), allocatable :: history_files(:)
  type(track), pointer :: t=>NULL(), s=>NULL()
  logical :: do_phases = .true.
  logical :: do_he_star = .false.
  real(dp) :: initial_Y, initial_Z, Fe_div_H, v_div_vcrit, alpha_div_Fe,zsun

  namelist /eep_controls/ do_phases, center_gamma_limit, &
       center_carbon_limit, log_center_T_limit, high_mass_limit, &
       very_low_mass_limit, weight_center_rho_T_by_Xc, Teff_scale, &
       logL_scale, age_scale, Tc_scale, Rhoc_scale, make_bin_tracks, &
       eep_input_file, do_he_star

  ierr=0

  if(command_argument_count()<1) then
     write(*,*) '   make_eeps                  '
     write(*,*) '   usage: ./make_eeps <input> '
     stop       '   no command line argument   '
  endif

  !read input file, set up columns, eeps, format specs
!  call read_input(ierr)
  call read_input_list(ierr)
  if(ierr/=0) stop 'make_eeps: failed in read_input'

  ! allocate tracks for history files, read them in and convert to eep
  do i=1,num
     call alloc_track(history_files(i),t)
     call read_history_file(t,ierr)
     write(*,*) '-------------------------------------------'
     write(*,*) trim(t% filename), t% neep, t% ntrack, t% star_type
     !now set header info
     t% he_star = do_he_star
!     t% initial_Y = initial_Y
!     t% initial_Z = initial_Z
!     t% Fe_div_H  = Fe_div_H
     t% Fe_div_H = log10(t% initial_Z/zsun)
     t% alpha_div_Fe = 0.d0 !alpha_div_Fe
     t% v_div_vcrit = v_div_vcrit
     t% version_string = '1.1'
     !hack!
     t% merger = .false.
     !hack!
     if(ierr/=0) then
        write(0,*) 'make_eep: problem reading!'
        cycle
     endif
     call primary_eep(t)
     write(*,'(99i8)') t% eep
     if( all(t% eep == 0) ) then
        write(*,*) ' PROBLEM WITH TRACK: NO EEPS DEFINED '
     else
        call alloc_track(t% filename,s)
        call secondary_eep(t,s)
        print*, 'before GB check',s% eep
        s% he_star = t% he_star
        call check_for_bgb(s)
        if (s% has_BGB) print*, 'after GB check',s% eep
        if(do_phases) then
            if (s% has_BGB) then
                call set_track_phase_BGB(s)
            else
                call set_track_phase(s)
            endif
        endif
        
        s% filename = trim(eep_dir) // '/' // trim(s% filename) // '.eep'
        call write_track(s)
        deallocate(s)
        nullify(s)
     endif
     deallocate(t)
     nullify(t)
  enddo

contains

  subroutine read_input(ierr)
    integer, intent(out) :: ierr

    ierr=0

    version_string = adjustr(version_string)

    call get_command_argument(1,input_file)

    open(newunit=io,file=trim(input_file),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading ', trim(input_file)
       return
    endif
    read(io,*) !skip comment
    read(io,'(a8)') version_string
    read(io,*) !skip comment
    read(io,*) initial_Y, initial_Z, Fe_div_H, alpha_div_Fe, v_div_vcrit
    read(io,*) !skip comment
    read(io,'(a)') history_dir
    read(io,'(a)') eep_dir
    read(io,'(a)') iso_dir
    read(io,*) !skip comment
    read(io,'(a)') controls_file
    read(io,*) !skip comment
    read(io,'(a)') history_columns_list
    read(io,*) !skip comment
    read(io,*) num
    allocate(history_files(num))
    do i=1,num
       read(io,'(a)',iostat=ierr) history_files(i)
       if(ierr/=0) exit
    enddo
    close(io)
        
    open(newunit=io,file=trim(controls_file), action='read', status='old', iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading input.nml '
       return
    endif
    read(io, nml=eep_controls, iostat=ierr)
    close(io)

    !set number of secondary EEPs between each primary EEP
    call set_eep_interval(ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading eep input file'
       write(0,*) '            setting default EEPs     '
    endif
    !read history file format specs

    open(newunit=io,file='input.format',status='old',action='read',iostat=ierr)
    if(ierr/=0)then
       write(0,*) ' make_eeps: problem reading input.format'
       return
    endif
    read(io,*) 
    read(io,*) head
    read(io,*) main
    read(io,*) xtra
    close(io)
    !set up columns to be used
    call setup_columns(history_columns_list,ierr)
  end subroutine read_input


subroutine read_input_list(ierr)
    integer, intent(out) :: ierr

    ierr=0

    version_string = adjustr(version_string)

    call get_command_argument(1,input_file)
    print*, 'reading inputs from',trim(input_file)

    open(newunit=io,file=trim(input_file),status='old',action='read',iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading ', trim(input_file)
       return
    endif
    read(io,*) !skip comment
    read(io,*) zsun, v_div_vcrit
    read(io,*) !skip comment
    read(io,'(a)') history_dir
    read(io,'(a)') eep_dir
    read(io,*) !skip comment
    read(io,'(a)') controls_file
    close(io)
        
    ! get actual controls_file
    if (trim(controls_file) == 'hydrogen') then
        controls_file = 'input.nml'
    else
        controls_file = 'input_he.nml'
    endif
    
    open(newunit=io,file=trim(controls_file), action='read', status='old', iostat=ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading input.nml '
       return
    endif
    read(io, nml=eep_controls, iostat=ierr)
    close(io)

    ! get history tracks
    call get_files_from_path(history_dir,'.data',history_files,ierr)
    
    num = size(history_files)
    history_columns_list = 'inputs/my_history_columns.list'
    
    !set number of secondary EEPs between each primary EEP
    call set_eep_interval(ierr)
    if(ierr/=0) then
       write(0,*) ' make_eeps: problem reading eep input file'
       write(0,*) '            setting default EEPs     '
    endif
    !read history file format specs

    open(newunit=io,file='input.format',status='old',action='read',iostat=ierr)
    if(ierr/=0)then
       write(0,*) ' make_eeps: problem reading input.format'
       return
    endif
    read(io,*)
    read(io,*) head
    read(io,*) main
    read(io,*) xtra
    close(io)
    !set up columns to be used
    call setup_columns(history_columns_list,ierr)

  end subroutine read_input_list
  
  subroutine get_files_from_path(path,extension,file_list,ierr)
    character(LEN=strlen), intent(in) :: path
    character(LEN=*), intent(in) :: extension
    character(LEN=strlen), allocatable :: file_list(:)

    integer, intent (out) ::  ierr

    character(LEN=strlen) :: str,find_cmd
    integer :: n,i, io
    
        ierr = 0
        
!        find_cmd = 'find '//trim(path)//'/*'//trim(extension)//' -maxdepth 1 > .file_name.txt'
        find_cmd = 'ls -1 '//trim(path)//' > .file_name.txt'
        call system(find_cmd,ierr)
        
        if (ierr/=0) return

        io = alloc_iounit(ierr)
        open(io,FILE='.file_name.txt',action="read")

        !count the number of tracks
        n = 0
        do while(.true.)
            read(io,*,iostat=ierr)
            if(ierr/=0) exit
            n = n+1
        end do

        allocate(file_list(n))
        rewind(io)
        ierr = 0
        do i = 1,n
            read(io,'(a)',iostat=ierr)str
            if (ierr/=0) exit
            file_list(i) = trim(str)
        end do
        
        close(io)
        call free_iounit(io)
        
    end subroutine get_files_from_path
end program make_eeps
