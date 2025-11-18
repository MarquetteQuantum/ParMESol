PROGRAM RK_Solution
! Solve the differential equation dC/dt = ..  using rk4 subroutine
	use barrier_energies_module
! Barrier energy module is used to get J,K-Dependent Barrier Energies
	use mpi
! Message Passing Interface(MPI) is used for parallelization	
	implicit none
! Parameters for Barrier Energies Module	
	integer, parameter :: barrier_rows=21
	integer, parameter :: barrier_columns = 65	
	real*8, dimension(barrier_rows, barrier_columns) :: barrier_energies_matrix
	character(len=200) :: barrier_directory, barrier_filename, barrier_filepath
! Standard Code Related Parameters
	character(len=200) :: directory, filename, filepath
	character(100) :: output, output1, output2, output3, output4, output5, output6, output7, output8, output9, output10, output11, &
						output14
	integer :: unit, aunit, bunit, cunit, dunit, funit, gunit, iunit, junit, kunit, munit, nunit, punit, qunit
	integer :: j, n, k, i, ios
! Propogation Parameters	
	real*8 :: t0, t_final, t, h, h_fs
	real*8 :: C_tot, dCdt_tot, dCdt_vdw, dCdt_cov, dCdt_tot_check, prev_krec, tolerance_low, tolerance_high, conv_precision
	real*8, dimension(:), allocatable:: C, dCdt, Cout, C_block, dCdt_block, C_equilibrium_block, C_check, dCdt_check, C_redo
	integer :: t_final_ps, conv, t_reset, conv_hit
! Parameters for Parallelization	
	integer :: pe_num, ierr_a, myid, chunk_size, variable_chunk_size, short_proc
	integer :: counter, a, w
	character(len=10) :: Js_string, Ks_string
	real*8 :: w_1
! Pressure related parameters	
	real*8 :: C_O_per_m3, C_O2_per_m3, M_per_m3, ref_pressure_per_m3, C_initial_O3_per_m3, pressure_std
! Tempertaure related parameters	
	real*8 :: temp_k, kt_energy_cm, kt_energy_j
! Conversions and Constants	
	real*8 :: h_j_per_s, hbar_js, c_cm_per_s, m_per_a0, j_per_k, j_per_cm, cm_per_k, pi
! Stabilization Rate Parameters	
	real*8 :: k0_m3_per_s, sigma0_m2, sigma_a2
! Window Parameters		
	real*8 :: lower_energy_lim, upper_energy_lim, upper_gamma_lim, lower_gamma_lim, localization_prob_cov, localization_prob_vdw, &
				G_filter
! Filtering Parameters
	real*8 :: E_ref
! Logical Statements	
	logical :: print_detail, truncation, K_dependent, print_transition_matrix, up_to_twenty_K_cells, temperature_dependence, &
				pressure_dependent_lindemann
	logical :: till_matrix, till_lindemann, K_blocks, transitions_only
! Printing Parameters	
	integer :: iteration_counter, print_freq
! Scanning Parameters	
	real*8 :: kij_min, bound_min, relative_error, relative_change, relative_error_prev
! Timing Parameters	
	real*8 :: start_time, sorting_start_time, sorting_end_time, start_time_matrix, end_time_matrix, start_time_master_eq, &
				end_time_master_eq, &
			  pre_prop_start_time, pre_prop_end_time, lind_start_time, lind_end_time, end_time
! Parameters for kinetics
	character(len=3) :: o3_molecule
    character(len=2) :: o2_molecule, o_atom
	character(len=200), dimension(:, :, :), allocatable :: filepath_JKsym
	integer :: Js, J_min, J_max, Ks, K_min, K_final, K_exact, K_max, Ks_indep, K_step
	integer :: vib_sym_well_start, vib_sym_well_last, vib_sym_well
	integer :: num_states, num_counter, broad_resonance_counter, narrow_resonance_counter
	integer :: J_rot_start, threshold_j
	integer, dimension(:), allocatable ::  K_value, J_value, sym_value, Resonance, state_tmp
	integer, dimension(:, :), allocatable :: state_matrix, state_matrix_tmp
	integer, dimension(:, :, :), allocatable :: num_states_J_K_sym, num_counter_J_K_sym, unit_J_K_sym, ios_J_K_sym
	real*8 :: k_rec
	real*8 :: dE_down
	real*8 :: dJ
	real*8 :: part_funcs_o2_per_m3, threshold_E
	real*8 :: Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
	real*8 :: tmp_energies, tmp_gammas, tmp_covalent_all, tmp_vdw_all, tmp_infinity_all, tmp_k_value, tmp_resonance, &
	tmp_sym_value, tmp_j_value, tmp_continium
	real*8 :: sum_eq_conc, K_eqs_tot
	real*8, dimension(2) :: dE
	real*8, dimension(:), allocatable :: Energies, Gammas, Covalent_All, Vdw_All, Infinity_All, Part_O3_cov, Part_O3_vdw, &
										Part_O3_inf, continium
	real*8, dimension(:), allocatable::  equilibrium_constants_m3, kdecs_per_s, channel_part_func, eq_concs, eq_conc_avg
	real*8, dimension(:, :), allocatable :: transition_matrix, Pfunc_matrix, krec_matrix
	real*8, dimension(:), allocatable :: First_term, Second_term, sum_rows, sum_direct_decay_rate, threshold_Energies_K, &
											Fourth_term, Fifth_term
	real*8, dimension(:, :, :), allocatable :: transition_matrix_3D
	real*8, dimension(:, :, :), allocatable :: threshold_Energies_J_K_sym
! Temperature Dependence at Low Pressure
	integer :: lunit
	character(100) :: output12
	real*8 :: QO3_tot, center, steepness, direct_decay_rate, steepness_factor
	real*8 :: dptmp, tmp_up, dpkt_energy_cm
	real*8, dimension(:), allocatable :: QO3
! Average E
	real*8 :: E_sum, E_avg, K_eqs_res
	real*8, dimension(:), allocatable :: E_bit, O_star
! Lindemann Limit	
	real*8 :: M_per_m3_Lindeman_first, M_per_m3_Lindeman_last, M_per_m3_Lindeman, krec3, krec3_low_p_lim, krec2_high_p_lim, &
	krec3_Lindeman, krec3_low, krec2_high_p_lim_ratio, weight, Lindeman_ratio, krec2_high, krec2_prime
	real*8, dimension(:), allocatable :: kstab, kstab_cov
	real*8, dimension(:, :), allocatable ::  krec2_Lindeman_high_p_lim_collected, krec3_Lindeman_collected, &
											krec3_Lindeman_low_p_lim_collected
! Idepth and Ilength
	integer :: ounit, trans_element_first, prevdepth, maxdepth, bound_counter, continium_counter
	character(100) :: output13
	real*8 :: ptran, matrix_el, procs_load_min, lindemann, statistic, E_lim, L_factor
	real*8, dimension(:), allocatable :: matrix_tmp, channel_factor
	integer, dimension(:), allocatable :: ilength, idepth
	real*8, dimension(:, :), allocatable ::trans_element, trans_element_tmp, lin_j_matrix, lin_matrix
	real*8, dimension(:, :), allocatable ::stat_j_matrix, stat_matrix
	real*8 :: dCdt_bound, dCdt_broad, dCdt_narrow, C_bound, C_narrow, C_broad, C_left, global_dis_energy, part_funcs_o2_per_m3_indp
!---------------------------------------------------------------------------------------------------------------------------	
	call MPI_INIT(ierr_a)
	call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr_a)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, pe_num, ierr_a)
	if (myid == 0) print *, pe_num, "Number of processors"	
!	write(myid*1000+1, *) 'I am processor #' , myid+1
		
	if (myid == 0) call cpu_time(start_time)
	
	barrier_directory = '/global/u2/n/nkr11/EQBMTOT'
	barrier_filename = 'barrier_energies.csv'
	barrier_filepath = trim(barrier_directory) // '/' // trim(barrier_filename)
	barrier_energies_matrix = read_matrix(barrier_filepath)
	
!   Print the matrix
!	if (myid == 0) then
!	do i = 1, barrier_rows
!			write(*, '(65(1x,E21.14))') barrier_energies_matrix(i, :)
!	end do
!	end if
	
! Conversion factors	
	m_per_a0 = 5.2917721092e-11
	j_per_k = 1.3806488e-23
	j_per_cm = 1.986445682933742e-23
	cm_per_k = j_per_k / j_per_cm
		
! Universal Constants
	h_j_per_s = 6.62607015e-34
	hbar_js = 1.054571726e-34
	c_cm_per_s = 2.9979245800e+10
	pi = 4.D0*DATAN(1.D0)
	
! Pressure related parameters
	pressure_std = 1.00d0
	ref_pressure_per_m3 = 6.44e+24
	steepness_factor = 6.d0
	
! Input Gas Mixture Parameters for reagents and bath gas	
	C_O_per_m3 = 6.44e+18
	C_O2_per_m3 = 6.44e+20
	M_per_m3 = pressure_std * ref_pressure_per_m3
	C_initial_O3_per_m3 = 0

! Temperature Parameters	
	temp_k = 298.d0
	kt_energy_j = temp_k * j_per_k;
	kt_energy_cm = temp_k * cm_per_k
	
! Sigma value in Bohr
	sigma_a2 = 97.088677d0
	sigma0_m2 = sigma_a2 * m_per_a0**2

! Parameters for energy transfer model in wavenumbers	
	dE_down = -240.5556d0
	dJ = 10.7894d0
		
! Energy Spectrum Parameters in wave numbers	
	lower_energy_lim = -10000.d0
	upper_energy_lim = 1000.d0
	upper_gamma_lim = 200.d0
	lower_gamma_lim = 1E-4
	localization_prob_cov = 1E-3
	localization_prob_vdw = 0.d0
	G_filter = 200.d0
	E_lim = kt_energy_cm
	L_factor = 1.d0
	
! Ozone isotope Labels
	o3_molecule = "666"
	o2_molecule = "66"
	o_atom = "6"

! Specify the directory and file name  
	filename = "state_properties.fwc"

	num_states = 0	
	
	J_min = 0
	J_max = 64
	K_exact = 20
	K_min = 0
	
	K_step = 1
	K_final = 20
	
	vib_sym_well_start = 0
	vib_sym_well_last = 1

! Lindemann Parameters
	M_per_m3_Lindeman = 6.44e17
	M_per_m3_Lindeman_last = 6.44e35
	
! Computing k0 = v(T)*sigma
	k0_m3_per_s = get_k0_2(o3_molecule, temp_k, sigma0_m2)

! Initial and final time and time step in picosecond and femtosecond units respectively
	if (M_per_m3 <= 0.001*6.44e+24) then
			h_fs = 50
	else if (M_per_m3 <= 0.01*6.44e+24) then
			h_fs = 25
	else if (M_per_m3 <= 0.1*6.44e+24) then
			h_fs = 10
	else if (M_per_m3 <= 1*6.44e+24) then
			h_fs = 5
	else if (M_per_m3 <= 100*6.44e+24) then
			h_fs = 1
	else if (M_per_m3 <= 1000*6.44e+24) then
			h_fs = 5e-1
	else if (M_per_m3 <= 10000*6.44e+24) then
			h_fs = 5e-2
	else if (M_per_m3 <= 100000*6.44e+24 ) then
			h_fs = 5e-3
	endif
	
	t0 = 0.0d0
	t_final = 1e-6	
	h = h_fs * 1e-15

! Convergance precision is determined for different circumstances	
	if (pressure_std > 10) then
		conv_precision = 1e-4
	else
		conv_precision = 1e-6
	end if	
	
	conv_hit = 250
	
! Physical Paramteres			
	K_dependent = .False.	
	up_to_twenty_K_cells = .False.
	kij_min = 1.d-3
	bound_min = 5e-1
	tolerance_low = 1e-6
	tolerance_high = 1.d0

! Printing Parameters
	print_freq = 1
	print_detail = .False.
	truncation = .False.
	print_transition_matrix = .False.
	temperature_dependence = .False.
	pressure_dependent_lindemann = .False.
	till_matrix = .False.
	till_lindemann = .False.
	K_blocks = .False.
	transitions_only = .False.
	
! Printing Initial Parameters into input file
	if (myid == 0) call cpu_time(sorting_start_time)
	
	if (myid == 0) then
		write(*, *) "Pressure :", pressure_std, "Std"
		write(*, *) "Temperature :", temp_k, "Kelvin"
		write(*, *) "Sigma :", sigma_a2, "Angstrom Square"
		write(*, *) "DE :", dE_down, "cm-1"
		write(*, *) "dJ :", dJ
		write(*, *) "Energy window :", upper_energy_lim, "cm-1 to", lower_energy_lim, "cm-1"
		write(*, *) "Gamma window :", upper_gamma_lim, "cm-1 to", lower_gamma_lim, "cm-1"
		write(*, *) "Covalent Probability above  :", localization_prob_cov
		write(*, *) "Van Der Wals Probability above  :", localization_prob_vdw
		write(*, *) "J values :", J_min, "to", J_max
		write(*, *) "K values :", K_min, "to", K_exact
		write(*, *) "Initial Timestep :", h, "seconds"
		write(*, *) "k0 is :", k0_m3_per_s, "m3/s"
		write(*, *)
	end if

! Allocating global dissociation minimum
	call get_threshold_energy_K(o3_molecule, o2_molecule, 1, threshold_E, threshold_j)
	global_dis_energy = threshold_E * j_per_cm
	
	if (myid == 0) write(*, *) "Global Dissociation Energy is :", global_dis_energy, "J or", threshold_E, "cm-1"
	
! All the data is read and filtered and the number of states in each J,K is counted		
	allocate(num_states_J_K_sym(J_max+1, K_final+1, vib_sym_well_last+1))
	num_states_J_K_sym = 0
	allocate(threshold_Energies_J_K_sym(J_max+1, K_final+1, vib_sym_well_last+1))
	threshold_Energies_J_K_sym = 0
	allocate(filepath_JKsym(J_max+1, K_final+1, vib_sym_well_last+1))
	
	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
		do Js = J_min, J_max
			write(Js_string, '(I0)') Js

! Opening all files for reading			
			if (Js < 20) then
				if (up_to_twenty_K_cells .eq. .True.) then
					K_max = Js
				else
					K_max = K_exact
				end if
			else
				if (up_to_twenty_K_cells .eq. .True.) then
					K_max = K_final
				else
					K_max = K_exact
				end if
			end if
			
			do Ks = K_min, K_max, K_step
			write(Ks_string, '(I0)') Ks
			directory = "/global/u2/n/nkr11/ozone_kinetics/data/resonances/mol_666/"				
			if (vib_sym_well == 0) then
			! vib_sym_well = 0
				if (mod(Ks, 2) == 0) then
				! K is even
				call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
				threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = threshold_E
!				write(*, *) Js, Ks, threshold_E
				! Reference energy is the highest value from the barrier energy and the threshold energy
				E_ref = max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1))
				directory = trim(directory) // "half_integers/J_" // trim(adjustl(Js_string)) // "/K_" &
				// trim(adjustl(Ks_string)) // "/symmetry_1"
				! Construct the full file path
				filepath = trim(directory) // '/' // trim(filename)
				filepath_JKsym(Js+1, Ks+1, vib_sym_well+1) = filepath					
					! Open and read the file with kinetics data
					open(newunit=unit, file=filepath, status='old', action='read', iostat=ios)
					! Skip the first line
					read(unit, *, iostat=ios)
					! Count the number of data points
						do
							read(unit, *, iostat=ios) Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
							if (ios /= 0) exit
								
							if (Energy_tmp > E_ref + lower_energy_lim .and. Energy_tmp < E_ref + upper_energy_lim &
							.and. Total_Gamma_tmp <= upper_gamma_lim .and. &
							Covalent_All_tmp > localization_prob_cov .and. Vdw_All_tmp > localization_prob_vdw) then
									
							num_states = num_states + 1
							num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1
							end if
						end do
						if (myid == 0) then
						write(*, *) "vib_sym_well = 0, Number of states with K = ", Ks, "is: ", num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)
						end if
						close(unit)						
				end if
			else	
				! vib_sym_well = 1 and ! K is odd
				if (mod(Ks, 2) /= 0) then	
				call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
				threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = threshold_E	
!				write(*, *) Js, Ks, threshold_E
				E_ref = max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1))
				directory = trim(directory) // "half_integers/J_" // trim(adjustl(Js_string)) // "/K_" &
				// trim(adjustl(Ks_string)) // "/symmetry_1"
				! Construct the full file path
				filepath = trim(directory) // '/' // trim(filename)
				filepath_JKsym(Js+1, Ks+1, vib_sym_well+1) = filepath
						
					! Open and read the file with kinetics data
					open(newunit=unit, file=filepath, status='old', action='read', iostat=ios)
					! Skip the first line
					read(unit, *, iostat=ios)
					! Count the number of data points
					do
						read(unit, *, iostat=ios) Energy_tmp, Total_Gamma_tmp, Covalent_All_tmp, Vdw_All_tmp, Infinity_All_tmp
						if (ios /= 0) exit
								
						if (Energy_tmp > E_ref + lower_energy_lim .and. Energy_tmp < E_ref + upper_energy_lim &
						.and. Total_Gamma_tmp <= upper_gamma_lim .and. &
						Covalent_All_tmp > localization_prob_cov .and. Vdw_All_tmp > localization_prob_vdw) then
									
						num_states = num_states + 1
						num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1									
						end if
					end do
					if (myid == 0) then
					write(*, *) "vib_sym_well = 1, Number of states with K = ", Ks, "is: ", num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)
					end if
					close(unit)
					end if
			end if
			end do
		end do
	end do
	 
	if (myid == 0) then
	write(*, *) "Total number of states", num_states
	end if
	
! Allocate arrays based on determined number of states and store filtered data
	allocate(Energies(num_states))
	allocate(Gammas(num_states))
	allocate(Covalent_All(num_states))
	allocate(Vdw_All(num_states))
	allocate(Infinity_All(num_states))
	allocate(J_value(num_states))
	allocate(K_value(num_states))
	allocate(sym_value(num_states))
	allocate(Resonance(num_states))
	allocate(continium(num_states))	
	
	Energies = 0
	Gammas = 0
	Covalent_All = 0
	Vdw_All = 0
	Infinity_All = 0
	J_value = 0
	K_value = 0
	sym_value = 0
	Resonance = 0
	bound_counter = 0
	narrow_resonance_counter = 0
	broad_resonance_counter = 0
	
! Temporary variables for sorting	
	tmp_energies = 0
	tmp_gammas = 0
	tmp_covalent_all = 0
	tmp_vdw_all = 0
	tmp_infinity_all = 0
	tmp_k_value = 0
	tmp_resonance = 0	
	tmp_continium = 0
	
	allocate(num_counter_J_K_sym(J_max+1, K_final+1, vib_sym_well_last+1))
	allocate(unit_J_K_sym(J_max+1, K_final+1, vib_sym_well_last+1))
	allocate(ios_J_K_sym(J_max+1, K_final+1, vib_sym_well_last+1))
	num_counter_J_K_sym = 1
	num_counter = 1
	
! Opening all files for reading
	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
	
		do Js = J_min, J_max
					
			if (Js < 20) then
				if (up_to_twenty_K_cells .eqv. .True.) then
					K_max = Js
				else
					K_max = K_exact
				end if
			else
				if (up_to_twenty_K_cells .eqv. .True.) then
					K_max = K_final
				else
					K_max = K_exact
				end if
			end if
			
			do Ks = K_min, K_max, K_step
				open(newunit=unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1), file=filepath_JKsym(Js+1, Ks+1, vib_sym_well+1), &
				status='old', iostat=ios_J_K_sym(Js+1, Ks+1, vib_sym_well+1))
				read(unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1), *, iostat=ios_J_K_sym(Js+1, Ks+1, vib_sym_well+1))			
			end do
		end do
	end do
	
! Read, sort and store filtered data
	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
	
		do Js = J_min, J_max
					
			if (Js < 20) then
				if (up_to_twenty_K_cells .eqv. .True.) then
					K_max = Js
				else
					K_max = K_exact
				end if
			else
				if (up_to_twenty_K_cells .eqv. .True.) then
					K_max = K_final
				else
					K_max = K_exact
				end if
			end if
			
! Read and store the filtered data
		do
		if (myid == 0) then
!		print*, 'entering most outer do loop (sorting step)'
		end if
			do Ks = K_min, K_final, K_step	
				do 
					if (num_counter_J_K_sym(Js+1, Ks+1, vib_sym_well+1) > num_states_J_K_sym(Js+1, Ks+1, vib_sym_well+1)) exit
					read(unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1), *, iostat=ios_J_K_sym(Js+1, Ks+1, vib_sym_well+1)) &
					Energies(num_counter), Gammas(num_counter), &
					Covalent_All(num_counter), Vdw_All(num_counter), Infinity_All(num_counter)
					
					E_ref = max(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1), barrier_energies_matrix(Ks+1, Js+1))
					if (Energies(num_counter) > E_ref + lower_energy_lim .and. Energies(num_counter) < E_ref + upper_energy_lim &
					.and. Gammas(num_counter) <= upper_gamma_lim .and. &
					Covalent_All(num_counter) > localization_prob_cov .and. Vdw_All(num_counter) > localization_prob_vdw) then					
! Neglect Gamma if state is below threshold and call it a bound state, otherwise call it a resonance
! Reference Enegeries that involve barrier energies can be used instead of threshold energies		
						if (Gammas(num_counter) < lower_gamma_lim) Gammas(num_counter) = 0.d0
						if (Energies(num_counter) < threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) ) then
							bound_counter = bound_counter + 1							
						else if (threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) < barrier_energies_matrix(Ks+1, Js+1) .and. &
						Energies(num_counter) > threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) .and. &
						Energies(num_counter) < barrier_energies_matrix(Ks+1, Js+1)) then
							Resonance(num_counter) = 1
							narrow_resonance_counter = narrow_resonance_counter + 1
						else if ((threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) > barrier_energies_matrix(Ks+1, Js+1) .and. &
						Energies(num_counter) > threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1)) .or. &
						(threshold_Energies_J_K_sym(Js+1, Ks+1, vib_sym_well+1) < barrier_energies_matrix(Ks+1, Js+1) .and. &
						Energies(num_counter) > barrier_energies_matrix(Ks+1, Js+1))) then
							broad_resonance_counter = broad_resonance_counter + 1
							Resonance(num_counter) = 2
						end if
						if (Gammas(num_counter) > G_filter) then
							continium(num_counter) = 1
							continium_counter = continium_counter + 1
						end if	
						J_value(num_counter) = Js
						K_value(num_counter) =  Ks
						sym_value(num_counter) = vib_sym_well
						if (num_counter > 1) then
							do k = num_counter, 2, -1
								if (Energies(k) < Energies(k-1)) then
									! Sorting Energies
									tmp_energies = Energies(k)
									Energies(k) = Energies(k-1)
									Energies(k-1) = tmp_energies
									! Re-assign Gammas
									tmp_gammas = Gammas(k)
									Gammas(k) = Gammas(k-1)
									Gammas(k-1) = tmp_gammas
									! Re-assign cov probs
									tmp_covalent_all = Covalent_All(k)
									Covalent_All(k) = Covalent_All(k-1)
									Covalent_All(k-1) = tmp_covalent_all
									! Re-assign vdw probs
									tmp_vdw_all = Vdw_All(k)
									Vdw_All(k) = Vdw_All(k-1)
									Vdw_All(k-1) = tmp_vdw_all
									! Re-assign inf probs
									tmp_infinity_all = Infinity_All(k)
									Infinity_All(k) = Infinity_All(k-1)
									Infinity_All(k-1) = tmp_infinity_all
									! Re-assign J_value array
									tmp_j_value = J_value(k)
									J_value(k) = J_value(k-1)
									J_value(k-1) = tmp_j_value
									! Re-assign K_value array
									tmp_k_value = K_value(k)
									K_value(k) = K_value(k-1)
									K_value(k-1) = tmp_k_value
									! Re-assign sym_value array
									tmp_sym_value = sym_value(k)
									sym_value(k) = sym_value(k-1)
									sym_value(k-1) = tmp_sym_value
									! Re-assign Bound/Resonance array
									tmp_resonance = Resonance(k)
									Resonance(k) = Resonance(k-1)
									Resonance(k-1) = tmp_resonance	
									! Re-assign Continium array
									tmp_continium = continium(k)
									continium(k) = continium(k-1)
									continium(k-1) = tmp_continium
								else
									exit
								end if						
							end do
						end if						
						num_counter_J_K_sym(Js+1, Ks+1, vib_sym_well+1) = num_counter_J_K_sym(Js+1, Ks+1, vib_sym_well+1) + 1
						num_counter = num_counter + 1
					end if
				end do
			end do
			if (num_counter > num_states) goto 10
			exit
 		end do
						
		end do
	end do
	
! Closing all files	
10	do vib_sym_well = vib_sym_well_start, vib_sym_well_last
		do Js = J_min, J_max

			if (Js < 20) then
				if (up_to_twenty_K_cells .eqv. .True.) then
					K_max = Js
				else
					K_max = K_exact
				end if
			else
				if (up_to_twenty_K_cells .eqv. .True.) then
					K_max = K_final
				else
					K_max = K_exact
				end if
			end if

			do Ks =  K_min, K_max, K_step
				close(unit_J_K_sym(Js+1, Ks+1, vib_sym_well+1))
			end do
		end do
	end do
	
	if (myid == 0) write(*, *) "Counter over number of states:", num_counter-1
	if (myid == 0) write(*, *) "Counter over number of bound states:", bound_counter
	if (myid == 0) write(*, *) "Counter over number of narrow resonances :", narrow_resonance_counter
	if (myid == 0) write(*, *) "Counter over number of broad resonances :", broad_resonance_counter
	if (myid == 0) write(*, *) "Counter over number of continium states :", continium_counter
	if (myid == 0) call cpu_time(sorting_end_time)
	if (myid == 0) write(*, *) "Time for the second reading and sorting:", sorting_end_time-sorting_start_time, "seconds"

! Determing chunk size in derivs for parallelization of the matrix

	chunk_size = int(num_states / pe_num) + 1
	if (mod(num_states, pe_num) == 0) then
		variable_chunk_size = chunk_size - 1
	else
		short_proc = (chunk_size * pe_num) - num_states
			if ((myid + 1) > (pe_num - short_proc)) then
				variable_chunk_size = chunk_size - 1
			else
				variable_chunk_size = chunk_size 
			end if	
	end if

! Check that sorting is correct	
	do i = 2, num_states
		if (Energies(i) < Energies(i-1)) then
			print *, 'error in sorting'
		end if
	end do

! K-Independent Partition Function 
	call get_threshold_energy_K(o3_molecule, o2_molecule, 0, threshold_E, threshold_j)
	J_rot_start = threshold_j
	part_funcs_o2_per_m3_indp = calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, 0)
	
! Printing partition functions of O+O2 system in K blocks	
	allocate(channel_part_func(K_final + 1))
	allocate(channel_factor(K_final + 1))
	channel_part_func = 0
	channel_factor = 0
	do Ks = K_min, K_final, K_step
		call get_threshold_energy_K(o3_molecule, o2_molecule, Ks, threshold_E, threshold_j)
		J_rot_start = threshold_j
		channel_part_func(Ks + 1) = calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, Ks)
		channel_factor(Ks + 1) = part_funcs_o2_per_m3_indp / channel_part_func(Ks + 1)
		if (myid == 0) then
			write(*, *) "Channael factor of K ", Ks, "is: ", channel_factor(Ks + 1)
		end if
	end do
	
! For each state compute O2 Partition function, Equilibrium constant and Decay rate
	allocate(equilibrium_constants_m3(num_states))
	allocate(kdecs_per_s(num_states))

	do  i = 1, num_states
		Js = J_value(i)
		Ks = K_value(i)
			if (K_dependent .eqv. .False.) then
				Ks_indep = 0
			else
				Ks_indep = Ks
			end if
		call get_threshold_energy_K(o3_molecule, o2_molecule, Ks_indep, threshold_E, threshold_j)
		J_rot_start = threshold_j
		part_funcs_o2_per_m3 = calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, Ks)
		equilibrium_constants_m3(i) = calculate_formation_decay_equilibrium(Energies(i), &
		temp_k, part_funcs_o2_per_m3, Js, Ks)
		kdecs_per_s(i) = channel_factor(K_value(i) + 1) * Gammas(i) * j_per_cm / hbar_js
	end do
	
! Sum of individual equilibrium constants to get the total equilibrium constant 
	allocate(Part_O3_cov(num_states))
	allocate(Part_O3_vdw(num_states))
	allocate(Part_O3_inf(num_states))
	allocate(pfunc_matrix(barrier_rows, barrier_columns))
	allocate(krec_matrix(barrier_rows, barrier_columns))
	
	K_eqs_tot = 0
	K_eqs_res = 0
	do  i = 1, num_states
		K_eqs_tot = K_eqs_tot + equilibrium_constants_m3(i)
		Part_O3_cov(i) = ((2 * J_value(i)) + 1) * Covalent_All(i) * exp(-(Energies(i)/kt_energy_cm))
		Part_O3_vdw(i) = ((2 * J_value(i)) + 1) * Vdw_All(i) * exp(-(Energies(i)/kt_energy_cm))
		Part_O3_inf(i) = ((2 * J_value(i)) + 1) * Infinity_All(i) * exp(-(Energies(i)/kt_energy_cm))
		Pfunc_matrix(K_value(i) + 1, J_value(i) + 1) = Pfunc_matrix(K_value(i) + 1, J_value(i) + 1) + Part_O3_cov(i)
		if(Resonance(i) /= 0) then
			K_eqs_res = K_eqs_res + equilibrium_constants_m3(i) * Covalent_All(i)
		end if	
	end do
	
	if (myid == 0) then
	write(*, *) " The total Resonance Equilibrium constant is", K_eqs_res, "m3"
	end if
	
! Writing spectrum properties of each state to the file		
	if (myid == 0) then
	output = 'spectrum_info.txt'
	open(newunit=iunit, file=output, status='replace')
	do j = 1, num_states
		write(iunit, &
		'((I8, 1x, E19.12, 1x, E19.12, 1x, I8, 1x, I8, 1x, I8, 1x, I8, 1x, E19.12, 1x, E19.12, 1x, E19.12, 1x, E19.12, 1x, E19.12))') j, &
			Energies(j), Gammas(j), J_value(j), K_value(j), sym_value(j), Resonance(j), continium(j), Covalent_All(j), &
!			(Energies(j) - threshold_Energies_J_K_sym(1, 1, 1)), (Energies(j)- threshold_Energies_J_K_sym(J_value(j) +1, K_value(j) +1, sym_value(j) +1)), &
!			(Energies(j) - barrier_energies_matrix(K_value(j) +1, J_value(j) +1))
			Vdw_All(j), Infinity_All(j), equilibrium_constants_m3(j)	
!		'((E19.12, 1x, E19.12, 1x, E19.12, 1x, E19.12, 1x))') Gammas(j), Part_O3_cov(j), Part_O3_vdw(j), Part_O3_inf(j)
			flush(iunit)
	end do
	close(iunit)
	end if
	
! Compute and print transition matrix, parallelization is inside the function
! Matrix truncation by finding ifinish and istart using two different criteria
	if (myid == 0) print*, 'Matrix calculations begin'
	if (myid == 0) call cpu_time(start_time_matrix)	
	
	allocate(matrix_tmp(num_states))
	allocate(state_tmp(num_states))
	allocate(idepth(chunk_size))
	allocate(sum_direct_decay_rate(chunk_size))

	idepth = 0
	sum_direct_decay_rate = 0
! ifinish is the last element which satisfies the criteria "kij_min" 
! (found going from last element of the row towards the diagonal)
!	do i = myid+1, num_states, pe_num !Last and First Row if condition
	do a = 1, variable_chunk_size
		i = (myid+1) + (a-1)*pe_num
! Make sure that this array is zero, before starting new row of the matrix
		matrix_tmp = 0.d0	
! istart is the first element which satisfies the criteria "kij_min* K_eq(i)/K_eq(j)" 
! (found going from first element of the row towards the diagonal)		
		do j = 1, i-1
		if(i == j) exit
! Transition Probability			
		ptran = Covalent_All(i) * Covalent_All(j) !+ Vdw_All(i) * Vdw_All(j) ! + Inf(i) * Inf(j)
		matrix_el = ptran * exp((Energies(i) - Energies(j)) / dE_down)&
			* exp(-abs(real(J_value(j))-real(J_value(i)))/dJ) * (equilibrium_constants_m3(i)/equilibrium_constants_m3(j))
			if (matrix_el .gt. kij_min*(equilibrium_constants_m3(i)/equilibrium_constants_m3(j))) then
				matrix_tmp(idepth(a) + 1) = matrix_el
				state_tmp(idepth(a) + 1) = j
				idepth(a) = idepth(a) + 1
			end if
		end do	
				
! Loop over columns	of matrix above the diagonal (going backwards)	
		do j = i+1, num_states
		if(i == j) exit
! Transition Probability			  
		ptran = Covalent_All(i) * Covalent_All(j) !+ Vdw_All(i) * Vdw_All(j) ! + Inf(i) * Inf(j)
		matrix_el = ptran * exp((Energies(j) - Energies(i)) / dE_down) * exp(-abs(real(J_value(j))-real(J_value(i)))/dJ)
			if (matrix_el .gt. kij_min) then
				matrix_tmp(idepth(a) + 1) = matrix_el
				state_tmp(idepth(a) + 1) = j
				idepth(a) = idepth(a) + 1			
			end if
		end do
		
! The initial matrix is allocated with maximum length				
		if (a == 1) then
			maxdepth = idepth(1)
			allocate(trans_element(maxdepth, chunk_size))
			allocate(state_matrix(maxdepth, chunk_size))
		end if
! Matrix is resized when longer than previous maximum length		
		if (idepth(a) > maxdepth) then
			prevdepth = maxdepth
			maxdepth = idepth(a)		
! Temporary matrix is allocated and previous data is stored		
			allocate(trans_element_tmp(prevdepth, chunk_size))
			allocate(state_matrix_tmp(prevdepth, chunk_size))
			trans_element_tmp = trans_element
			state_matrix_tmp = state_matrix
			deallocate(trans_element)
			deallocate(state_matrix)
! New matrix is allocated			
			allocate(trans_element(maxdepth, chunk_size))
			allocate(state_matrix(maxdepth, chunk_size))
				do k = 1, a-1
					do j = 1, prevdepth
					trans_element(j,k) = trans_element_tmp(j,k)
					state_matrix(j,k) = state_matrix_tmp(j,k)
					end do
				end do	
			deallocate(trans_element_tmp)
			deallocate(state_matrix_tmp)
		end if
		
! Elements are written into each matrix		
		do j = 1, idepth(a)
			trans_element(j, a) = matrix_tmp(j)
			state_matrix(j,a) = state_tmp(j)
		end do
		
! Sum of all Non-continium to continium states	
		do j = 1, idepth(a)
			if (continium(i) == 0 .and. continium(state_matrix(j, a)) == 1) then
				sum_direct_decay_rate(a) = sum_direct_decay_rate(a) + (trans_element(j, a) &
				*(equilibrium_constants_m3(state_tmp(j))/equilibrium_constants_m3(i)))		
			end if	
		end do	
	end do
	
	if (myid == 0) call cpu_time(end_time_matrix)
	if (myid == 0) write(*, *) "Time for calculating matricies :", &
	end_time_matrix-start_time_matrix, "seconds"	
	
! Load for each processor: collect the number of elements between istart and ifinish from each row of every prcessor
	procs_load_min = 0	
	do a = 1, variable_chunk_size
		procs_load_min = procs_load_min + idepth(a)
	end do
	
! Checking compression and memory saving for root processor
	if (myid == 0) then
	write(*, *) " The approximate length savings are", (1 - (real(maxdepth)/real(num_states))) * 100, "percent"
	write(*, *) " The approximate speed savings are", (1 - (real(procs_load_min)/real(chunk_size*num_states))) * 100, "percent"
	end if

! Stop if lindemann or Master Equation are not required
	if (till_matrix .eqv. .True.) stop
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Old Way of Calculating the transition matrix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rewrite upper part to calculate Ifinish and Istart

! Equilibrium concentrations of individual states
!	allocate(eq_concs(num_states))
!	do i = 1, num_states
!		eq_concs(i) = equilibrium_constants_m3(i) * C_O_per_m3 * C_O2_per_m3
!	end do	

!	allocate(transition_matrix(num_states, chunk_size))
!	transition_matrix = calculate_transition_matrix_unitless(Energies, Covalent_All, Vdw_All, Infinity_All, dE_down)

! Writing transition matrix elements to a file with row indices and column indices
!	if (myid == 0) then
!		if (print_transition_matrix .eqv. .True.) then
!			output2 = 'transition_matrix.csv'
!			open(newunit=kunit, file=output2, status='replace')
!			! The column index
!			do j = 1, size(transition_matrix, 1)
!				write(kunit, '(1x,I8,a1)', advance='no') j
!			end do
!				write (kunit, *)
!			
!			! The row index and matrix element if column=1, and only matrix element if column != 1
!			do i = 1, size(transition_matrix, 2)
!				do j = 1, size(transition_matrix, 1)
!				 if (j == 1) then
!					write(kunit, '(I8,a1,1x,E19.12,a1,1x)', advance='no') i, ',', transition_matrix(j, i)
!				 else 
!						write(kunit, '(E19.12,a1,1x)', advance='no') transition_matrix(j, i)
!				end if
!				end do
!				write (kunit, *)
!			end do
!			close(kunit)
!		end if
!	end if
	
! Build the overall transition matrix as a 3D matrix, then convert into 2D and print for debugging
!	if (myid .eqv. 0) then
!		allocate(transition_matrix_3D(num_states, chunk_size, pe_num))
!		transition_matrix_3D = 0
!	end if

! Gather local transition matrix from all processes to global transition matrix
!	call MPI_GATHER(transition_matrix, chunk_size*num_states, MPI_REAL8, &
!                  transition_matrix_3D(:, :, myid+1), chunk_size*num_states, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
!	call MPI_BARRIER( MPI_COMM_WORLD, ierr_a )				  

! Print global_matrix to a file on the master process
!	if (myid == 0) then
!	  call write_global_matrix_to_file(transition_matrix_3D, "output_global_matrix.csv", pe_num)
!	end if			  
!	if (myid .eqv. 0) deallocate(transition_matrix_3D)	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculating krec from Lindeman mechanism
	if (myid == 0) call cpu_time(lind_start_time)
	
	if (pressure_dependent_lindemann .eqv. .True.) then
	
	allocate(kstab(num_states))
	allocate(kstab_cov(num_states))
	allocate(krec3_Lindeman_collected(1, pe_num))
	allocate(krec3_Lindeman_low_p_lim_collected(1, pe_num))	
	allocate(krec2_Lindeman_high_p_lim_collected(1, pe_num))

	allocate(O_star(num_states))
	allocate(lin_j_matrix(barrier_rows, barrier_columns))
	allocate(lin_matrix(barrier_rows, barrier_columns))
	allocate(stat_j_matrix(barrier_rows, barrier_columns))
	allocate(stat_matrix(barrier_rows, barrier_columns))

! Pre Pressure Loop calculations
	kstab = 0.d0
	kstab_cov = 0.d0
	krec3_low_p_lim = 0.d0
	krec2_high_p_lim = 0.d0	
		do a = 1, variable_chunk_size
			i = (myid+1) + (a-1)*pe_num				
			Lindeman_ratio = 0.d0
			if (Resonance(i) == 1) then
				do j = 1, idepth(a)
					if (Resonance(state_matrix(j, a)) == 0) then
						kstab(i) = kstab(i) &
							     + trans_element(j, a) *(equilibrium_constants_m3(state_matrix(j, a))/equilibrium_constants_m3(i))
						kstab_cov(i) = kstab_cov(i) + Covalent_All(state_matrix(j, a)) &
									 * trans_element(j, a) *(equilibrium_constants_m3(state_matrix(j, a))/equilibrium_constants_m3(i))	 
						end if
					end do
				if (kstab(i) > 0) then
					Lindeman_ratio = kstab_cov(i) / kstab(i)
				end if
				krec3_low_p_lim = krec3_low_p_lim + equilibrium_constants_m3(i) * k0_m3_per_s* kstab_cov(i) 
				krec2_high_p_lim = krec2_high_p_lim + equilibrium_constants_m3(i) * kdecs_per_s(i)* Lindeman_ratio
				O_star(i) = ((2 * J_value(i)) + 1) * exp(-(Energies(i)/kt_energy_cm)) * Covalent_All(i)			 
				lin_j_matrix(K_value(i) + 1, J_value(i) + 1) = lin_j_matrix(K_value(i) + 1, J_value(i) + 1) + (O_star(i)/Covalent_All(i)) &
				* kstab_cov(i)	
				stat_j_matrix(K_value(i) + 1, J_value(i) + 1) = stat_j_matrix(K_value(i) + 1, J_value(i) + 1) + O_star(i) * kstab(i) 
				end if	
			end do

! Values of limits are gathered from all processors	
		call MPI_AllREDUCE(lin_j_matrix, lin_matrix, 1365, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr_a)
		call MPI_AllREDUCE(stat_j_matrix, stat_matrix, 1365, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr_a)
		call MPI_GATHER(krec3_low_p_lim, 1, MPI_REAL8, krec3_Lindeman_low_p_lim_collected, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)			
		call MPI_GATHER(krec2_high_p_lim, 1, MPI_REAL8, krec2_Lindeman_high_p_lim_collected, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
		call MPI_BARRIER( MPI_COMM_WORLD, ierr_a )

! Values of limits over all processors are summed up							
		if (myid == 0) then
			krec2_prime = 0.d0
			krec3_low = sum(krec3_Lindeman_low_p_lim_collected)
			krec2_high = sum(krec2_Lindeman_high_p_lim_collected)
			do i = 1, num_states
				if(Resonance(i) == 1) then
					krec2_prime = krec2_prime + equilibrium_constants_m3(i) * kdecs_per_s(i)
				end if
			end do	
		end if
			krec2_Lindeman_high_p_lim_collected = 0
			krec3_Lindeman_low_p_lim_collected = 0

! Values are made available to all processors				
		call MPI_BCAST(krec3_low, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
		call MPI_BCAST(krec2_high, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
		call MPI_BARRIER(MPI_COMM_WORLD, ierr_a)
		
! loop over pressure starts
		do while(M_per_m3_Lindeman <= M_per_m3_Lindeman_last) 
! krec is calculated by each processor
		krec3 = 0
		do a = 1, variable_chunk_size
			i = (myid+1) + (a-1)*pe_num		
			if (Resonance(i) == 1) then			
				weight = kdecs_per_s(i)/(kdecs_per_s(i) + M_per_m3_Lindeman * kstab(i)* k0_m3_per_s)
				krec3 = krec3 + equilibrium_constants_m3(i) * weight * k0_m3_per_s * kstab_cov(i)
			end if
		end do

! Values of krec are gathered from all processors		
		call MPI_GATHER(krec3, 1, MPI_REAL8, krec3_Lindeman_collected, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)			
		call MPI_BARRIER( MPI_COMM_WORLD, ierr_a )

! Values of krec over all processors are summed up							
		if (myid == 0) then
			krec3_Lindeman = sum(krec3_Lindeman_collected)
		end if
			krec3_Lindeman_collected = 0

! Values of krec are made available to all processors							
		call MPI_BCAST(krec3_Lindeman, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
		call MPI_BARRIER(MPI_COMM_WORLD, ierr_a)
		
! Results are printed into the file by root processor		
		if (myid == 0) then
		open(unit=1, file='LindemannPropogation.txt', position='append')
		write(1, '(*(E19.12,1x))') M_per_m3_Lindeman, krec3_Lindeman, krec3_low, krec2_high, krec2_prime
		close(1)
		end if		

! Pressure is increased and krec3 is reset for next itertion
		krec3_Lindeman = 0
		M_per_m3_Lindeman = M_per_m3_Lindeman * 10

! Printing everything
		output7 = 'Lindeman_matrix.txt'
		output8 = 'Statistic_matrix.txt'
		output9 = 'Pfunc_matrix.txt'
		open(newunit=punit, file=output7, status='replace')
		open(newunit=qunit, file=output8, status='replace')
		open(newunit=cunit, file=output9, status='replace')
		do i = 1, barrier_rows
			do j = 1, barrier_columns
				write(punit, '((E19.12, 1x))', advance = 'no') lin_matrix(i, j)
				write(qunit, '((E19.12, 1x))', advance = 'no') stat_matrix(i, j)
				write(cunit, '((E19.12, 1x))', advance = 'no') pfunc_matrix(i, j)
				lindemann = lindemann + ((lin_matrix(i, j) * k0_m3_per_s)/part_funcs_o2_per_m3)
				statistic = statistic + ((stat_matrix(i, j) * k0_m3_per_s)/part_funcs_o2_per_m3)
			end do
			write(punit, *)
			write(qunit, *)
			write(cunit, *)
		end do	
		close(punit)
		close(qunit)
		close(cunit)
	
		if (myid == 0) write(*, *) "Lindemann low pressure limit is:", lindemann	
		if (myid == 0) write(*, *) "Statisical low pressure limit is:", statistic
	
		end do		
	end if		

	if (myid == 0) call cpu_time(lind_end_time)
	if (myid == 0) write(*, *) "Time for Lindemann Calculation :", &
	lind_end_time-lind_start_time, "seconds"		
	
! Temperature Dependence of krec
	if (temperature_dependence .eqv. .True.) then
	
	allocate(QO3(num_states))
	allocate(E_bit(num_states))
	
	if (myid == 0) then
	output12 = 'lowp_tmpdp.txt'
	open(newunit=lunit, file=output12, status='replace')
	end if

! The partition function at different temperature is calculated	
	dptmp = 100
	tmp_up = ((log10(900.d0) - log10(dptmp))/19)
! Loop over temperatures	
	do while (dptmp <= 900) 
! Parameters arereset for each iteration	
	dpkt_energy_cm = dptmp * cm_per_k
	QO3 = 0
	QO3_tot = 0
	E_bit = 0
	E_sum = 0
	E_avg = 0
		do i = 1, num_states
			if (Resonance(i) == 1) then
				QO3(i) = Covalent_All(i) * (2 * (J_value(i) + 1)) * exp(-Energies(i)/dpkt_energy_cm)
				E_bit(i) = Energies(i) * QO3(i)
				QO3_tot = QO3_tot + QO3(i)
				E_sum = E_sum + E_bit(i)
			end if
		end do
	E_avg = E_sum/QO3_tot	
	if (myid == 0) write(lunit, '(*(E19.12,1x))') QO3_tot, E_avg
	dptmp = 10**(log10(dptmp) + tmp_up)
		
	end do
	close(lunit)
	end if
	
! Stop if master equation propogation is not required
	if (till_lindemann .eqv. .True.) stop
	
! Scaling of the unitless matrices
	trans_element = trans_element * k0_m3_per_s * M_per_m3
	
! Checking the reversibility principle and writing the result to a file
!	output11 = 'distribution_of_concentrations_over_kij_kji.txt'
!	open(newunit=gunit, file=output11, status='replace')
!	do i = 1, num_states
!		do j = 1, num_states
!			if (i /= j) then
!				write(gunit, '(I8, 1x, I8, 1x, E19.12, 1x, E19.12, 1x, E19.12, 1x, E19.12)') i, j, &
!				eq_concs(i)*transition_matrix(i, j), eq_concs(j)*transition_matrix(j, i), &
!				(eq_concs(i)*transition_matrix(i, j)-eq_concs(j)*transition_matrix(j, i))/&
!				((0.5)*(eq_concs(i)*transition_matrix(i, j)+eq_concs(j)*transition_matrix(j, i)))
!			end if
!		end do
!	end do
!	close(gunit)
	
!  Compute Constant terms for Master-equations
	if (myid == 0) call cpu_time(pre_prop_start_time)
		
	allocate(First_term(chunk_size))
	allocate(Second_term(chunk_size))
	allocate(Fourth_term(chunk_size))
	allocate(Fifth_term(chunk_size))
	allocate(sum_rows(chunk_size))
	
! Computing the "First term"
	First_term = 0
	Second_term = 0
	Fourth_term = 0
	sum_rows = 0
	do a = 1, variable_chunk_size
	i = (myid+1) + (a-1) * pe_num
		if (continium(i) == 0) then
			do j = 1, idepth(a)
			if (continium(state_matrix(j, a)) == 0) sum_rows(a) = sum_rows(a) + trans_element(j, a)&
			*(equilibrium_constants_m3(state_matrix(j, a))/equilibrium_constants_m3(i))	
			end do
			E_ref = max(threshold_Energies_J_K_sym(J_value(i)+1, K_value(i)+1, sym_value(i)+1),&
			barrier_energies_matrix(K_value(i)+1, J_value(i)+1))
			First_term(a) = kdecs_per_s(i)*equilibrium_constants_m3(i) * C_O_per_m3 * C_O2_per_m3
			Second_term(a) = kdecs_per_s(i) + sum_rows(a)
			center = ((E_ref + E_lim + (E_ref - E_lim))/2) + (L_factor*kt_energy_cm)
			steepness = 1.0986/(E_ref + E_lim - (L_factor*kt_energy_cm))
			direct_decay_rate = 0.5 * (tanh(steepness*(Energies(i) - center)) + 1) * (1 - Covalent_All(i))
			if (Energies(i) < (E_ref - E_lim)) then
				direct_decay_rate = 0.d0
			end if
			Fourth_term(a) = M_per_m3 * direct_decay_rate * k0_m3_per_s
			Fifth_term(a) = M_per_m3 * C_O_per_m3* C_O2_per_m3* direct_decay_rate * equilibrium_constants_m3(i) * k0_m3_per_s
!			write((myid+1)*1000, '(*(E19.12,1x))') real(i), Energies(i), First_term(a), Second_term(a), sum_direct_decay_rate(a), Fourth_term(a), Fifth_term(a)
		end if	
	end do

! Setting up and printing the initial conditions for concentrations at t=0
	allocate(C(num_states))
	allocate(dCdt(num_states))
	
	t = t0
	C = C_initial_O3_per_m3
	iteration_counter = 0

!   Uncomment if want to use equilibrium concentrations as initial ones	
!	do i = 1, num_states
!		C(i) = equilibrium_constants_m3(i)*C_O_per_m3*C_O2_per_m3
!	end do

! If Print detail is True concentration at every moment of time is calculated
	if (print_detail .eqv. .True.) then
		if (myid == 0) then
			output1 = 'propagated_concentration.txt'
			open(newunit=junit, file=output1, status='replace')
			write(junit, '(*(E19.12,1x))') t, C
		end if
	end if	

! Compute dCdt, total dCdt and krec at the initial moment of time
	call derivs(t, C, dCdt)
	C_tot = 0
	C_bound = 0 
	C_narrow = 0
	C_broad = 0
	C_left = 0
	dCdt_tot = 0
		do  i = 1, num_states
		if (continium(i) == 0) then
		   C_tot = C_tot + C(i)
		   dCdt_tot = dCdt_tot + (Covalent_All(i)) * dCdt(i)
			if (Resonance(i) == 0) then
				C_bound = C_bound + C(i)
			else if (Resonance(i) == 1) then
				C_narrow = C_narrow + C(i)
			else if (Resonance(i) == 2) then
				C_broad = C_broad + C(i)
			end if	
		end if	
		end do
	k_rec = dCdt_tot / (M_per_m3*(C_O_per_m3*C_O2_per_m3 - C_tot/K_eqs_tot))
	
! Compute sum of C, dCdt and C_eq  within each K block at the initial moment of time 
	if (K_blocks .eqv. .True.) then
	
	allocate(dCdt_block(K_final+1))
	allocate(C_block(K_final+1))
	allocate(C_equilibrium_block(K_final+1))
	dCdt_block = 0
	C_block = 0
	C_equilibrium_block = 0
	relative_error = 0
	relative_change = 0
	  	  
! Summation of dCdt, C and C_eq within each K block	 
! Change this part of the code into 2D for J,K blocks
	do i = 1, num_states
	if (continium(i) == 0) then
			dCdt_block(K_value(i) + 1) = dCdt_block(K_value(i) + 1) &
				+ (Covalent_All(i) + Vdw_All(i))*dCdt(i)
			C_block(K_value(i) + 1) = C_block(K_value(i) + 1) &
				+ C(i)
			C_equilibrium_block(K_value(i) + 1) = C_equilibrium_block(K_value(i) + 1) &
				+equilibrium_constants_m3(i) * C_O_per_m3 * C_O2_per_m3
	end if
	end do
	
		if(myid == 0) then
			output6 = 'equilibrium_concentrations_in_K_blocks'
				open(newunit=bunit, file=output6, status='replace')
				write(bunit, '(*(E19.12,1x))') h, 1.d0, C_equilibrium_block
				write(bunit, '(*(E19.12,1x))') t_final, 1.d0, channel_part_func
				close(bunit)
	
		end if
	end if	
! Writing the different output files		
			if (myid == 0) then
				output3 = 'recombination_coefficient_and_dCdt_tot.txt'
				open(newunit=nunit, file=output3, status='replace')
				write(nunit, '(*(E19.12,1x))') t, h, relative_error, C_tot, C_bound, C_narrow, C_broad, &
				dCdt_tot, k_rec, relative_change
			
				if (print_detail .eqv. .True.) then
					output4 = 'total_concentrations.txt'
					open(newunit=munit, file=output4, status='replace')
					write(munit, '(*(E19.12,1x))') t, C_tot
				
					output5 = 'equilibrium_concentrations.txt'
					open(newunit=aunit, file=output5, status='replace')
					write(aunit, '(*(E19.12,1x))') h, equilibrium_constants_m3 * C_O_per_m3 * C_O2_per_m3
					write(aunit, '(*(E19.12,1x))') t_final, equilibrium_constants_m3 * C_O_per_m3 * C_O2_per_m3
					close(aunit)
				end if
	
! Final equilibrium distribution of concentrations in a test when formation and decay are disabled
! And the initial concentration of each state is chosen to be C_initial_O3_per_m3
		if (transitions_only .eqv. .True. ) then
		
				output10 = 'equilibrium_concentrations_transitions_only.txt'
				open(newunit=funit, file=output10, status='replace')
				sum_eq_conc = 0
				
				do i = 1, num_states
					if (continium(i) == 0) sum_eq_conc = sum_eq_conc + exp(-Energies(i)/kt_energy_cm)
				end do
				
				allocate(eq_conc_avg(num_states))
				do i = 1, num_states
					if (continium(i) == 0) eq_conc_avg(i) = num_states * C_initial_O3_per_m3 * exp(-Energies(i)/kt_energy_cm) / sum_eq_conc
				end do
				write(funit, '(*(E19.12,1x))') h, eq_conc_avg
				write(funit, '(*(E19.12,1x))') t_final, eq_conc_avg
				close(funit)		
		end if	

! Parameters before propogation loop timing		
				if (myid == 0) call cpu_time(pre_prop_end_time)
				write(*, *) "All parameters before main propagation loop:", pre_prop_end_time-pre_prop_start_time, "seconds"	
			end if
	
	if (myid == 0) call cpu_time(start_time_master_eq)
		
!  The main time-propogation loop
	allocate(Cout(num_states))
	allocate(C_check(num_states))
	allocate(dCdt_check(num_states))
	allocate(C_redo(num_states))
	conv = 0
	t_reset = 0
	C_check = 0.d0
	dCdt_check = 0.d0
	C_redo = 0.d0
	prev_krec = 0.d0
	relative_error_prev = 0.d0
	do while (t<=t_final)
!			  if (iteration_counter == 100) stop 		  
			  C_redo = C
! Solve the differential equation using rk4 subroutine and update variables
!			  call rk4(C, dCdt, num_states, t, h, Cout, derivs)
!			  C = Cout 	 
			  
! Euler method is used her instead of Rungi_Kutta
! Two steps are done with smaller timestep and one step is done with larger timestep		    
			  do i = 1, num_states
			  if (continium(i) == 0) then
			  C_check(i) = C(i) + dCdt(i) * h * 2
			  C(i) = C(i) + dCdt(i) * h
			  end if
			  end do
			  
			  t = t + h
			  iteration_counter = iteration_counter + 1
			  call derivs(t, C, dCdt)
			  
			  do i = 1, num_states
			  if (continium(i) == 0) C(i) = C(i) + dCdt(i) * h
			  end do
			  
			  t = t + h
			  iteration_counter = iteration_counter + 1
			  call derivs(t, C, dCdt)
			  call derivs(t, C_check, dCdt_check)

! Values of krec and error at previous iteration are taken	  
			  prev_krec = k_rec	 
			  relative_error_prev = relative_error
! Processing and printing of the results: Individual and/or Total dCdt and Krec
			if (mod(iteration_counter, print_freq) == 0) then  
				  C_tot = 0
                  C_bound = 0 
                  C_narrow = 0
                  C_broad = 0
				  C_left = 0
                  dCdt_tot = 0
				  dCdt_tot_check = 0
                  do  i = 1, num_states
				  if (continium(i) == 0) then
					C_tot = C_tot + C(i)
                  	dCdt_tot = dCdt_tot + (Covalent_All(i)) * dCdt(i)
					dCdt_tot_check = dCdt_tot_check + (Covalent_All(i)) * dCdt_check(i)
                  	if (Resonance(i) == 0) then
						C_bound = C_bound + C(i)
                  	else if (Resonance(i) == 1) then
                  		C_narrow = C_narrow + C(i)
                  	else if (Resonance(i) == 2) then
                  		C_broad = C_broad + C(i)
                  	end if
				  end if				
                  end do	
                  
! Relative error is calculated				 
			relative_error = abs(dCdt_tot_check-dCdt_tot)/dCdt_tot

! Iteration is skipped if Error is higher than tolerance or if breakdown d=begins to occur			
			if (relative_error > tolerance_high .or. (relative_error > relative_error_prev & 
			.and. relative_error_prev > tolerance_low .and. relative_error_prev /= 0)) go to 68
				
				  k_rec = dCdt_tot / (M_per_m3 * (C_O_per_m3*C_O2_per_m3 - C_tot/K_eqs_tot))

! Relative change is calculated to check convergance				  
				  relative_change = abs(prev_krec - k_rec)/k_rec
				  
				  if (myid == 0) write(nunit, '(*(E19.12,1x))') t, h, relative_error, C_tot, C_bound, C_narrow, C_broad, &
					dCdt_tot, k_rec, relative_change, real(conv)

				  if (print_detail .eqv. .True.) then
					if (myid == 0) write(junit, '(*(E19.12,1x))') t, C
				  end if

! If error for different timesteps are low timestep is increased				  
				if (relative_error < tolerance_low) h = h * 2
! If krec change is low for sustained period converagnce is said to have been reached and propogation stopped				
				if (relative_change < conv_precision) conv = conv + 1
				if (conv == conv_hit) goto 69
				
			end if		  
			
! Contribution of individual J,K blocks into the rate coefficient
!				  dCdt_block = 0
!			      C_block = 0			 				  
!				  do i = 1, num_states
!					 dCdt_block(K_value(i) + 1) = dCdt_block(K_value(i) + 1) + (Covalent_All(i) + Vdw_All(i)) * dCdt(i)
!					 C_block(K_value(i) + 1) = C_block(K_value(i) + 1) + C(i)
!				  end do
				  				  
!				  if (myid == 0) write(munit, '(*(E19.12,1x))') t, C_block, dCdt_block
!				  flush(nunit)
!				  flush(munit)

! If the error is high or breakdown begins we reset to conditions before this step and continue with lower timestep
68			if (relative_error > tolerance_high .or. (relative_error > relative_error_prev & 
			.and. relative_error_prev > tolerance_low .and. relative_error_prev /= 0)) then
				t = t - (2 * h)
				iteration_counter = iteration_counter - 2
				C = C_redo
				relative_error = relative_error_prev
				h = h / 2
			! If multiple zig-zags occur lower tolerance is reduced	
				t_reset = t_reset + 1
				call derivs(t, C, dCdt)
				if (t_reset == 10) then
					tolerance_low = tolerance_low / 10
!					t_reset = 0
				end if	
			end if	
			
	end do

! Printing of the final results at the final moment of time				 
				C_tot = 0
				C_bound = 0 
				C_narrow = 0
				C_broad = 0
				C_left = 0
				dCdt_tot = 0
				do  i = 1, num_states
				if (continium(i) == 0) then
					C_tot = C_tot + C(i)
					dCdt_tot = dCdt_tot + (Covalent_All(i)) * dCdt(i)
					krec_matrix(K_value(i) + 1, J_value(i) + 1) = krec_matrix(K_value(i) + 1, J_value(i) + 1) + &
					dCdt(i) / (M_per_m3 * (C_O_per_m3*C_O2_per_m3 - C_tot/K_eqs_tot))
					if (Resonance(i) == 0) then
						C_bound = C_bound + C(i)
					else if (Resonance(i) == 1) then
						C_narrow = C_narrow + C(i)
					else if (Resonance(i) == 2) then
						C_broad = C_broad + C(i)
					end if
				end if		
				end do
			  
				  k_rec = dCdt_tot / (M_per_m3 * (C_O_per_m3*C_O2_per_m3 - C_tot/K_eqs_tot))
				  
69				  if (myid == 0) then
				  write(nunit, '(*(E19.12,1x))') t, h, relative_error, C_tot, C_bound, C_narrow, C_broad, &
					dCdt_tot, k_rec, relative_change, real(conv)
				  write(*, *) "The final concentration of Ozone is", C_tot, "m-3"
				  write(*, *) "The rate coefficient after equilibriation is", k_rec, "m6/s"
				  end if
! Print details
				  if (print_detail .eqv. .True.) then
					if (myid == 0) write(junit, '(*(E19.12,1x))') t, C
				  end if
				  
	close(junit)
	close(nunit)
	close(munit)

! Print Krec Contribution
	output14 = 'krec_matrix.txt'
	open(newunit=dunit, file=output14, status='replace')
	do i = 1, barrier_rows
		do j = 1, barrier_columns
			write(dunit, '((E19.12, 1x))', advance = 'no') krec_matrix(i, j)
		end do
		write(dunit, *)
	end do	
	close(dunit)

! Compute and Print Final J,K distributions here !

	if (myid == 0) call cpu_time(end_time_master_eq)
	if (myid == 0) write(*, *) "Master equation propagation:", end_time_master_eq-start_time_master_eq, "seconds"
	
	if (myid == 0) call cpu_time(end_time)
	if (myid == 0) write(*, *) "CPU Time:", end_time-start_time, "seconds"

	call MPI_Finalize(ierr_a)
	
! Main part of the code ends here
!-------------------------------------------------------------------------------------------------------------
	contains
!-------------------------------------------------------------------------------------------------------------		  
	function calculate_transition_matrix_unitless(Energies, Cov, Vdw, Inf, dE_down) result(matrix)
!-------------------------------------------------------------------------------------------------------------	
! Calculates unitless state-to-state transition matrix (matrix(j, i) = kappa j->i)
! dE = [down, up]
		real*8, dimension(:), intent(in) :: Energies(:)
		real*8, dimension(:), intent(in) :: Cov(:), Vdw(:), Inf(:)
		real*8 :: dE_down, dE_up, dJ
		real*8 :: ptran, matrixtmp
		real*8, allocatable, dimension(:, :) :: matrix
		integer :: i, j, a, chunk_size_last

! Allocate a 2D array for a part of a matrix on each processor
		allocate(matrix(num_states, chunk_size))
		matrix = 0.d0
	
		dJ = 10.7894
! Remanant of old way of parallelization	
!		do i = myid*chunk_size+1, myid*chunk_size+chunk_size
!		a = i - myid*chunk_size

! New way of parallelization
	do i = myid+1, num_states, pe_num
		a = (i - (myid+1))/pe_num + 1		
		if (i > num_states) exit				
		do j = 1, num_states
			  if (j == i) cycle
! Transition Probability			  
			  ptran = Cov(i)*Cov(j) + Vdw(i)*Vdw(j)! + Inf(i)*Inf(j)
! Quenching is above the diagonal			  
			  if (j > i) then
				matrixtmp = ptran * exp((Energies(j) - Energies(i)) / dE_down) * exp(-abs(real(J_value(j))-real(J_value(i)))/dJ)
				if (matrixtmp > kij_min) matrix(j, a) = matrixtmp
			  else
! Excitation is below the diagonal, reversibility is used
				matrixtmp = ptran * exp((Energies(i) - Energies(j)) / dE_down)&
				* exp(-abs(real(J_value(j))-real(J_value(i)))/dJ) * eq_concs(i)/eq_concs(j)
				if (matrixtmp > kij_min* eq_concs(i)/eq_concs(j)) matrix(j, a) = matrixtmp
! Old way of including reversibility through excitation dE_up
!			 	dE_up = get_dE_up(dE_down, temp_k, Energies(j), Energies(i), equilibrium_constants_m3(j), equilibrium_constants_m3(i))				
!				matrix(j, a) = ptran * exp((Energies(i) - Energies(j)) / dE_up)
			  end if
		   end do   
 		end do
		
	end function calculate_transition_matrix_unitless
	   
!-------------------------------------------------------------------------------------------------------------	
	subroutine write_global_matrix_to_file(matrix, filename, num_procs)
!-------------------------------------------------------------------------------------------------------------
	  real*8, dimension(:, :, :) :: matrix
	  character(len=*), intent(in) :: filename
	  integer, intent(in) :: num_procs
	  integer :: unit, i, j, k, chunk_size_last

	  ! Open the file for writing
	  open(newunit=unit, file=filename, status='replace')

	  ! Print header
	  do j = 1, size(matrix, 1)
		write(unit, '(1x,I8,a1)', advance='no') j, ','
	  end do
	  write(unit, *)

	  ! Loop over each processor
	  do k = 1, num_procs
		! Loop over each row of the matrix and write to the file
!		do i = (k-1)*chunk_size+1, (k-1)*chunk_size+chunk_size
!		a = i - (k-1)*chunk_size
		do i = k, num_states, pe_num
		a = (i - k)/pe_num + 1
		
		if (myid == pe_num-1) then
			chunk_size_last = num_states - (pe_num-1)*chunk_size
			if (a > chunk_size_last) exit
		end if
		
!		if (i .gt. num_states) exit
		
		  ! Loop over each column of the matrix and write to the file
		  do j = 1, size(matrix, 1)
			if (j == 1) then
			  write(unit, '(I8,a1,1x,E19.12,a1,1x)', advance='no') i, ',', matrix(j, a, k), ','
			else
			  write(unit, '(E19.12,a1,1x)', advance='no') matrix(j, a, k), ','
			end if
		  end do
		  write(unit, *)
		end do
	  end do

	  ! Close the file
	  close(unit)
	end subroutine write_global_matrix_to_file

!-------------------------------------------------------------------------------------------------------------		 	
	subroutine derivs(t, C, dCdt)
!-------------------------------------------------------------------------------------------------------------
		use mpi
		integer :: myid_counter, counter
		real*8, intent(in) :: t, C(num_states)
		real*8, intent(out) :: dCdt(num_states)
		real*8 :: dCdt_myid(chunk_size), dCdt_myid2(chunk_size, pe_num)
		real*8 Third_term(chunk_size), C_myid(chunk_size)
		integer :: chunk_size_last

!The number of rows of the matrix on the last processor
		chunk_size_last = num_states - (pe_num-1)*chunk_size

! Remanant from the old way of parallelization by continous chunks	
!		do i = myid*chunk_size+1, myid*chunk_size+chunk_size
!	        a = i - myid*chunk_size

! New way of Parallelization by processor number, variable "a" is the equation number
! Loop over final states, same as equation numbers:
		Third_term = 0.d0
		do a = 1, variable_chunk_size
		i = (myid+1) + (a-1)*pe_num
		if (continium(i) == 0) then
!		   if (i > num_states) exit	   
! Remanant
!		do i = 1, chunk_size
! New way of computing third term by each processor
			do j = 1, idepth(a)  
				if (continium(state_matrix(j,a)) == 0) Third_term(a) = Third_term(a) + trans_element(j, a) * C(state_matrix(j,a))
			end do
		end if	
		end do
		
! Load a part of array "C" of concentration into a smaller array for each processor
		do a = 1, variable_chunk_size
		i = (myid+1) + (a-1)*pe_num
!			if (i > num_states) exit	
			C_myid(a) = C(i)
		end do

! Each processor computes RHS for its part of the system of equations	
		dCdt_myid = 0 	
		do a = 1, variable_chunk_size
			i = (myid+1) + (a-1)*pe_num
			dCdt_myid(a) = First_term(a) - (Second_term(a) + Fourth_term(a))*C_myid(a) + Third_term(a) + Fifth_term(a) 
! For debugging
!			write(myid*1000+1,'(I8,1x,E20.12,1x,E20.12,1x,E20.12,1x,E20.12,1x,E20.12)') a, dCdt_myid(a), First_term(a), Second_term(a), C_myid(a), Third_term(a)			
!			write(myid*1000+1,'(I8,1x,I8,1x,E20.12)') i, a, dCdt_myid(a)
		end do

! For debugging
!		write(myid*1000+1,'(9(E20.12))') (dCdt_myid(j), j = 1, num_states)

! Collecting RHS from all processesors into a 2D array.	
		call MPI_GATHER(dCdt_myid, chunk_size, MPI_REAL8, dCdt_myid2, chunk_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)			
		call MPI_BARRIER( MPI_COMM_WORLD, ierr_a )

! Straightening of a 2D array into 1D array of RHSs								 
		if (myid == 0 ) then 
			do i = 1, chunk_size
				do j = 1, pe_num
					counter = (i-1)*pe_num + j
					if (counter .gt. num_states) goto 90
					dCdt(counter) = dCdt_myid2(i, j)
!					write(1000, '(I8, 1x, E20.12)') counter, dCdt(counter)
				end do
			 end do
90 		end if 

! Alternative way of straightening 
!		if (myid == 0 ) then 
!			do counter = 1, num_states
!				j = int(counter/chunk_size) + 1
!				i = counter - chunk_size*(j-1)
!				dCdt(counter) = dCdt_myid2(i, j)
!			 end do
!		end if 

! Root processor broadcasts RHS array to all processors				
		call MPI_BCAST(dCdt, num_states, MPI_REAL8, 0, MPI_COMM_WORLD, ierr_a)
		call MPI_BARRIER(MPI_COMM_WORLD, ierr_a)
							
	end subroutine derivs

!-------------------------------------------------------------------------------------------------------------	
	subroutine rk4(y,dydx,n,x,h,yout,derivs)
!-------------------------------------------------------------------------------------------------------------
! Taken from Numerical Recipes
	   integer n,NMAX  
	   real*8 h,x,dydx(n),y(n),yout(n)  
	   EXTERNAL derivs  
	   PARAMETER (NMAX=1000000) ! Increase if the number of states is greater than 1M
	   integer i  
	   real*8 h6, hh, xh, dym(NMAX), dyt(NMAX), yt(NMAX)

	   hh=h*0.5d0  
	   h6=h/6.d0  
	   xh=x+hh  

	   do 11 i=1,n  
		 yt(i)=y(i)+hh*dydx(i)  
	11    continue  

	   call derivs(xh,yt,dyt)  
	   
	   do 12 i=1,n  
		 yt(i)=y(i)+hh*dyt(i)  
	12    continue  

	   call derivs(xh,yt,dym)  
	   
	   do 13 i=1,n  
		 yt(i)=y(i)+h*dym(i)  
		 dym(i)=dyt(i)+dym(i)  
	13    continue  

	   call derivs(x+h,yt,dyt)  

	   do 14 i=1,n  
		 yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))  
	14    continue  
	   return  
	end subroutine rk4

!-------------------------------------------------------------------------------------------------------------	
	subroutine get_threshold_energy_K(o3_molecule, o2_molecule, K, threshold_E, threshold_j) 
!-------------------------------------------------------------------------------------------------------------	
		character(len=3), intent(in) :: o3_molecule
		character(len=2), intent(in) ::  o2_molecule
		integer, intent(in) :: K
		real*8 :: threshold_energy_j
		real*8, intent(out) :: threshold_E
		integer, intent(out) :: threshold_j
		!real*8 :: channel_shift_j
		!channel_shift_j = get_channel_shift(o3_molecule, o2_molecule)
		threshold_j = get_o2_threshold_j(o2_molecule, K)
		threshold_energy_j = rigid_rotor_energy(threshold_j, get_inertia_moment_2(o2_molecule)) !+ channel_shift_j
		threshold_E = threshold_energy_j / j_per_cm
	end subroutine get_threshold_energy_K

!-------------------------------------------------------------------------------------------------------------		
	function get_o2_threshold_j(o2_molecule, K) result(o2_threshold_j)
!-------------------------------------------------------------------------------------------------------------	
		character(len=2)  o2_molecule
		integer :: K
		integer :: o2_threshold_j
		
		if (mod(K, 2) == 0 .and. is_monoisotopic(o2_molecule)) then
		  o2_threshold_j = K + 1
		else 
		  o2_threshold_j = K
		end if
		
	end function get_o2_threshold_j

!-------------------------------------------------------------------------------------------------------------		
	function get_inertia_moment_2(o2_molecule) result(I_kg_m2)
!-------------------------------------------------------------------------------------------------------------	
		character(len=2), intent(in) :: o2_molecule
		real*8 :: mu_rot_kg, I_kg_m2
		mu_rot_kg = get_mu_rot(o2_molecule)
		I_kg_m2 = get_inertia_moment(mu_rot_kg)
	end function get_inertia_moment_2

!-------------------------------------------------------------------------------------------------------------		
	function calculate_formation_decay_equilibrium(Energy, temp_k, part_funcs_o2_per_m3, J_value, K_value) result(Kfds_m3)
!-------------------------------------------------------------------------------------------------------------	
		real*8 :: temp_k, part_funcs_o2_per_m3, threshold_E, threshold_energy_j
		real*8 :: Energy, Total_Energy_j, J, j_per_k, kt_energy_j, Kfds_m3
		integer :: J_value, K_value
		
		call get_threshold_energy_K(o3_molecule, o2_molecule, K_value, threshold_E, threshold_j)
		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k;
		Total_Energy_j = Energy * j_per_cm
		threshold_energy_j = threshold_E * j_per_cm
			
		Kfds_m3 = (2*J_value+1)*exp(-(Total_Energy_j - global_dis_energy) / kt_energy_j) / part_funcs_o2_per_m3
		
	end function calculate_formation_decay_equilibrium
!-------------------------------------------------------------------------------------------------------------		
	function get_dE_up(dE_down, temp_k, Energy_i, Energy_j, K_eq_i, K_eq_j) result(dE_up)
!-------------------------------------------------------------------------------------------------------------	
! Calculates dE_up corresponding to dE_down to satisfy the reversibility principle
		real*8, intent(in) :: dE_down, temp_k, Energy_i, Energy_j, K_eq_i, K_eq_j
		real*8 :: j_per_k, kt_energy_j, dE_up , dE_down_j

		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k
		
!		dE_down_j = dE_down * j_per_cm
!		dE_up = (dE_down_j/j_per_cm) / (dE_down_j / kt_energy_j - 1)
		dE_up = 1.d0 / ( log(K_eq_j/K_eq_i)/(Energy_i-Energy_j) - 1.d0/dE_down )
		
	end function get_dE_up

!-------------------------------------------------------------------------------------------------------------		
	function calc_part_func_O2_per_m3_elec(temp_k) result(pfunc)
!-------------------------------------------------------------------------------------------------------------	
! returns electronic partition function of O2 (+O?) system
! The expression is taken from:
! Hathorn, B. C.; Marcus, R. A. 
! An Intramolecular Theory of the Mass-Independent Isotope Effect for Ozone. 
! II. Numerical Implementation at Low Pressures Using a Loose Transition State. 
! J. Chem. Phys. 2000, 113 (21), 9497–9509. 
! https://doi.org/10.1063/1.480267
		real*8, intent(in) :: temp_k
		real*8 :: pfunc
		real*8 :: j_per_k, j_per_cm, cm_per_k, kt_energy_cm
		j_per_k = get_j_per_k()
		j_per_cm = get_j_per_cm()
		cm_per_k = j_per_k / j_per_cm
		kt_energy_cm = temp_k * cm_per_k
		pfunc = 15 + 9 * exp(-158.5 / kt_energy_cm) + 3 * exp(-226.5 / kt_energy_cm)
	end function calc_part_funC_O2_per_m3_elec

!-------------------------------------------------------------------------------------------------------------		
	function get_j_per_k() result(j_per_k)
!-------------------------------------------------------------------------------------------------------------	
		real*8 :: j_per_k
		j_per_k = 1.3806488e-23
	end function get_j_per_k
		
!-------------------------------------------------------------------------------------------------------------			
	function get_j_per_cm() result(j_per_cm)
!-------------------------------------------------------------------------------------------------------------	
		real*8 :: j_per_cm
		j_per_cm = 1.986445682933742e-23
	end function get_j_per_cm
	
!-------------------------------------------------------------------------------------------------------------		
	function calc_part_func_O2_per_m3_vib(zpe_j, kt_energy_j) result(pfunc)
!-------------------------------------------------------------------------------------------------------------			
		real*8 :: zpe_j, kt_energy_j, pfunc
		pfunc = 1.0 / (1.0 - exp(-2.0 * zpe_j / kt_energy_j))
	end function calc_part_func_O2_per_m3_vib

!-------------------------------------------------------------------------------------------------------------		
	function calc_part_func_O2_per_m3_rot(mu_rot_kg, J_start, J_step, kt_energy_j, K) result(pfunc)
!-------------------------------------------------------------------------------------------------------------	
! mu_rot_kg - reduced mass of the O2 system
		integer, intent(in) :: J_start, J_step
		real*8, intent(in) :: mu_rot_kg, kt_energy_j
		real*8 :: eps, I_kg_m2, threshold_energy_j, energy_j, new_pfunc, pfunc, threshold_E
		integer :: J, threshold_j, K
		eps = 1e-10
		I_kg_m2 = get_inertia_moment(mu_rot_kg)
!		threshold_energy_j = rigid_rotor_energy(J_start, I_kg_m2)
		call get_threshold_energy_K(o3_molecule, o2_molecule, K, threshold_E, threshold_j)
		threshold_energy_j = threshold_E * j_per_cm
		
		pfunc = 0.0
		J = J_start
		do while (.true.)
		energy_j = rigid_rotor_energy(J, I_kg_m2)
		new_pfunc = pfunc + (2*J + 1) * exp(-(energy_j - global_dis_energy) / kt_energy_j)
		if (new_pfunc - pfunc < eps) then
		  exit
		end if
		pfunc = new_pfunc
		J = J + J_step
		end do
	end function calc_part_func_O2_per_m3_rot

!-------------------------------------------------------------------------------------------------------------		
	function get_inertia_moment(mu_rot_kg) result(I_kg_m2)
!-------------------------------------------------------------------------------------------------------------		
! Returns moment of inertia of a given O2 molecule
		real*8, intent(in) :: mu_rot_kg
		real*8 :: r0_m, I_kg_m2
		r0_m = get_o2_distance_a0() * get_m_per_a0()
		I_kg_m2 = mu_rot_kg * r0_m**2
	end function get_inertia_moment

!-------------------------------------------------------------------------------------------------------------		
	function get_o2_distance_a0() result(dist_a0)
!-------------------------------------------------------------------------------------------------------------		
! Returns equilibrium distance between atoms in an O2 molecule
		real*8 :: dist_a0
		dist_a0 = 2.2819
	end function get_o2_distance_a0

!-------------------------------------------------------------------------------------------------------------		
	function get_m_per_a0() result(m_per_a0)
!-------------------------------------------------------------------------------------------------------------		
! Returns number of meters in one Bohr
		real*8 :: m_per_a0
		m_per_a0 = 5.2917721092e-11
	end function get_m_per_a0

!-------------------------------------------------------------------------------------------------------------	
	function rigid_rotor_energy(J, I_kg_m2) result(energy_J)
!-------------------------------------------------------------------------------------------------------------	
! returns energy (in J) of a rigid rotor corresponding to given J and I (moment of inertia)
		integer J
		real*8 :: I_kg_m2, energy_J
		real*8 :: hbar_js
		hbar_js = get_hbar_js()
		energy_J = J * (J + 1) / (2 * I_kg_m2) * hbar_js**2
	end function rigid_rotor_energy

!-------------------------------------------------------------------------------------------------------------			
	function get_hbar_js() result(hbar_js)
!-------------------------------------------------------------------------------------------------------------	
		real*8 :: hbar_js
		hbar_js = 1.054571726d-34
	end function get_hbar_js

!-------------------------------------------------------------------------------------------------------------		
	function calc_part_func_O2_per_m3_trans(mu_trans_kg, kt_energy_j) result(pfunc)
!-------------------------------------------------------------------------------------------------------------	
! mu_trans_kg - reduced mass of the O+O2 system
		real*8 :: mu_trans_kg, kt_energy_j, hbar_js, de_broglie_wl, pfunc
		hbar_js = get_hbar_js()
		de_broglie_wl = sqrt(2 * pi / (mu_trans_kg * kt_energy_j)) * hbar_js
		pfunc = de_broglie_wl**(-3)
	end function calc_part_func_O2_per_m3_trans

!-------------------------------------------------------------------------------------------------------------							
	function calc_part_func_O2_per_m3_total(temp_k, zpe_j, mu_rot_kg, J_rot_start, J_rot_step, mu_trans_kg, K) result(pfunc)
!-------------------------------------------------------------------------------------------------------------	
! Calculates the overall partition function of O+O2 system
		real*8 :: temp_k, zpe_j, mu_rot_kg, mu_trans_kg
		real*8 :: j_per_k, kt_energy_j, pfunc_elec, pfunc_vib, pfunc_rot, pfunc_trans, pfunc
		integer J_rot_start, J_rot_step, K
		j_per_k = get_j_per_k()
		kt_energy_j = temp_k * j_per_k
		pfunc_elec = calc_part_func_O2_per_m3_elec(temp_k)
		pfunc_vib = calc_part_func_O2_per_m3_vib(zpe_j, kt_energy_j)
		pfunc_rot = calc_part_func_O2_per_m3_rot(mu_rot_kg, J_rot_start, J_rot_step, kt_energy_j, K)
		pfunc_trans = calc_part_func_O2_per_m3_trans(mu_trans_kg, kt_energy_j)
		pfunc = pfunc_elec * pfunc_vib * pfunc_rot * pfunc_trans
	end function calc_part_func_O2_per_m3_total

!-------------------------------------------------------------------------------------------------------------		
	function calc_part_func_O2_per_m3_total_mol(temp_k, o2_molecule, o_atom, J_rot_start, K) result(pfunc_m_3)
!-------------------------------------------------------------------------------------------------------------		
		real*8 :: temp_k, pfunc_m_3
		character(len=*), intent(in) :: o2_molecule, o_atom
		real*8 :: zpe_j, mu_rot_kg, mu_trans_kg
		integer J_rot_start, J_rot_step, K
		
		zpe_j = get_channel_zpe(o2_molecule)
		mu_rot_kg = get_mu_rot(o2_molecule)
		mu_trans_kg = get_mu_trans(o2_molecule, o_atom)
		if (is_monoisotopic(o2_molecule)) then
			J_rot_step = 2
			else
			J_rot_step = 1
		end if
		
		pfunc_m_3 = calc_part_func_O2_per_m3_total(temp_k, zpe_j, mu_rot_kg, J_rot_start, J_rot_step, mu_trans_kg, K)
	end function calc_part_func_O2_per_m3_total_mol

!-------------------------------------------------------------------------------------------------------------		
	function is_monoisotopic(molecule) result(res)
!-------------------------------------------------------------------------------------------------------------		
	  character(len=*), intent(in) :: molecule
	  logical :: res
	  integer :: i, molecule_len

	  molecule_len = len(molecule)
	  res = .true.
	  do i = 2, molecule_len
		if (molecule(i:i) /= molecule(1:1)) then
		  res = .false.
		  exit
		end if
	  end do
	end function is_monoisotopic
	
!-------------------------------------------------------------------------------------------------------------		
	function get_channel_zpe(o2_molecule) result(zpe_j)
!-------------------------------------------------------------------------------------------------------------		
! Returns zero-point energy of a given channel
		character(len=*), intent(in) :: o2_molecule
		real*8 :: zpe_cm, zpe_j

		if (o2_molecule == "66") then
		zpe_cm = 7.916382691754641e+02
		elseif (o2_molecule == "67") then
		zpe_cm = 7.799050607081324e+02
		elseif (o2_molecule == "68") then
		zpe_cm = 7.693708543295361e+02
		elseif (o2_molecule == "77") then
		zpe_cm = 7.679904668774815e+02
		elseif (o2_molecule == "78") then
		zpe_cm = 7.572885872707559e+02
		elseif (o2_molecule == "88") then
		zpe_cm = 7.464315071358510e+02
		endif

		zpe_j = zpe_cm * get_j_per_cm()
	end function get_channel_zpe

!-------------------------------------------------------------------------------------------------------------				
	function get_mu_rot(o2_molecule) result(mu_rot_kg)
!-------------------------------------------------------------------------------------------------------------		
! Returns reduced mass of a given O2 molecule
		character(len=*), intent(in) :: o2_molecule
		real*8 :: atom_masses_kg(2), mu_rot_kg

		atom_masses_kg = get_atom_masses(o2_molecule)
		mu_rot_kg = product(atom_masses_kg) / sum(atom_masses_kg)
	end function get_mu_rot

!-------------------------------------------------------------------------------------------------------------				
	function get_atom_masses(molecule) result(atom_masses)
!-------------------------------------------------------------------------------------------------------------	
		character(len=*), intent(in) :: molecule
		real*8, dimension(:), allocatable :: atom_masses
		real*8, dimension(3) :: oxygen_kg
		integer :: i, len

		oxygen_kg = get_oxygen_mass_amu() * get_kg_per_amu()
		len = len_trim(molecule)
		allocate(atom_masses(len))

		do i = 1, len
		select case (molecule(i:i))
		case ('6')
		  atom_masses(i) = oxygen_kg(1)
		case ('7')
		  atom_masses(i) = oxygen_kg(2)
		case ('8')
		  atom_masses(i) = oxygen_kg(3)
		end select
		end do
	end function get_atom_masses

!-------------------------------------------------------------------------------------------------------------					
	function get_oxygen_mass_amu() result(oxygen_amu)
!-------------------------------------------------------------------------------------------------------------	
! Returns masses of oxygen isotopoes (in amu)
		real*8, dimension(3) :: oxygen_amu

		oxygen_amu = [15.99491461956, 16.9991317, 17.9991596129]
	end function get_oxygen_mass_amu

!-------------------------------------------------------------------------------------------------------------		
	function get_nitrogen_mass_amu() result(nitrogen_amu)
!-------------------------------------------------------------------------------------------------------------	
! Returns masses of nitrogen isotopoes (in amu)
		real*8,dimension(2) :: nitrogen_amu
		
		nitrogen_amu = [14.0030740048, 15.0001088982]
	end function get_nitrogen_mass_amu

!-------------------------------------------------------------------------------------------------------------							
	function get_kg_per_amu() result(kg_per_amu)
!-------------------------------------------------------------------------------------------------------------	
		real*8 :: kg_per_amu

		kg_per_amu = 1.660538921d-27
	end function get_kg_per_amu

!-------------------------------------------------------------------------------------------------------------			
	function get_k0_2(o3_molecule, temp_k, sigma0_m2) result(k0_m3_per_s)
!-------------------------------------------------------------------------------------------------------------	
		real*8, intent(in) :: temp_k, sigma0_m2
		character(len=3), intent(in) :: o3_molecule
		real*8 :: ozone_mass_kg, third_body_mass_kg, kt_energy_j, k0_m3_per_s, kg_per_amu
		real*8, dimension(2) :: nitrogen_amu
		
		j_per_k = get_j_per_k();
		kg_per_amu = get_kg_per_amu();
		nitrogen_amu = get_nitrogen_mass_amu();
		
		ozone_mass_kg = 7.96805347915539e-26
		third_body_mass_kg = nitrogen_amu(1) * 2 * kg_per_amu
		kt_energy_j = temp_k * j_per_k
		k0_m3_per_s = get_k0(ozone_mass_kg, third_body_mass_kg, kt_energy_j, sigma0_m2)
		
	end function get_k0_2
	
!-------------------------------------------------------------------------------------------------------------		
	function get_k0(ozone_mass_kg, third_body_mass_kg, kt_energy_j, sigma_stab_m2) result(k_stab_m3_per_s)
!-------------------------------------------------------------------------------------------------------------	
! Calculates kstab corresponding to sigma_stab (stabilization cross-section)
! ozone_mass is the sum of masses of individual atoms, NOT the reduced mass (mu)
		real(8), intent(in) :: ozone_mass_kg, third_body_mass_kg, kt_energy_j, sigma_stab_m2
		real(8) :: mu_stab_kg, velocity_m_per_s, k_stab_m3_per_s

		mu_stab_kg = ozone_mass_kg * third_body_mass_kg / (ozone_mass_kg + third_body_mass_kg)
		velocity_m_per_s = sqrt(8 * kt_energy_j / (pi * mu_stab_kg))
		k_stab_m3_per_s = sigma_stab_m2 * velocity_m_per_s
	end function get_k0
	
!-------------------------------------------------------------------------------------------------------------		
	function get_mu_trans(o2_molecule, o_atom) result(mu_trans_kg)
!-------------------------------------------------------------------------------------------------------------		
		character(len=*), intent(in) :: o2_molecule, o_atom
		real*8 :: mu_trans_kg
		real*8, dimension(:), allocatable :: o2_atom_masses_kg, o_mass_kg
		real*8 :: o2_mass_kg, o_mass_sum_kg

		o2_atom_masses_kg = get_atom_masses(o2_molecule)
		o2_mass_kg = sum(o2_atom_masses_kg)

		o_mass_kg = get_atom_masses(o_atom)
		o_mass_sum_kg = sum(o_mass_kg)

		mu_trans_kg = o2_mass_kg * o_mass_sum_kg / (o2_mass_kg + o_mass_sum_kg)
	end function get_mu_trans
		
end PROGRAM RK_solution
