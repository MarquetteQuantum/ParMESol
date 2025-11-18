module barrier_energies_module
	implicit none
	
	integer, parameter :: rows=21
	integer, parameter :: columns = 65
	integer :: ios_barrier
	real*8, dimension(rows, columns) :: barrier_energies
	
contains

function read_matrix(filepath) result(matrix)

	character(len=*), intent(in) :: filepath
	real*8, dimension(rows, columns) :: matrix
	
	integer :: i,j
	real*8 :: temp
	
	open(10, file=filepath, status='old', action='read', iostat=ios_barrier)
	
	do i = 1, rows
		do j = 1, columns
			read(10, '(E21.14)', advance='no', iostat=ios_barrier) temp
			matrix(i, j) = temp
		end do
		read(10, *, iostat=ios_barrier)
	end do
	
	close(10)
	
	end function read_matrix
	
end module barrier_energies_module
	