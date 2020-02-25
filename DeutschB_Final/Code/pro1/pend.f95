module setup
	use numtype
	implicit none
	
	integer, parameter :: n_eq = 3
	real(dp):: a, b, c

end module setup

program pendulum
	use setup
	implicit none
	real(dp) :: t, dt, tmax
	real(dp), dimension(n_eq) :: y
	
	a = 0.2_dp			!FIXED
	b = 0.2_dp			!FIXED
	c = 5.7_dp       	!Part 1 alternating around 5.7
	
	tmax    = 100
	t       = 0._dp
	dt      = 0.01_dp
	
	!intial 
	y(1)   = -1._dp			!X		Part 2 hold (ABC) fixed 
	y(2)   = 0._dp 	 		!Y		alter (xyz) intial conditions
	y(3)   = 0._dp 	 		!Z
	
	do while (t<tmax) 

	
		write(11,*) t, y(2), y(1)
		write(22,*) t, y(3)
		write(33,*) y(2), y(1)
		write(44,*) y(2), y(1), y(3)
		
		call rk4step(t, y, dt)
	
	end do


end program pendulum