module mod
	use numtype
	use integr
	implicit none
	
	real(dp), parameter :: w    = 1._dp
	real(dp), parameter :: w2   = 2._dp
	real(dp), parameter :: hbar = 1._dp
	real(dp), parameter :: mass = 1._dp
	real(dp), parameter :: tt    = 1._dp     ! Varied at 1 and 5 
 	integer :: i, k, j, m, info
	integer, parameter :: ndim = 5, lwork= 16*ndim 
	real(dp) :: a, b, c, d, dia, work(lwork) 
	real(dp) :: res, resdos
	integer :: nint, ifail
	real(dp), dimension(maxint) :: xint, wint

	real(dp), dimension (0:ndim, 0:ndim) :: E, BigU, BigV, T, H, P
	real(dp), dimension (1:ndim+1, 1:ndim+1) :: PH, PH2, COM, TE
	complex(dp), dimension (1:ndim+1, 1:ndim+1) :: ZZ, EXPH
	real(dp), dimension (1:ndim+1) :: O

end module mod

Program BUharmo
	use mod
	real(dp), external :: psiguy1, psiguy2
	real(dp), external :: U, V

!MAKING EXPLICIT E and P MATRIX
		do  m = 0, ndim
			do  j = 0, ndim 
				dia = 0.5_dp
				if ( m==j ) then                                                                                 
                	E(m,j) = (dia + j)*w
            	else 
					E(m,j) = 0._dp
				end if 
			end do 
		end do
		do m = 0, ndim
			do j = 0, ndim 
				if (abs(m-j)==5) then 
				P(m,j) = (hbar*w2)/2
				else 
				P(m,j) = 0._dp
				end if
			end do
		end do 

print *, '--------------E MATRIX-----------------'  		
	 do i =0, ndim
	 	print '(''(''11(f8.3, 3x)'')'')', E(i, 0:ndim)
	 end do 

print *, '--------------P MATRIX-----------------'  		
	 do i =0, ndim
	 	print '(''(''11(f8.3, 3x)'')'')', P(i, 0:ndim)
	 end do 

!INTEGRATION
	nint = 100
	a = -10._dp
    b = 10._dp
	c = 0._dp
	d = 0._dp

    call d01bcf(0,a,b,c,d,nint,wint,xint,ifail)

	do i = 0, ndim
		do j= 0, ndim 
			res = 0._dp
			resdos = 0._dp
			do k = 0, nint
        		res = res + wint(k)*(psiguy1(xint(k),i)) * psiguy1(xint(k), j) * U(xint(k))  
				resdos = resdos + wint(k)*(psiguy2(xint(k),i)) * psiguy2(xint(k), j) * V(xint(k))  
			end do
			bigU(i, j) = res
			bigV(i, j) = resdos
		end do
	end do 	

!CREATION OF BIG ELEMENTS OF EQUATIONS
print *, '--------------BIG U MATRIX-----------------'	
 	do i =0, ndim
		print '(''(''11(f8.3, 3x)'')'')', bigU(i, 0:ndim)
	end do 

print *, '--------------Big V MATRIX-----------------'	
 	do i =0, ndim
		print '(''(''11(f8.3, 3x)'')'')', bigV(i, 0:ndim)
	end do 

!MATRIX EQUATIONS
T(0:ndim, 0:ndim)= E(0:ndim, 0:ndim)- bigU(0:ndim, 0:ndim)
H(0:ndim, 0:ndim)= T(0:ndim, 0:ndim)+ bigV(0:ndim, 0:ndim)
PH(1:ndim+1, 1:ndim+1) = H(0:ndim, 0:ndim) + P(0:ndim, 0:ndim)
PH2(1:ndim+1, 1:ndim+1) = PH(1:ndim+1, 1:ndim+1)

 print *, '--------------H MATRIX---------------------'	
 	do i =0, ndim
		print '(''(''11(f8.3, 3x)'')'')', H(i, 0:ndim)
	end do 

print *, '--------------PH MATRIX---------------------'	
 	do i =1, ndim+1
		print '(''(''11(f8.3, 3x)'')'')', PH(i, 1:ndim+1)
	end do 

print *, '--------------Eigenvalues--------------------'
call dsyev('v', 'u', ndim+1, PH2, ndim+1, O, work, lwork, info)
	do i = 1, ndim+1	
			print *, O(i)
	end do

print *, '---------------Eigenvector-------------------'	
 	do i =1, ndim+1
		print '(''(''11(f8.3, 3x)'')'')', PH2(i, 1:ndim+1)
	end do

print *, '--------------Completeness-------------------'
COM(1:ndim+1, 1:ndim+1) = matmul( PH2(1:ndim+1, 1:ndim+1), transpose(PH2(1:ndim+1,1:ndim+1)))
	do i = 1, ndim+1
		print '(''(''11(f8.3, 3x)'')'')', COM(i, 1:ndim+1)
	end do

print *, '---------------TEigenvector-------------------'
TE(1:ndim+1, 1:ndim+1)= transpose(PH(1:ndim+1, 1:ndim+1))
	do i =1, ndim+1
		print '(''(''11(f8.3, 3x)'')'')', TE(i, 1:ndim+1)
	end do

print *, '---------------ExpyDiaEigenguys-------------------'
do i = 1, ndim+1
ZZ(i,i) = exp(-iic * tt * O(i)) * COM(i, i)
print '(6(''('',f10.5,'','',f10.5,'')'',2x))', ZZ(i, 1:ndim+1)
end do 

print *, '-----------------ExpyH-------------------------------'
EXPH(1:ndim+1, 1:ndim+1) = matmul(PH2(1:ndim+1, 1:ndim+1), matmul(ZZ(1:ndim+1, 1:ndim+1), transpose(PH2(1:ndim+1, 1:ndim+1))))
print '(6(''('',f10.5,'','',f10.5,'')'',2x))', EXPH(1:ndim+1, 1:ndim+1)

End program BUharmo

!EXTERNAL FUNCTIONS---------------------------------------
Function psiguy1(x,n)
	use mod 
 	integer :: n,l
	real(dp) :: x, psiguy1
	real(dp), dimension(0:n) :: psi
	
		psi(0) = (mass*w/(hbar*pi))**(0.25_dp) * exp(-mass*w*x**2/(2*hbar))
		
	do l = 0, n-1	
 		psi(l+1) = (2*sqrt(mass*w/hbar)*x)/(2._dp*(l+1))**(0.5_dp)*psi(l)& 
 				 		   -(float(l)/(l+1))**(0.5_dp)*psi(l-1) 	 
	end do
			
		psiguy1 = psi(n)

End function 

Function psiguy2(x,n)
	use mod 
 	integer :: n,l
	real(dp) :: x, psiguy2
	real(dp), dimension(0:n) :: psi
	
		psi(0) = (mass*w2/(hbar*pi))**(0.25_dp) * exp(-mass*w2*x**2/(2*hbar))
		
	do l = 0, n-1	
 		psi(l+1) = (2*sqrt(mass*w2/hbar)*x)/(2._dp*(l+1))**(0.5_dp)*psi(l)& 
 				 		   -(float(l)/(l+1))**(0.5_dp)*psi(l-1) 	 
	end do
			
		psiguy2 = psi(n)
End function 

Function U(x) 
	use mod 
	real(dp) :: U, x

		U= 0.5_dp * mass * (w**2) * x**2 
End function 

Function V(x) 
	use mod
	real(dp) :: V, x

		V= 0.5_dp * mass * (w2**2) * x**2 
End function 