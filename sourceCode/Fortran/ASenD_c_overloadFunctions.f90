 module ASenD_c_overloadFunctions
     use ASenD_constantVals
 
     contains
 	
 	subroutine c_greater(gt,inp1,inp2)
 	    implicit none
 		
 		integer, intent(out) :: gt
 		complex*16, intent(in) :: inp1, inp2
 		
 		complex*16 :: i2shft
 		
 		i2shft = inp2 + 1e-1*abs(inp1-inp2)*c_1
 		
 		if(abs(inp1-i2shft) .lt. abs(inp1-inp2)) then
 		    gt = 1
 		else
 		    gt = 0
 		endif
 		
 	end subroutine c_greater
 
     subroutine c_sqrt(sqt,inpt)
 	    implicit none
 		
 		complex*16, intent(out) :: sqt
 		complex*16, intent(in) :: inpt
 		
 		complex*16 :: del, fVal, dF
 		integer :: ct
 		
 		sqt = c_1
 		del = c_1
 		ct = 0
 		do while(abs(del) .gt. 1e-13 .and. ct .lt. 50)
 		    fVal = sqt*sqt - inpt
 			dF = c_2*sqt
 			del = -fVal/dF
 			sqt = sqt + del
 			ct = ct + 1
 		enddo
 		
 	end subroutine c_sqrt
 
     subroutine c_sin0_90(sT,theta)
 	    implicit none
 		
 		complex*16, intent(out) :: sT
 		complex*16, intent(in) :: theta
 		
 		complex*16 :: coeffs(7), thetSq
 		integer :: i1
 		
 		
 		if(abs(theta) .lt. 0.7853981633974483) then
 		    thetSq = theta*theta
 		    coeffs(1) = c_1
 			coeffs(2) = -c_1*1.666666666666666e-1
 			coeffs(3) = c_1*8.333333333333333e-3
 			coeffs(4) = -c_1*1.984126984126984e-4
 			coeffs(5) = c_1*2.755731922398589e-6
 			coeffs(6) = -c_1*2.505210838544171e-8
 			coeffs(7) = c_1*1.6059043836821614e-10
 			sT = coeffs(7)
 			do i1 = 6, 1, -1
 			    sT = sT*thetSq + coeffs(i1)
 			enddo
 			sT = sT*theta
 		else
 		    thetSq = (theta-0.5d0*c_pi)*(theta-0.5d0*c_pi)
 		    coeffs(1) = c_1
 			coeffs(2) = -c_p5
 			coeffs(3) = c_1*4.1666666666666666e-2
 			coeffs(4) = -c_1*1.388888888888889e-3
 			coeffs(5) = c_1*2.4801587301587301e-5
 			coeffs(6) = -c_1*2.755731922398589e-7
 			coeffs(7) = c_1*2.0876756987868098e-9
 			sT = coeffs(7)
 			do i1 = 6, 1, -1
 			    sT = sT*thetSq + coeffs(i1)
 			enddo
 		endif
 		
 	end subroutine c_sin0_90
 	
 	subroutine c_sin(sT,theta)
 	    implicit none
 		
 		complex*16, intent(out) :: sT
 		complex*16, intent(in) :: theta
 		
 		integer :: gt, i1
 		complex*16 :: mapTheta1, mapTheta2
 
         mapTheta1 = theta	
 	    call c_greater(gt,theta,2d0*c_pi)
 	    if(gt .eq. 1) then
 		    i1 = abs(c_p5*theta/c_pi)
 			mapTheta1 = theta-i1*2d0*c_pi
 		endif
 		
 		call c_greater(gt,c_0,theta)
 		if(gt .eq. 1) then
 		    i1 = abs(c_p5*theta/c_pi)
 			mapTheta1 = theta+(i1+1)*2d0*c_pi
 		endif
 		
 		call c_greater(gt,mapTheta1,1.5d0*c_pi)
 		if(gt .eq. 1) then
 			mapTheta2 = 2d0*c_pi - mapTheta1
 			call c_sin0_90(sT,mapTheta2)
 			sT = -sT
 			return
 		endif
 		
 		call c_greater(gt,mapTheta1,c_pi)
 		if(gt .eq. 1) then
 			mapTheta2 = mapTheta1 - c_pi
 			call c_sin0_90(sT,mapTheta2)
 			sT = -sT
 			return
 		endif
 		
 		call c_greater(gt,mapTheta1,0.5d0*c_pi)
 		if(gt .eq. 1) then
 			mapTheta2 = c_pi - mapTheta1
 			call c_sin0_90(sT,mapTheta2)
 			return
 		endif
 		
 		call c_sin0_90(sT,mapTheta1)
 		
 		
 	end subroutine c_sin
 	
 	subroutine c_cos(cT,theta)
 	    implicit none
 		
 		complex*16, intent(out) :: cT
 		complex*16, intent(in) :: theta
 		
 		complex*16 :: thetaAdj
 		
 		thetaAdj = theta + 0.5d0*c_pi
 		
 		call c_sin(cT,thetaAdj)
 		
 	end subroutine c_cos
 	
 	subroutine c_asin(aS,sn)
 	    implicit none
 		
 		complex*16, intent(in) :: sn
 		complex*16, intent(out) :: aS
 
         complex*16 :: dt, fVal, dF
 		integer :: i1
 		
 		aS = c_0
 		dt = c_1
 		i1 = 0
 		do while(abs(dt) .gt. 1e-13 .and. i1 .lt. 20)
             call c_sin(fVal,aS)
 			fVal = fVal - sn
 			call c_cos(dF,aS)
 			dt = -fVal/dF
 			aS = aS + dt
 			i1 = i1 + 1
 		enddo
 		
 	end subroutine c_asin
 	
 	subroutine c_acos(aC,cs)
 	    implicit none
 		
 		complex*16, intent(in) :: cs
 		complex*16, intent(out) :: aC
 
         complex*16 :: dt, fVal, dF
 		integer :: i1
 		
 		aC = 0.5d0*c_pi
 		dt = c_1
 		i1 = 0
 		do while(abs(dt) .gt. 1e-13 .and. i1 .lt. 20)
             call c_cos(fVal,aC)
 			fVal = fVal - cs
 			call c_sin(dF,aC)
 			dF = -dF
 			dt = -fVal/dF
 			aC = aC + dt
 			i1 = i1 + 1
 		enddo
 		
 	end subroutine c_acos
 
 end module ASenD_c_overloadFunctions
