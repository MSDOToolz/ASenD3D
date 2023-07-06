 module ASenD_c_designPropertyFunctions
     use ASenD_globalData
 	use ASenD_constantVals
 	use ASenD_c_overloadFunctions
 	
     contains
 	
     subroutine c_getMaterialStiffness(CMat,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: CMat(6,6)
 		
 		complex*16 :: eProps(9), stiffMat(21)
 		complex*16 :: SMat(6,6), fVec(6), zVec(6)
 		integer :: secNum, matNum
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		matNum = sectionMatId(secNum)
 		eProps(:) = c_1*materialElastic(:,matNum)
 		stiffMat(:) = c_1*materialStiffMat(:,matNum)
 		
 		if(abs(stiffMat(1)) .gt. 0d0) then
 		    do i1 = elToDRange(element-1) + 1, elToDRange(element)
 			    i2 = elToD(i1)
 				if(dCategory(i2) .eq. 'stiffnessMat') then
 				    i3 = dComponent(i2)
 					stiffMat(i3) = stiffMat(i3) + elToCoef(i1)*c_dVec(i2)
 				endif
 			enddo
 			CMat(1,:) = stiffMat(1:6)
 			CMat(2,2:6) = stiffMat(7:11)
 			CMat(3,3:6) = stiffMat(12:15)
 			CMat(4,4:6) = stiffMat(16:18)
 			CMat(5,5:6) = stiffMat(19:20)
 			CMat(6,6) = stiffMat(21)
 			do i1 = 2, 6
 			    do i2 = 1, i1-1
 				    CMat(i1,i2) = CMat(i2,i1)
 				enddo
 			enddo
         else	
 			do i1 = elToDRange(element-1) + 1, elToDRange(element)
 				i2 = elToD(i1)
 				if(dCategory(i2) .eq. 'modulus') then
 					i3 = dComponent(i2)
 					eProps(i3) = eProps(i3) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'poissonRatio') then
 					i3 = dComponent(i2) + 3
 					eProps(i3) = eProps(i3) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'shearModulus') then
 					i3 = dComponent(i2) + 6
 					eProps(i3) = eProps(i3) + elToCoef(i1)*c_dVec(i2)
 				endif
 			enddo
 			
 			SMat(:,:) = c_0
 			SMat(1,1) = c_1/eProps(1)
 			SMat(1,2) = -eProps(4)/eProps(1)
 			SMat(1,3) = -eProps(5)/eProps(1)
 			SMat(2,2) = c_1/eProps(2)
 			SMat(2,3) = -eProps(6)/eProps(2)
 			SMat(3,3) = c_1/eProps(3)
 			SMat(4,4) = c_1/eProps(7)
 			SMat(5,5) = c_1/eProps(8)
 			SMat(6,6) = c_1/eProps(9)
 			
 			do i1 = 2, 6
 				do i2 = 1, i1-1
 					SMat(i1,i2) = SMat(i2,i1)
 				enddo
 			enddo
 			
 			CMat(:,:) = c_0
 			
 			call c_GetInvLU(CMat,SMat,fVec,zVec,6,1,6)
 		endif
 	
 	end subroutine c_getMaterialStiffness
 	
 	subroutine c_getMaterialDensity(den,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: den
 		
 		integer :: secNum, matId
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		matId = sectionMatId(secNum)
 		
 		den = c_1*materialDensity(matId)
 		
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 		    i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'density') then
 			    den = den + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo
 		
 	end subroutine c_getMaterialDensity
 	
 	subroutine c_getMaterialThermExp(tExp,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: tExp(6)
 		
 		integer :: secNum, matId
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		matId = sectionMatId(secNum)
 		
 		tExp(:) = c_1*materialThermExp(:,matId)
 		
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 		    i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'thermalExp') then
 			    i3 = dComponent(i2)
 			    tExp(i3) = tExp(i3) + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo
 		
 	end subroutine c_getMaterialThermExp
 	
 	subroutine c_getMaterialThermCond(tCond,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: tCond(3,3)
 		
 		complex*16 :: tCVec(6)
 		integer :: secNum, matId
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		matId = sectionMatId(secNum)
 		
 		tCVec(:) = c_1*materialThermCond(:,matId)
 		
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 		    i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'thermalCond') then
 			    i3 = dComponent(i2)
 			    tCVec(i3) = tCVec(i3) + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo
 		
 		tCond(1,1) = tCVec(1)
 		tCond(2,2) = tCVec(2)
 		tCond(3,3) = tCVec(3)
 		tCond(1,2) = tCVec(4)
 		tCond(2,1) = tCVec(4)
 		tCond(1,3) = tCVec(5)
 		tCond(3,1) = tCVec(5)
 		tCond(2,3) = tCVec(6)
 		tCond(3,2) = tCVec(6)
 		
 	end subroutine c_getMaterialThermCond
 	
 	subroutine c_getMaterialSpecHeat(specHeat,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: specHeat
 		
 		integer :: secNum, matId
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		matId = sectionMatId(secNum)
 		
 		specHeat = c_1*materialSpecHeat(matId)
 		
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 		    i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'specHeat') then
 			    specHeat = specHeat + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo
 		
 	end subroutine c_getMaterialSpecHeat
 	
 	subroutine  c_getOrientation(alpha,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: alpha(3,3)
 		
 		complex*16 :: alOrig(3,3), locRot(3), globRot(3)
 		integer :: secNum
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		
 		locRot(:) = c_0
 		
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
             i2 = elToD(i1)
             if(dCategory(i2) .eq. 'orientation') then
 			    i3 = dComponent(i2)
 				locRot(i3) = locRot(i3) + elToCoef(i1)*c_dVec(i2)
             endif			
 		enddo
 		
 		locRot(:) = c_pi180*locRot(:)
 		
 		i1 = secNum*3
 		alOrig(:,:) = c_1*sectionOrient(:,i1-2:i1)
 		
 		globRot(:) = c_0
 		do i1 = 1,3
 		    do i2 = 1,3
 			    globRot(i1) = globRot(i1) + alOrig(i2,i1)*locRot(i2)
 			enddo
 		enddo
 		
 		call c_rotateAlpha(alpha,alOrig,globRot,0,0)
 		
 	end subroutine c_getOrientation
 	
 	subroutine c_getElThickOffset(thick,zOff,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: thick, zOff
 		
 		complex*16 :: layerThick
 		integer :: secNum, layerNum
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		
 		thick = c_0
 		layerNum = 0
 		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum)
 		    layerNum = layerNum + 1
 			layerThick = c_1*layupThickness(i1)
 			do i2 = elToDRange(element-1)+1, elToDRange(element)
 			    i3 = elToD(i2)
 			    if(dCategory(i3) .eq. 'thickness' .and. dLayer(i3) .eq. layerNum) then
 				    layerThick = layerThick + elToCoef(i2)*c_dVec(i3)
 				endif
 			enddo
 			thick = thick + layerThick
 		enddo
 		
 		zOff = c_1*sectionZOffset(secNum)
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 		    i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'zOffset') then
 			    zOff = zOff + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo
 		
 	end subroutine c_getElThickOffset
 	
 	subroutine c_getShellStiffness(ABDMat,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: ABDMat(9,9)
 		
 		complex*16 :: secThick, layerThick, layerProps(6)
 		complex*16 :: zCrd, zNext, zOff
 		complex*16 :: SMat(3,3), QMat(3,3), fVec(3), zVec(3)
 		complex*16 :: thetaRad, a11, a12, a21, a22
 		complex*16 :: Ts(3,3), Te(3,3), TeInv(3,3)
 		integer :: secNum, layerNum
 		integer :: i1, i2, i3, i4
 		
 		secNum = elementSection(element)
 		
 		call c_getElThickOffset(secThick,zOff,element)
 		
 		ABDMat(:,:) = c_0
 		zCrd = -c_p5*secThick*(zOff+c_1)  !! initialize to zMin
 		layerNum = 0
 		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum)
 		    layerNum = layerNum + 1
 			i2 = layupMatId(i1)
 			layerProps(1) = c_1*materialElastic(1,i2)
 			layerProps(2) = c_1*materialElastic(2,i2)
 			layerProps(3) = c_1*materialElastic(4,i2)
 			layerProps(4) = c_1*materialElastic(7,i2)
 			layerProps(5) = c_1*layupThickness(i1)
 			layerProps(6) = c_1*layupAngle(i1)
 			do i2 = elToDRange(element-1)+1, elToDRange(element)
 			    i3 = elToD(i2)
 			    if(dLayer(i3) .eq. layerNum) then
 				    if(dCategory(i3) .eq. 'modulus') then
 					    i4 = dComponent(i3)
 						layerProps(i4) = layerProps(i4) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'poissonRatio') then
 					    layerProps(3) = layerProps(3) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'shearModulus') then
 					    layerProps(4) = layerProps(4) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'thickness') then
 					    layerProps(5) = layerProps(5) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'angle') then
 					    layerProps(6) = layerProps(6) + elToCoef(i2)*c_dVec(i3)
 				    endif
 				endif
 			enddo
 			zNext = zCrd + layerProps(5)
 			SMat(:,:) = c_0
 			SMat(1,1) = c_1/layerProps(1)
 			SMat(1,2) = -layerProps(3)/layerProps(1)
 			SMat(2,2) = c_1/layerProps(2)
 			SMat(3,3) = c_1/layerProps(4)
 			SMat(2,1) = SMat(1,2)
 			! write(lfUnit,*) 'SMat: '
 			! write(lfUnit,*) SMat
 			call c_GetInvLU(QMat,SMat,fVec,zVec,3,1,3)
 			! write(lfUnit,*) 'QMat: '
 			! write(lfUnit,*) QMat
 			thetaRad = c_pi180*layerProps(6)
 			call c_cos(a11,thetaRad)
             call c_sin(a21,thetaRad)
 			! write(lfUnit,*) 'a11: ', a11, 'a21: ', a21
             a12 = -a21
             a22 = a11
             Ts(1,1) = a11*a11
             Ts(1,2) = a12*a12
             Ts(1,3) = c_2*a11*a12
             Ts(2,1) = a21*a21
             Ts(2,2) = a22*a22
             Ts(2,3) = c_2*a22*a21
             Ts(3,1) = a11*a21
             Ts(3,2) = a12*a22
             Ts(3,3) = a11*a22 + a12*a21
             Te(:,:) = Ts(:,:)
             Te(:,3) = c_p5*Te(:,3)
             Te(3,:) = c_2*Te(3,:)
 			!write(lfUnit,*) 'Te: '
 			!write(lfUnit,*) Te
             call c_GetInvLU(TeInv,Te,fVec,zVec,3,1,3)
 			!write(lfUnit,*) 'TeInv: '
 			!write(lfUnit,*) TeInv
             Te(:,:) = c_0
             do i2 = 1, 3
 			    do i3 = 1, 3
 				    do i4 = 1, 3
 					    Te(i2,i3) = Te(i2,i3) + QMat(i2,i4)*TeInv(i4,i3)					
 				    enddo
 				enddo
             enddo
             QMat(:,:) = c_0
             do i2 = 1, 3
 			    do i3 = 1, 3
 				    do i4 = 1, 3
 					    QMat(i2,i3) = QMat(i2,i3) + Ts(i2,i4)*Te(i4,i3)					
 				    enddo
 				enddo
             enddo		
 			ABDMat(1:3,1:3) = ABDMat(1:3,1:3) + (zNext-zCrd)*QMat
 			ABDMat(1:3,4:6) = ABDMat(1:3,4:6) + c_p5*(zNext*zNext - zCrd*zCrd)*QMat
 			ABDMat(4:6,4:6) = ABDMat(4:6,4:6) + c_1o3*(zNext*zNext*zNext - zCrd*zCrd*zCrd)*QMat
 			zCrd = zNext
 		enddo
 		
 		do i1 = 2,6
 		    do i2 = 1, i1-1
 			    ABDMat(i1,i2) = ABDMat(i2,i1)
 			enddo
 		enddo
 		
 		ABDMat(7,7) = ABDMat(3,3)
 		ABDMat(8,8) = ABDMat(3,3)
 		ABDMat(9,9) = ABDMat(3,3)
 		
 	end subroutine c_getShellStiffness
 	
 	subroutine c_getShellExpLoad(expLoad,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: expLoad(6)
 		
 		complex*16 :: secThick, layerThick, layerExp(6), layerProps(6)
 		complex*16 :: zCrd, zNext, zOff
 		complex*16 :: SMat(3,3), QMat(3,3), fVec(3), zVec(3)
 		complex*16 :: thetaRad, a11, a12, a21, a22
 		complex*16 :: Ts(3,3), Te(3,3), TeInv(3,3)
 		integer :: secNum, layerNum
 		integer :: i1, i2, i3, i4
 		
 		secNum = elementSection(element)
 		
 		call c_getElThickOffset(secThick,zOff,element)
 		
 		expLoad(:) = c_0
 		zCrd = -c_p5*secThick*(zOff+c_1)  !! initialize to zMin
 		layerNum = 0
 		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum)
 		    layerNum = layerNum + 1
 			i2 = layupMatId(i1)
 			layerExp(:) = c_1*materialThermExp(:,i2)
 			layerProps(1) = c_1*materialElastic(1,i2)
 			layerProps(2) = c_1*materialElastic(2,i2)
 			layerProps(3) = c_1*materialElastic(4,i2)
 			layerProps(4) = c_1*materialElastic(7,i2)
 			layerProps(5) = c_1*layupThickness(i1)
 			layerProps(6) = c_1*layupAngle(i1)
 			do i2 = elToDRange(element-1)+1, elToDRange(element)
 			    i3 = elToD(i2)
 			    if(dLayer(i3) .eq. layerNum) then
 				    if(dCategory(i3) .eq. 'thermalExp') then
 					    i4 = dComponent(i3)
 						layerExp(i4) = layerExp(i4) + elToCoef(i2)*c_dVec(i3)
 				    elseif(dCategory(i3) .eq. 'modulus') then
 					    i4 = dComponent(i3)
 						layerProps(i4) = layerProps(i4) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'poissonRatio') then
 					    layerProps(3) = layerProps(3) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'shearModulus') then
 					    layerProps(4) = layerProps(4) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'thickness') then
 					    layerProps(5) = layerProps(5) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'angle') then
 					    layerProps(6) = layerProps(6) + elToCoef(i2)*c_dVec(i3)
 				    endif
 				endif
 			enddo
 			layerExp(3) = layerExp(4)
 			zNext = zCrd + layerProps(5)
 			SMat(:,:) = c_0
 			SMat(1,1) = c_1/layerProps(1)
 			SMat(1,2) = -layerProps(3)/layerProps(1)
 			SMat(2,2) = c_1/layerProps(2)
 			SMat(3,3) = c_1/layerProps(4)
 			SMat(2,1) = SMat(1,2)
 			call c_GetInvLU(QMat,SMat,fVec,zVec,3,1,3)
 			thetaRad = c_pi180*layerProps(6)
 			call c_cos(a11,thetaRad)
             call c_sin(a21,thetaRad)
             a12 = -a21
             a22 = a11
             Ts(1,1) = a11*a11
             Ts(1,2) = a12*a12
             Ts(1,3) = c_2*a11*a12
             Ts(2,1) = a21*a21
             Ts(2,2) = a22*a22
             Ts(2,3) = c_2*a22*a21
             Ts(3,1) = a11*a21
             Ts(3,2) = a12*a22
             Ts(3,3) = a11*a22 + a12*a21
             Te(:,:) = Ts(:,:)
             Te(:,3) = c_p5*Te(:,3)
             Te(3,:) = c_2*Te(3,:)
             call c_GetInvLU(TeInv,Te,fVec,zVec,3,1,3)
             Te(:,:) = c_0
             do i2 = 1, 3
 			    do i3 = 1, 3
 				    do i4 = 1, 3
 					    Te(i2,i3) = Te(i2,i3) + QMat(i2,i4)*TeInv(i4,i3)					
 				    enddo
 				enddo
             enddo
             QMat(:,:) = c_0
             do i2 = 1, 3
 			    do i3 = 1, 3
 				    do i4 = 1, 3
 					    QMat(i2,i3) = QMat(i2,i3) + Ts(i2,i4)*Te(i4,i3)					
 				    enddo
 				enddo
             enddo
             do i2 = 1, 3
                 do i3 = 1, 3
 				    expLoad(i2) = expLoad(i2) + (zNext-zCrd)*Qmat(i2,i3)*layerExp(i3)
 					expLoad(i2+3) = expLoad(i2+3) + c_p5*(zNext*zNext - zCrd*zCrd)*Qmat(i2,i3)*layerExp(i3)
 				enddo
             enddo
 			zCrd = zNext
 		enddo
 		
 	end subroutine c_getShellExpLoad
 	
 	subroutine c_getShellMass(mMat,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: mMat(6,6)
 		
 		complex*16 :: secThick, layerThick, layerDen
 		complex*16 :: zCrd, zNext, zOff
 		complex*16 :: Ident(3,3), unitB(3,3)
 		integer :: secNum, layerNum
 		integer :: i1, i2, i3, i4
 		
 		ident(:,:) = c_0
 		ident(1,1) = c_1
 		ident(2,2) = c_1
 		ident(3,3) = c_1
 		
 		unitB(:,:) = c_0
 		unitB(1,2) = c_1
 		unitB(2,1) = -c_1
 		
 		secNum = elementSection(element)
 		
 		call c_getElThickOffset(secThick,zOff,element)
 		
 		mMat(:,:) = c_0
 		zCrd = -c_p5*secThick*(zOff+c_1)  !! initialize to zMin
 		layerNum = 0
 		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum)
 		    layerNum = layerNum + 1
 			i2 = layupMatId(i1)
 			layerDen = c_1*materialDensity(i2)
 			layerThick = c_1*layupThickness(i1)
 			do i2 = elToDRange(element-1)+1, elToDRange(element)
 			    i3 = elToD(i2)
 			    if(dLayer(i3) .eq. layerNum) then
 				    if(dCategory(i3) .eq. 'density') then
 				        layerDen = layerDen + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'thickness') then
 					    layerThick = layerThick + elToCoef(i2)*c_dVec(i3)
 					endif
 				endif
 			enddo
 			zNext = zCrd + layerThick
 			mMat(1:3,1:3) = mMat(1:3,1:3) + layerDen*(zNext-zCrd)*ident
 			mMat(1:3,4:6) = mMat(1:3,4:6) + c_p5*layerDen*(zNext*zNext - zCrd*zCrd)*unitB
 			mMat(4:6,4:6) = mMat(4:6,4:6) + c_1o3*layerDen*(zNext*zNext*zNext - zCrd*zCrd*zCrd)*ident
 			zCrd = zNext
 		enddo
 		
 		do i1 = 2,6
 		    do i2 = 1, i1-1
 			    mMat(i1,i2) = mMat(i2,i1)
 			enddo
 		enddo
 		
 	end subroutine c_getShellMass
 	
 	subroutine c_getShellThermCond(tCond,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: tCond(3,3)
 		
 		complex*16 :: secThick, layerThick, layerCond(6), layerAngle, thetaRad
 		complex*16 :: tCLoc(3,3), tCglob(3,3), prod(3,3)
 		complex*16 :: zCrd, zNext, zOff
 		complex*16 :: alpha(3,3)
 		integer :: secNum, layerNum
 		integer :: i1, i2, i3, i4
 		
 		secNum = elementSection(element)
 		
 		call c_getElThickOffset(secThick,zOff,element)
 		
 		tCond(:,:) = c_0
 		zCrd = -c_p5*secThick*(zOff+c_1)  !! initialize to zMin
 		layerNum = 0
 		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum)
 		    layerNum = layerNum + 1
 			i2 = layupMatId(i1)
 			layerCond(:) = c_1*materialThermCond(:,i2)
 			layerAngle = c_1*layupAngle(i1)
 			layerThick = c_1*layupThickness(i1)
 			do i2 = elToDRange(element-1)+1, elToDRange(element)
 			    i3 = elToD(i2)
 			    if(dLayer(i3) .eq. layerNum) then
 				    if(dCategory(i3) .eq. 'thermalCond') then
 					    i4 = dComponent(i3)
 						layerCond(i4) = layerCond(i4) + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'angle') then
 					    layerAngle = layerAngle + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'thickness') then
 					    layerThick = layerThick + elToCoef(i2)*c_dVec(i3)
 				    endif
 				endif
 			enddo
 			tCLoc(1,1) = layerCond(1)
 		    tCLoc(2,2) = layerCond(2)
 		    tCLoc(3,3) = layerCond(3)
 		    tCLoc(1,2) = layerCond(4)
 		    tCLoc(2,1) = layerCond(4)
 		    tCLoc(1,3) = layerCond(5)
 		    tCLoc(3,1) = layerCond(5)
 		    tCLoc(2,3) = layerCond(6)
 		    tCLoc(3,2) = layerCond(6)
 			alpha(:,:) = c_0
 			thetaRad = c_pi180*layerAngle
 			call c_cos(alpha(1,1),thetaRad)
             call c_sin(alpha(2,1),thetaRad)
             alpha(1,2) = -alpha(2,1)
             alpha(2,2) = alpha(1,1)
 			alpha(3,3) = c_1
             prod(:,:) = c_0
 			do i2 = 1, 3
 			    do i3 = 1, 3
 				    do i4 = 1, 3
 					    prod(i2,i3) = prod(i2,i3) + tCLoc(i2,i4)*alpha(i3,i4)
 					enddo
 				enddo
 			enddo
 			tCglob(:,:) = c_0
 			do i2 = 1, 3
 			    do i3 = 1, 3
 				    do i4 = 1, 3
 					    tCglob(i2,i3) = tCglob(i2,i3) + alpha(i2,i4)*prod(i4,i3)
 					enddo
 				enddo
 			enddo		
 			tCond(:,:) = tCond(:,:) + layerThick*tCglob(:,:)
 		enddo		
 		
 	end subroutine c_getShellThermCond
 	
     subroutine c_getShellSpecHeat(sH,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: sH
 		
 		complex*16 :: secThick, layerThick, layerSH, layerDen
 		complex*16 :: zCrd, zNext, zOff
 		integer :: secNum, layerNum
 		integer :: i1, i2, i3, i4
 		
 		secNum = elementSection(element)
 		
 		call c_getElThickOffset(secThick,zOff,element)
 		
 		sH = c_0
 		zCrd = -c_p5*secThick*(zOff+c_1)  !! initialize to zMin
 		layerNum = 0
 		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum)
 		    layerNum = layerNum + 1
 			i2 = layupMatId(i1)
 			layerSH = c_1*materialSpecHeat(i2)
 			layerDen = c_1*materialDensity(i2)
 			layerThick = c_1*layupThickness(i1)
 			do i2 = elToDRange(element-1)+1, elToDRange(element)
 			    i3 = elToD(i2)
 			    if(dLayer(i3) .eq. layerNum) then
 				    if(dCategory(i3) .eq. 'specHeat') then
 						layerSH = layerSH + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'density') then
 					    layerDen = layerDen + elToCoef(i2)*c_dVec(i3)
 					elseif(dCategory(i3) .eq. 'thickness') then
 					    layerThick = layerThick + elToCoef(i2)*c_dVec(i3)
 				    endif
 				endif
 			enddo
             sH = sH + layerThick*layerDen*layerSH
 		enddo
 		
 	end subroutine c_getShellSpecHeat
 	
 	subroutine c_getBeamSecProps(props,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: props(7)
 
 		integer :: secNum
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 
         props(:) = c_1*beamProperties(:,secNum)
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 			i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'area') then
 				props(1) = props(1) + elToCoef(i1)*c_dVec(i2)
 			elseif(dCategory(i2) .eq. 'areaMoment') then
 				i3 = dComponent(i2)
 				props(1+i3) = props(1+i3) + elToCoef(i1)*c_dVec(i2)
 			elseif(dCategory(i2) .eq. 'polarMoment') then
 				props(7) = props(7) + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo	
 		
 	end subroutine c_getBeamSecProps
 	
 	subroutine c_getBeamStiffness(Cmat,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: Cmat(6,6)
 		
 		complex*16 :: props(7), modulus, shearMod, stiffMat(21)
 		integer :: secNum
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		Cmat(:,:) = c_0
 		if(abs(beamStiffness(1,secNum)) .gt. 0d0) then
 		    stiffMat = c_1*beamStiffness(:,secNum)
 			do i1 = elToDRange(element-1)+1, elToDRange(element)
 			    i2 = elToD(i1)
 				if(dCategory(i2) .eq. 'stiffnessMat') then
 					i3 = dComponent(i2)
 					stiffMat(i3) = props(i3) + elToCoef(i1)*c_dVec(i2)
 				endif
 			enddo
 		    Cmat(1,1:6) = stiffMat(1:6)
 			Cmat(2,2:6) = stiffMat(7:11)
 			Cmat(3,3:6) = stiffMat(12:15)
 			Cmat(4,4:6) = stiffMat(16:18)
 			Cmat(5,5:6) = stiffMat(19:20)
 			Cmat(6,6) = stiffMat(21)
 			do i1 = 2, 6
 			    do i2 = 1, i1-1
 				    Cmat(i1,i2) = Cmat(i2,i1)
 				enddo
 			enddo
 		else
 		    props(:) = c_1*beamProperties(:,secNum)
 			i3 = sectionMatId(secNum)
 			modulus = c_1*materialElastic(1,i3)
 			shearMod = c_1*materialElastic(7,i3)
 		    do i1 = elToDRange(element-1)+1, elToDRange(element)
 			    i2 = elToD(i1)
 				if(dCategory(i2) .eq. 'area') then
 					props(1) = props(1) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'areaMoment') then
 					i3 = dComponent(i2)
 					props(1+i3) = props(1+i3) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'polarMoment') then
 					props(7) = props(7) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'modulus') then
 					modulus = modulus + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'shearModulus') then
 					shearMod = shearMod + elToCoef(i1)*c_dVec(i2)
 				endif
 			enddo
 			Cmat(1,1) = props(1)*modulus
 			Cmat(1,5) = props(2)*modulus
 			Cmat(1,6) = -props(3)*modulus
 			Cmat(2,2) = props(1)*shearMod
 			Cmat(2,4) = -props(2)*shearMod
 			Cmat(3,3) = props(1)*shearMod
 			Cmat(3,4) = props(3)*shearMod
 			Cmat(4,4) = props(7)*shearMod
 			Cmat(5,5) = props(4)*modulus
 			Cmat(5,6) = -props(6)*modulus
 			Cmat(6,6) = props(5)*modulus
 			do i1 = 2, 6
 			    do i2 = 1, i1-1
 				    Cmat(i1,i2) = Cmat(i2,i1)
 				enddo
 			enddo
 		endif
 		
 	end subroutine c_getBeamStiffness
 	
 	subroutine c_getBeamExpLoad(expLoad,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: expLoad(6)
 		
 		complex*16 :: CMat(6,6), cTE(6)
 		integer :: secNum
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		expLoad(:) = c_0
 		if(abs(beamExpLoadCoef(1,secNum)) .gt. 0d0) then
 		    expLoad(:) = c_1*beamExpLoadCoef(:,secNum)
 			do i1 = elToDRange(element-1)+1, elToDRange(element)
 			    i2 = elToD(i1)
 				if(dCategory(i2) .eq. 'thermExp') then
 					i3 = dComponent(i2)
 					expLoad(i3) = expLoad(i3) + elToCoef(i1)*c_dVec(i2)
 				endif
 			enddo
 		else
 		    call c_getBeamStiffness(CMat,element)
 			call c_getMaterialThermExp(cTE,element)
 			cTE(2) = cTE(4)
 			cTE(3) = cTE(5)
 			do i1 = 1, 6
 			    do i2 = 1, 3
 				    expLoad(i1) = expLoad(i1) + CMat(i1,i2)*cTE(i2)
 				enddo
 			enddo
 		endif
 		
 	end subroutine c_getBeamExpLoad
 	
 	subroutine c_getBeamMass(mMat,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: mMat(6,6)
 
 		complex*16 :: props(7), secDen, bMass(21)
         integer :: secNum, matId
 		integer :: i1, i2, i3
 		
 		secNum = elementSection(element)
 		matId = sectionMatId(secNum)
 		props(:) = c_1*beamProperties(:,secNum)
 		secDen = c_1*materialDensity(matId)
 		bMass(:) = c_1*beamMass(:,secNum)
 		
 		if(abs(bMass(1)) .gt. 0d0) then
 		    do i1 = elToDRange(element-1)+1, elToDRange(element)
 			    i2 = elToD(i1)
 				if(dCategory(i2) .eq. 'massMat') then
 				    i3 = dComponent(i2)
 					bMass(i3) = bMass(i3) + elToCoef(i1)*c_dVec(i2)
 				endif
 			enddo
 			mMat(1,1:6) = bMass(1:6)
 			mMat(2,2:6) = bMass(7:11)
 			mMat(3,3:6) = bMass(12:15)
 			mMat(4,4:6) = bMass(16:18)
 			mMat(5,5:6) = bMass(19:20)
 			mMat(6,6) = bMass(21)
 			do i1 = 2, 6
 			    do i2 = 1, i1-1
 				    mMat(i1,i2) = mMat(i2,i1)
 				enddo
 			enddo
 		else
 		    do i1 = elToDRange(element-1)+1, elToDRange(element)
 			    i2 = elToD(i1)
 				if(dCategory(i2) .eq. 'area') then
 					props(1) = props(1) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'areaMoment') then
 					i3 = dComponent(i2)
 					props(1+i3) = props(1+i3) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'polarMoment') then
 					props(7) = props(7) + elToCoef(i1)*c_dVec(i2)
 				elseif(dCategory(i2) .eq. 'density') then
 					secDen = secDen + elToCoef(i1)*c_dVec(i2)
 				endif
 			enddo
 			mMat(:,:) = c_0
 			mMat(1,1) = props(1)
 			mMat(1,5) = props(2)
 			mMat(1,6) = -props(3)
 			mMat(2,2) = props(1)
 			mMat(2,4) = -props(2)
 			mMat(3,3) = props(1)
 			mMat(3,4) = props(3)
 			mMat(4,4) = props(7)
 			mMat(5,5) = props(4)
 			mMat(5,6) = -props(6)
 			mMat(6,6) = props(5)
 			mMat = secDen*mMat
 			do i1 = 2, 6
 			    do i2 = 1, i1-1
 				    mMat(i1,i2) = mMat(i2,i1)
 				enddo
 			enddo
 		endif
 	
 	end subroutine c_getBeamMass
 	
 	subroutine c_getBeamTCond(tCond,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: tCond(3,3)
 		
 		complex*16 :: tCVec(6), area
 		integer :: secNum, matId
 		integer :: i1, i2, i3
 		
         secNum = elementSection(element)
         if(abs(beamThermCond(1,secNum)) .ne. 0d0) then
 		    tCVec(:) = c_1*beamThermCond(:,secNum)
 		else
 		    matId = sectionMatId(secNum)
 			tCVec(:) = c_1*materialThermCond(:,matId)
 			area = c_1*beamProperties(1,secNum)
 		endif
 		
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 			i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'thermalCond') then
 				i3 = dComponent(i2)
 				tCVec(i3) = tCVec(i3) + elToCoef(i1)*c_dVec(i2)
 			elseif(dCategory(i2) .eq. 'area') then
 				area = area + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo
 		tCond(1,1) = tCVec(1)
 		tCond(2,2) = tCVec(2)
 		tCond(3,3) = tCVec(3)
 		tCond(1,2) = tCVec(4)
 		tCond(2,1) = tCVec(4)
 		tCond(1,3) = tCVec(5)
 		tCond(3,1) = tCVec(5)
 		tCond(2,3) = tCVec(6)
 		tCond(3,2) = tCVec(6)
 		
 		if(abs(beamThermCond(1,secNum)) .eq. 0d0) then
 		    tCond = area*tCond
 		endif
 		
 	end subroutine c_getBeamTCond
 	
 	subroutine c_getBeamSpecHeat(sH,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: sH
 		
 		complex*16 :: area, den
 		integer :: secNum, matId
 		integer :: i1, i2, i3
 		
         secNum = elementSection(element)
         if(abs(beamSpecHeat(secNum)) .ne. 0d0) then
 		    sH = c_1*beamSpecHeat(secNum)
 		else
 		    matId = sectionMatId(secNum)
 			sH = c_1*materialSpecHeat(matId)
 			den = c_1*materialDensity(matId)
 			area = c_1*beamProperties(1,secNum)
 		endif
 		
 		do i1 = elToDRange(element-1)+1, elToDRange(element)
 			i2 = elToD(i1)
 			if(dCategory(i2) .eq. 'specHeat') then
 				sH = sH + elToCoef(i1)*c_dVec(i2)
 			elseif(dCategory(i2) .eq. 'area') then
 				area = area + elToCoef(i1)*c_dVec(i2)
 			elseif(dCategory(i2) .eq. 'density') then
 				den = den + elToCoef(i1)*c_dVec(i2)
 			endif
 		enddo
 		
 		if(abs(beamSpecHeat(secNum)) .eq. 0d0) then
 		    sH = sH*den*area
 		endif
 	
 	end subroutine c_getBeamSpecHeat
 	
 	subroutine c_getNodeCoord(nodes,element)
 	    implicit none
 		
 		integer, intent(in) :: element
 		complex*16, intent(out) :: nodes(3,10)
 		
 		integer :: i1, i2, i3, i4, i5
 		
 		nodes(:,:) = c_0
 		do i1 = 1,8
 		    i2 = elementList(i1,element)
 			if(i2 .gt. 0) then
 				nodes(:,i1) = c_1*nodeList(:,i2)
 				do i3 = ndToDRange(i2-1)+1, ndToDRange(i2)
 					i4 = ndToD(i3)
 					if(dCategory(i4) .eq. 'nodeCoord') then
 						i5 = dComponent(i4)
 						nodes(i5,i1) = nodes(i5,i1) + ndToCoef(i3)*c_dVec(i4)
 					endif
 				enddo
 			endif
 		enddo
 		
 	end subroutine c_getNodeCoord
 	
 	subroutine c_getNodeDLoad(time,ndLd,node)
 	    implicit none
 		
 		integer, intent(in) :: node
 		complex*16, intent(in) :: time
 		complex*16, intent(out) :: ndLd(7)
 		
 		integer :: i1, i2, i3, i4, i5, i6, check1, check2, nDof
 		
 		ndLd(:) = c_0
 		do i2 = ndToDRange(node-1)+1, ndToDRange(node)
 			i3 = ndToD(i2)
 			call c_greater(check1,time,c_1*dActTime(1,i3))
 			call c_greater(check2,c_1*dActTime(2,i3),time)
 			if(check1 .eq. 1 .and. check2 .eq. 1) then
 				if(dCategory(i3) .eq. 'thermalLoad') then
 					ndLd(7) = ndLd(7) + ndToCoef(i2)*c_dVec(i3)
 				elseif(dCategory(i3) .eq. 'elasticLoad') then
 					i4 = dComponent(i3)
 	                ndLd(i4) = ndLd(i4) + ndToCoef(i2)*c_dVec(i3)
 				endif
 			endif
 		enddo
 		
 	end subroutine c_getNodeDLoad
 	
 	subroutine c_GetInvLU(Ainv,Amat,fVec,zVec,adim,stRow,endRow)
 	    implicit none
 	
 	    integer, intent(in) :: adim, stRow, endRow
 	    complex*16, intent(out) :: Amat(adim,adim)
 	    complex*16, intent(out) :: Ainv(adim,adim), fVec(adim), zVec(adim)
 	
 	    integer :: i1
 
 	    call c_LUFactor(Amat, adim, stRow, endRow)
 	    do i1 = stRow, endRow
 	        fVec(:) = c_0
 	        fVec(i1) = c_1
 	        call c_SolveLUxb(Ainv(:,i1),Amat,fVec,zVec,adim,stRow,endRow)
 	    enddo
 		
 	    return
     end subroutine c_GetInvLU
 	
 	subroutine c_LUFactor(Amat, Adim, stRow, endRow)
         implicit none
 
 	    integer, intent(in) :: Adim, stRow, endRow
         complex*16, intent(out) :: Amat(Adim,Adim)
 
         complex*16 :: fact
         integer :: i1, i2, i3
 
         do i1 = stRow, endRow
             do i2 = i1, endRow
 			    do i3 = stRow, i1-1
 				    Amat(i1,i2) = Amat(i1,i2) - Amat(i1,i3)*Amat(i3,i2)
 				enddo
 			enddo
 			fact = c_1/Amat(i1,i1)
 			do i2 = i1+1, endRow
 			    do i3 = stRow, i1-1
 				    Amat(i2,i1) = Amat(i2,i1) - Amat(i2,i3)*Amat(i3,i1)
 				enddo
 				Amat(i2,i1) = fact*Amat(i2,i1)
 			enddo
         end do
 
         return
     end subroutine c_LUFactor
 
     !SolveLUxb finds the solution for the linear n X n system LUx = b given the L-U factorization LU, storing the solution in x
     subroutine c_SolveLUxb(x, LU, b, z, n, stRow, endRow)
         implicit none
 
 	    integer, intent(in) :: n, stRow, endRow
         complex*16, intent(in) :: LU(n,n), b(n)
         complex*16, intent(out) :: x(n), z(n)
 
         complex*16 :: fact, dp
 
         integer i1, i2
 
         z(stRow) = b(stRow)
         do i1 = stRow+1, endRow
             z(i1) = b(i1)
             do i2 = stRow, i1-1
                 z(i1) = z(i1) - LU(i1,i2)*z(i2)
             enddo
         end do
 
         x(endRow) = z(endRow)/LU(endRow,endRow)
         do i1 = endRow-1, stRow, -1
             fact = c_1/LU(i1,i1)
             dp = c_0
             do i2 = i1+1, endRow
                 dp = dp + LU(i1,i2)*x(i2)
             enddo
             x(i1) = fact*(z(i1) - dp)
         end do
 
         return
     end subroutine c_SolveLUxb
 	
 	subroutine c_rotateAlpha(aGI,aGL,rot,dv1,dv2)
         implicit none
 
         integer, intent(in) :: dv1, dv2
         complex*16, intent(in) :: aGL(3,3) , rot(3)
         complex*16, intent(out) :: aGI(3,3)
 
         integer :: i1, i2, i3, nd1, var1, nd2, var2
         complex*16 :: rSum(3), rLoc(3), Rmag, nRot(3), sR, cR, a1(3,3), a2(3,3), aTemp(3,3), tVal
         complex*16 :: drot1(3), drSum1(3), drLoc1(3), dRmag1, dnRot1(3), dsR1, dcR1, da11(3,3), da21(3,3), daTemp1(3,3), dtVal1
 		complex*16 :: drot2(3), drSum2(3), drLoc2(3), dRmag2, dnRot2(3), dsR2, dcR2, da12(3,3), da22(3,3), daTemp2(3,3), dtVal2
 		complex*16 :: d2rot(3), d2rSum(3), d2rLoc(3), d2Rmag, d2nRot(3), d2sR, d2cR, d2a1(3,3), d2a2(3,3), d2aTemp(3,3), d2tVal
 
         aGI(:,:) = c_0
 
         if(dv1 + dv2 .eq. 0) then
 		    tVal = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3)
 			call c_sqrt(Rmag,tVal)
 			if(abs(Rmag) .lt. 1e-5) then
 				rLoc = c_0
 				do i1 = 1, 3
 					do i2 = 1, 3
 						rLoc(i1) = rLoc(i1) + aGL(i1,i2)*rot(i2)
 					enddo
 				enddo
 				tVal = c_1o6*(rLoc(3)*rLoc(3) + rLoc(1)*rLoc(1) + rLoc(2)*rLoc(2))
 				aTemp(1,1) = c_1 - 0.5d0*(rLoc(2)*rLoc(2) + rLoc(3)*rLoc(3))
 				aTemp(2,1) = 0.5d0*rLoc(1)*rLoc(2) - rLoc(3) + tVal*rLoc(3)
 				aTemp(3,1) = 0.5d0*rLoc(1)*rLoc(3) + rLoc(2) - tVal*rLoc(2)
 				aTemp(1,2) = 0.5d0*rLoc(1)*rLoc(2) + rLoc(3) - tVal*rLoc(3)
 				aTemp(2,2) = c_1 - 0.5d0*(rLoc(1)*rLoc(1) + rLoc(3)*rLoc(3))
 				aTemp(3,2) = 0.5d0*rLoc(2)*rLoc(3) - rLoc(1) + tVal*rLoc(1)
 				aTemp(1,3) = 0.5d0*rLoc(1)*rLoc(3) - rLoc(2) + tVal*rLoc(2)
 				aTemp(2,3) = 0.5d0*rLoc(2)*rLoc(3) + rLoc(1) - tVal*rLoc(1)
 				aTemp(3,3) = c_1 - 0.5d0*(rLoc(1)*rLoc(1) + rLoc(2)*rLoc(2))
 				do i1 = 1, 3
 					do i2 = 1, 3
 						do i3 = 1, 3
 							aGI(i1,i2) = aGI(i1,i2) + aTemp(i1,i3)*aGL(i3,i2)
 						enddo
 					enddo
 				enddo
 			else
 			    call c_sin(sR,Rmag)
 				call c_cos(cR,Rmag)
 
 				nRot = (c_1/Rmag)*rot
 				a1(1,:) = nRot(:)
 				i1 = 1
 				if(abs(nRot(2)) .lt. abs(nRot(1))) then
 					i1 = 2
 				endif
 				if(abs(nRot(3)) .lt. abs(nRot(i1))) then
 					i1 = 3
 				endif
 				a1(3,i1) = c_0
 				tVal = c_1 - nRot(i1)*nRot(i1)
 				call c_sqrt(a1(2,i1),tVal)
 				do i2 = 1, 3
 					if(i2 .ne. i1) then
 						a1(2,i2) = -a1(1,i1)*a1(1,i2)/a1(2,i1)
 					endif
 				enddo
 				a1(3,1) = a1(1,2)*a1(2,3) - a1(1,3)*a1(2,2)
 				a1(3,2) = a1(1,3)*a1(2,1) - a1(1,1)*a1(2,3)
 				a1(3,3) = a1(1,1)*a1(2,2) - a1(1,2)*a1(2,1)
 
 				a2(1,:) = a1(1,:)
 				a2(2,1) = cR*a1(2,1) - sR*a1(3,1)
 				a2(2,2) = cR*a1(2,2) - sR*a1(3,2)
 				a2(2,3) = cR*a1(2,3) - sR*a1(3,3)
 				a2(3,1) = sR*a1(2,1) + cR*a1(3,1)
 				a2(3,2) = sR*a1(2,2) + cR*a1(3,2)
 				a2(3,3) = sR*a1(2,3) + cR*a1(3,3)
 
 				aTemp(:,:) = c_0
 				do i1 = 1, 3
 					do i2 = 1, 3
 						do i3 = 1, 3
 							aTemp(i1,i2) = aTemp(i1,i2) + a2(i3,i1)*a1(i3,i2)
 						enddo
 					enddo
 				enddo
 
 				aGI(:,:) = c_0
 				do i1 = 1, 3
 					do i2 = 1, 3
 						do i3 = 1, 3
 							aGI(i1,i2) = aGI(i1,i2) + aGL(i1,i3)*aTemp(i3,i2)
 						enddo
 					enddo
 				enddo
 			endif
         elseif(dv1*dv2 .eq. 0) then
             drot1 = c_0
             drot1(dv1) = c_1
 			tVal = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3)
 			call c_sqrt(Rmag,tVal)
             if(abs(Rmag) .lt. 1e-5) then
                 rLoc = c_0
                 drLoc1 = c_0
                 do i1 = 1, 3
                     do i2 = 1, 3
                         rLoc(i1) = rLoc(i1) + aGL(i1,i2)*rot(i2)
                         drLoc1(i1) = drLoc1(i1) + aGL(i1,i2)*drot1(i2)
                     enddo
                 enddo
                 tVal = c_1o6*(rLoc(3)*rLoc(3) + rLoc(1)*rLoc(1) + rLoc(2)*rLoc(2))
                 dtVal1 = c_1o3*(drLoc1(3)*rLoc(3) + drLoc1(1)*rLoc(1) + drLoc1(2)*rLoc(2))
                 daTemp1(1,1) = -(rLoc(2)*drLoc1(2) + rLoc(3)*drLoc1(3))
                 daTemp1(2,1) = 0.5d0*(drLoc1(1)*rLoc(2) + rLoc(1)*drLoc1(2)) - drLoc1(3) + dtVal1*rLoc(3) + tVal*drLoc1(3)
                 daTemp1(3,1) = 0.5d0*(drLoc1(1)*rLoc(3) + rLoc(1)*drLoc1(3)) + drLoc1(2) - dtVal1*rLoc(2) - tVal*drLoc1(2)
                 daTemp1(1,2) = 0.5d0*(drLoc1(1)*rLoc(2) + rLoc(1)*drLoc1(2)) + drLoc1(3) - dtVal1*rLoc(3) - tVal*drLoc1(3)
                 daTemp1(2,2) = -(rLoc(1)*drLoc1(1) + rLoc(3)*drLoc1(3))
                 daTemp1(3,2) = 0.5d0*(drLoc1(2)*rLoc(3) + rLoc(2)*drLoc1(3)) - drLoc1(1) + dtVal1*rLoc(1) + tVal*drLoc1(1)
                 daTemp1(1,3) = 0.5d0*(drLoc1(1)*rLoc(3) + rLoc(1)*drLoc1(3)) - drLoc1(2) + dtVal1*rLoc(2) + tVal*drLoc1(2)
                 daTemp1(2,3) = 0.5d0*(drLoc1(2)*rLoc(3) + rLoc(2)*drLoc1(3)) + drLoc1(1) - dtVal1*rLoc(1) - tVal*drLoc1(1)
                 daTemp1(3,3) = -(drLoc1(1)*rLoc(1) + drLoc1(2)*rLoc(2))
                 !write(*,*) var1, daTemp1
                 do i1 = 1, 3
                     do i2 = 1, 3
                         do i3 = 1, 3
                             aGI(i1,i2) = aGI(i1,i2) + daTemp1(i1,i3)*aGL(i3,i2)
                         enddo
                     enddo
                 enddo
             else
 			    tVal = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3)
 			    call c_sqrt(dRmag1,tVal)
                 dRmag1 = (c_1/dRmag1)*(rot(dv1)*drot1(dv1))
 
 				call c_sin(sR,Rmag)
 				call c_cos(cR,Rmag)
 				dsR1 = cR*dRmag1
 				dcR1 = -sR*dRmag1
 
                 nRot = (c_1/Rmag)*rot
                 dnRot1 = (-c_1/(Rmag*Rmag))*dRmag1*rot + (c_1/Rmag)*drot1
                 a1(1,:) = nRot(:)
                 da11(1,:) = dnRot1(:)
                 i1 = 1
                 if(abs(nRot(2)) .lt. abs(nRot(1))) then
                     i1 = 2
                 endif
                 if(abs(nRot(3)) .lt. abs(nRot(i1))) then
                     i1 = 3
                 endif
                 a1(3,i1) = c_0
                 da11(3,i1) = c_0
 				tVal = c_1 - nRot(i1)*nRot(i1)
 				call c_sqrt(a1(2,i1),tVal)
                 da11(2,i1) = -nRot(i1)*dnRot1(i1)/a1(2,i1)
                 do i2 = 1, 3
                     if(i2 .ne. i1) then
                         a1(2,i2) = -a1(1,i1)*a1(1,i2)/a1(2,i1)
                         da11(2,i2) = -da11(1,i1)*a1(1,i2)/a1(2,i1) - a1(1,i1)*da11(1,i2)/a1(2,i1)
                         da11(2,i2) = da11(2,i2) + a1(1,i1)*a1(1,i2)*da11(2,i1)/(a1(2,i1)*a1(2,i1))
                     endif
                 enddo
                 a1(3,1) = a1(1,2)*a1(2,3) - a1(1,3)*a1(2,2)
                 da11(3,1) = da11(1,2)*a1(2,3) + a1(1,2)*da11(2,3) - da11(1,3)*a1(2,2) - a1(1,3)*da11(2,2)
                 a1(3,2) = a1(1,3)*a1(2,1) - a1(1,1)*a1(2,3)
                 da11(3,2) = da11(1,3)*a1(2,1) + a1(1,3)*da11(2,1) - da11(1,1)*a1(2,3) - a1(1,1)*da11(2,3)
                 a1(3,3) = a1(1,1)*a1(2,2) - a1(1,2)*a1(2,1)
                 da11(3,3) = da11(1,1)*a1(2,2) + a1(1,1)*da11(2,2) - da11(1,2)*a1(2,1) - a1(1,2)*da11(2,1)
 
                 a2(1,:) = a1(1,:)
                 da21(1,:) = da11(1,:)
                 a2(2,1) = cR*a1(2,1) - sR*a1(3,1)
                 da21(2,1) = dcR1*a1(2,1) + cR*da11(2,1) - dsR1*a1(3,1) - sR*da11(3,1)
                 a2(2,2) = cR*a1(2,2) - sR*a1(3,2)
                 da21(2,2) = dcR1*a1(2,2) + cR*da11(2,2) - dsR1*a1(3,2) - sR*da11(3,2)
                 a2(2,3) = cR*a1(2,3) - sR*a1(3,3)
                 da21(2,3) = dcR1*a1(2,3) + cR*da11(2,3) - dsR1*a1(3,3) - sR*da11(3,3)
                 a2(3,1) = sR*a1(2,1) + cR*a1(3,1)
                 da21(3,1) = dsR1*a1(2,1) + sR*da11(2,1) + dcR1*a1(3,1) + cR*da11(3,1)
                 a2(3,2) = sR*a1(2,2) + cR*a1(3,2)
                 da21(3,2) = dsR1*a1(2,2) + sR*da11(2,2) + dcR1*a1(3,2) + cR*da11(3,2)
                 a2(3,3) = sR*a1(2,3) + cR*a1(3,3)
                 da21(3,3) = dsR1*a1(2,3) + sR*da11(2,3) + dcR1*a1(3,3) + cR*da11(3,3)
 
                 daTemp1(:,:) = c_0
                 do i1 = 1, 3
                     do i2 = 1, 3
                         do i3 = 1, 3
                             daTemp1(i1,i2) = daTemp1(i1,i2) + da21(i3,i1)*a1(i3,i2) + a2(i3,i1)*da11(i3,i2)
                         enddo
                     enddo
                 enddo
 
                 aGI(:,:) = c_0
                 do i1 = 1, 3
                     do i2 = 1, 3
                         do i3 = 1, 3
                             aGI(i1,i2) = aGI(i1,i2) + aGL(i1,i3)*daTemp1(i3,i2)
                         enddo
                     enddo
                 enddo
             endif
         else
             drot1 = c_0
             drot2 = c_0
             drot1(dv1) = c_1
             drot2(dv2) = c_1
 			tVal = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3)
             call c_sqrt(Rmag,tVal) !!Rmag = (rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3))**0.5d0
             if(abs(Rmag) .lt. 1e-5) then
                 rLoc = c_0
                 drLoc1 = c_0
                 drLoc2 = c_0
                 do i1 = 1, 3
                     do i2 = 1, 3
                         rLoc(i1) = rLoc(i1) + aGL(i1,i2)*rot(i2)
                         drLoc1(i1) = drLoc1(i1) + aGL(i1,i2)*drot1(i2)
                         drLoc2(i1) = drLoc2(i1) + aGL(i1,i2)*drot2(i2)
                     enddo
                 enddo
                 tVal = c_1o6*(rLoc(3)*rLoc(3) + rLoc(1)*rLoc(1) + rLoc(2)*rLoc(2))
                 dtVal1 = c_1o3*(drLoc1(3)*rLoc(3) + drLoc1(1)*rLoc(1) + drLoc1(2)*rLoc(2))
                 dtVal2 = c_1o3*(drLoc2(3)*rLoc(3) + drLoc2(1)*rLoc(1) + drLoc2(2)*rLoc(2))
                 d2tVal = c_1o3*(drLoc1(3)*drLoc2(3) + drLoc1(1)*drLoc2(1) + drLoc1(2)*drLoc2(2))
                 d2aTemp(1,1) = -(drLoc2(2)*drLoc1(2) + drLoc2(3)*drLoc1(3))
                 d2aTemp(2,1) = 0.5d0*(drLoc1(1)*drLoc2(2) + drLoc2(1)*drLoc1(2)) + &
                                d2tVal*rLoc(3) + dtVal1*drLoc2(3) + dtVal2*drLoc1(3)
                 d2aTemp(3,1) = 0.5d0*(drLoc1(1)*drLoc2(3) + drLoc2(1)*drLoc1(3)) - &
                                d2tVal*rLoc(2) - dtVal1*drLoc2(2) - dtVal2*drLoc1(2)
                 d2aTemp(1,2) = 0.5d0*(drLoc1(1)*drLoc2(2) + drLoc2(1)*drLoc1(2)) - &
                                d2tVal*rLoc(3) - dtVal1*drLoc2(3) - dtVal2*drLoc1(3)
                 d2aTemp(2,2) = -(drLoc2(1)*drLoc1(1) + drLoc2(3)*drLoc1(3))
                 d2aTemp(3,2) = 0.5d0*(drLoc1(2)*drLoc2(3) + drLoc2(2)*drLoc1(3)) + &
                                d2tVal*rLoc(1) + dtVal1*drLoc2(1) + dtVal2*drLoc1(1)
                 d2aTemp(1,3) = 0.5d0*(drLoc1(1)*drLoc2(3) + drLoc2(1)*drLoc1(3)) + &
                                d2tVal*rLoc(2) + dtVal1*drLoc2(2) + dtVal2*drLoc1(2)
                 d2aTemp(2,3) = 0.5d0*(drLoc1(2)*drLoc2(3) + drLoc2(2)*drLoc1(3)) - &
                                d2tVal*rLoc(1) + dtVal1*drLoc2(1) + dtVal2*drLoc1(1)
                 d2aTemp(3,3) = -(drLoc1(1)*drLoc2(1) + drLoc1(2)*drLoc2(2))
                 do i1 = 1, 3
                     do i2 = 1, 3
                         do i3 = 1, 3
                             aGI(i1,i2) = aGI(i1,i2) + d2aTemp(i1,i3)*aGL(i3,i2)
                         enddo
                     enddo
                 enddo
             else
                 dRmag1 = c_1/Rmag !!(rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3))**(-0.5d0)
                 dRmag2 = dRmag1
                 if(dv1 .eq. dv2) then
                     d2Rmag = dRmag1*drot1(dv1)*drot2(dv1)
                 else
                     d2Rmag = c_0
                 endif
                 d2Rmag = d2Rmag - (dRmag1*dRmag1*dRmag1)*rot(dv1)*drot1(dv1)*rot(dv2)*drot2(dv2) !!d2Rmag = d2Rmag - dRmag1**3d0*rot(dv1)*drot1(dv1)*rot(dv2)*drot2(dv2)
                 dRmag1 = dRmag1*(rot(dv1)*drot1(dv1))
                 dRmag2 = dRmag2*(rot(dv2)*drot2(dv2))
                 call c_sin(sR,Rmag) !! sR = sin(Rmag)
 				call c_cos(cR,Rmag)
                 dsR1 = cR*dRmag1
                 dsR2 = cR*dRmag2
                 d2sR = -sR*dRmag2*dRmag1 + cR*d2Rmag
                 !!cR = cos(Rmag)
                 dcR1 = -sR*dRmag1!!-sin(Rmag)*dRmag1
                 dcR2 = -sR*dRmag2 !!-sin(Rmag)*dRmag2
                 d2cR = -cR*dRmag2*dRmag1 - sR*d2Rmag !!-cos(Rmag)*dRmag2*dRmag1 - sin(Rmag)*d2Rmag
 
                 nRot = (c_1/Rmag)*rot
                 dnRot1 = (-c_1/(Rmag*Rmag))*dRmag1*rot + (c_1/Rmag)*drot1
                 dnRot2 = (-c_1/(Rmag*Rmag))*dRmag2*rot + (c_1/Rmag)*drot2
                 d2nRot = (2d0*dRmag1*dRmag2/Rmag)*rot - d2Rmag*rot - dRmag1*dRot2 - dRmag2*dRot1
                 d2nRot = (c_1/(Rmag*Rmag))*d2nRot
                 a1(1,:) = nRot(:)
                 da11(1,:) = dnRot1(:)
                 da12(1,:) = dnRot2(:)
                 d2a1(1,:) = d2nRot(:)
                 i1 = 1
                 if(abs(nRot(2)) .lt. abs(nRot(1))) then
                     i1 = 2
                 endif
                 if(abs(nRot(3)) .lt. abs(nRot(i1))) then
                     i1 = 3
                 endif
                 a1(3,i1) = c_0
                 da11(3,i1) = c_0
                 da12(3,i1) = c_0
                 d2a1(3,i1) = c_0
                 !!a1(2,i1) = (c_1 - nRot(i1)*nRot(i1))**0.5d0
 				tVal = c_1 - nRot(i1)*nRot(i1)
 				call c_sqrt(a1(2,i1),tVal)
                 !!da11(2,i1) = -nRot(i1)*dnRot1(i1)*(c_1 - nRot(i1)*nRot(i1))**(-0.5d0)
 				da11(2,i1) = -nRot(i1)*dnRot1(i1)*(c_1/a1(2,i1))
                 !!da12(2,i1) = -nRot(i1)*dnRot2(i1)*(c_1 - nRot(i1)*nRot(i1))**(-0.5d0)
 				da12(2,i1) = -nRot(i1)*dnRot2(i1)*(c_1/a1(2,i1))
                 !!d2a1(2,i1) = -(dnRot2(i1)*dnRot1(i1) + nRot(i1)*d2nRot(i1))*(c_1 - nRot(i1)*nRot(i1))**(-0.5d0)
 				d2a1(2,i1) = -(dnRot2(i1)*dnRot1(i1) + nRot(i1)*d2nRot(i1))*(c_1/a1(2,i1))
 				tVal = (c_1/(a1(2,i1)*a1(2,i1)*a1(2,i1)))
                 !! d2a1(2,i1) = d2a1(2,i1) - nRot(i1)*dnRot1(i1)*nRot(i1)*dnRot2(i1)*(c_1 - nRot(i1)*nRot(i1))**(-1.5d0)
 				d2a1(2,i1) = d2a1(2,i1) - nRot(i1)*dnRot1(i1)*nRot(i1)*dnRot2(i1)*tVal
                 do i2 = 1, 3
                     if(i2 .ne. i1) then
                         a1(2,i2) = -a1(1,i1)*a1(1,i2)/a1(2,i1)
 
                         da11(2,i2) = -da11(1,i1)*a1(1,i2)/a1(2,i1) - a1(1,i1)*da11(1,i2)/a1(2,i1)
                         da11(2,i2) = da11(2,i2) + a1(1,i1)*a1(1,i2)*da11(2,i1)/(a1(2,i1)*a1(2,i1))
                         da12(2,i2) = -da12(1,i1)*a1(1,i2)/a1(2,i1) - a1(1,i1)*da12(1,i2)/a1(2,i1)
                         da12(2,i2) = da12(2,i2) + a1(1,i1)*a1(1,i2)*da12(2,i1)/(a1(2,i1)*a1(2,i1))
 
                         d2a1(2,i2) = -(d2a1(1,i1)*a1(1,i2) + da11(1,i1)*da12(1,i2))/a1(2,i1)
                         d2a1(2,i2) = d2a1(2,i2) - (da12(1,i1)*da11(1,i2) + a1(1,i1)*d2a1(1,i2))/a1(2,i1)
                         d2a1(2,i2) = d2a1(2,i2) + da11(1,i1)*a1(1,i2)*da12(2,i1)/(a1(2,i1)*a1(2,i1))
                         d2a1(2,i2) = d2a1(2,i2) + a1(1,i1)*da11(1,i2)*da12(2,i1)/(a1(2,i1)*a1(2,i1))
                         d2a1(2,i2) = d2a1(2,i2) + da12(1,i1)*a1(1,i2)*da11(2,i1)/(a1(2,i1)*a1(2,i1))
                         d2a1(2,i2) = d2a1(2,i2) + a1(1,i1)*da12(1,i2)*da11(2,i1)/(a1(2,i1)*a1(2,i1))
                         d2a1(2,i2) = d2a1(2,i2) + a1(1,i1)*a1(1,i2)*d2a1(2,i1)/(a1(2,i1)*a1(2,i1))
                         d2a1(2,i2) = d2a1(2,i2) - 2d0*a1(1,i1)*a1(1,i2)*da11(2,i1)*da12(2,i1)/(a1(2,i1)*a1(2,i1)*a1(2,i1))
                     endif
                 enddo
                 a1(3,1) = a1(1,2)*a1(2,3) - a1(1,3)*a1(2,2)
                 da11(3,1) = da11(1,2)*a1(2,3) + a1(1,2)*da11(2,3) - da11(1,3)*a1(2,2) - a1(1,3)*da11(2,2)
                 da12(3,1) = da12(1,2)*a1(2,3) + a1(1,2)*da12(2,3) - da12(1,3)*a1(2,2) - a1(1,3)*da12(2,2)
                 d2a1(3,1) = d2a1(1,2)*a1(2,3) + da11(1,2)*da12(2,3) + da12(1,2)*da11(2,3) + a1(1,2)*d2a1(2,3)
                 d2a1(3,1) = d2a1(3,1) - d2a1(1,3)*a1(2,2) - da11(1,3)*da12(2,2) - da12(1,3)*da11(2,2) - a1(1,3)*d2a1(2,2)
                 a1(3,2) = a1(1,3)*a1(2,1) - a1(1,1)*a1(2,3)
                 da11(3,2) = da11(1,3)*a1(2,1) + a1(1,3)*da11(2,1) - da11(1,1)*a1(2,3) - a1(1,1)*da11(2,3)
                 da12(3,2) = da12(1,3)*a1(2,1) + a1(1,3)*da12(2,1) - da12(1,1)*a1(2,3) - a1(1,1)*da12(2,3)
                 d2a1(3,2) = d2a1(1,3)*a1(2,1) + da11(1,3)*da12(2,1) + da12(1,3)*da11(2,1) + a1(1,3)*d2a1(2,1)
                 d2a1(3,2) = d2a1(3,2) - d2a1(1,1)*a1(2,3) - da11(1,1)*da12(2,3) - da12(1,1)*da11(2,3) - a1(1,1)*d2a1(2,3)
                 a1(3,3) = a1(1,1)*a1(2,2) - a1(1,2)*a1(2,1)
                 da11(3,3) = da11(1,1)*a1(2,2) + a1(1,1)*da11(2,2) - da11(1,2)*a1(2,1) - a1(1,2)*da11(2,1)
                 da12(3,3) = da12(1,1)*a1(2,2) + a1(1,1)*da12(2,2) - da12(1,2)*a1(2,1) - a1(1,2)*da12(2,1)
                 d2a1(3,3) = d2a1(1,1)*a1(2,2) + da11(1,1)*da12(2,2) + da12(1,1)*da11(2,2) + a1(1,1)*d2a1(2,2)
                 d2a1(3,3) = d2a1(3,3) - d2a1(1,2)*a1(2,1) - da11(1,2)*da12(2,1) - da12(1,2)*da11(2,1) - a1(1,2)*d2a1(2,1)
 
                 a2(1,:) = a1(1,:)
                 da21(1,:) = da11(1,:)
                 da22(1,:) = da12(1,:)
                 d2a2(1,:) = d2a1(1,:)
                 a2(2,1) = cR*a1(2,1) - sR*a1(3,1)
                 da21(2,1) = dcR1*a1(2,1) + cR*da11(2,1) - dsR1*a1(3,1) - sR*da11(3,1)
                 da22(2,1) = dcR2*a1(2,1) + cR*da12(2,1) - dsR2*a1(3,1) - sR*da12(3,1)
                 d2a2(2,1) = d2cR*a1(2,1) + dcR1*da12(2,1) + dcR2*da11(2,1) + cR*d2a1(2,1)
                 d2a2(2,1) = d2a2(2,1) - d2sR*a1(3,1) - dsR1*da12(3,1) - dsR2*da11(3,1) - sR*d2a1(3,1)
                 a2(2,2) = cR*a1(2,2) - sR*a1(3,2)
                 da21(2,2) = dcR1*a1(2,2) + cR*da11(2,2) - dsR1*a1(3,2) - sR*da11(3,2)
                 da22(2,2) = dcR2*a1(2,2) + cR*da12(2,2) - dsR2*a1(3,2) - sR*da12(3,2)
                 d2a2(2,2) = d2cR*a1(2,2) + dcR1*da12(2,2) + dcR2*da11(2,2) + cR*d2a1(2,2)
                 d2a2(2,2) = d2a2(2,2) - d2sR*a1(3,2) - dsR1*da12(3,2) - dsR2*da11(3,2) - sR*d2a1(3,2)
                 a2(2,3) = cR*a1(2,3) - sR*a1(3,3)
                 da21(2,3) = dcR1*a1(2,3) + cR*da11(2,3) - dsR1*a1(3,3) - sR*da11(3,3)
                 da22(2,3) = dcR2*a1(2,3) + cR*da12(2,3) - dsR2*a1(3,3) - sR*da12(3,3)
                 d2a2(2,3) = d2cR*a1(2,3) + dcR1*da12(2,3) + dcR2*da11(2,3) + cR*d2a1(2,3)
                 d2a2(2,3) = d2a2(2,3) - d2sR*a1(3,3) - dsR1*da12(3,3) - dsR2*da11(3,3) - sR*d2a1(3,3)
                 a2(3,1) = sR*a1(2,1) + cR*a1(3,1)
                 da21(3,1) = dsR1*a1(2,1) + sR*da11(2,1) + dcR1*a1(3,1) + cR*da11(3,1)
                 da22(3,1) = dsR2*a1(2,1) + sR*da12(2,1) + dcR2*a1(3,1) + cR*da12(3,1)
                 d2a2(3,1) = d2sR*a1(2,1) + dsR1*da12(2,1) + dsR2*da11(2,1) + sR*d2a1(2,1)
                 d2a2(3,1) = d2a2(3,1) + d2cR*a1(3,1) + dcR1*da12(3,1) + dcR2*da11(3,1) + cR*d2a1(3,1)
                 a2(3,2) = sR*a1(2,2) + cR*a1(3,2)
                 da21(3,2) = dsR1*a1(2,2) + sR*da11(2,2) + dcR1*a1(3,2) + cR*da11(3,2)
                 da22(3,2) = dsR2*a1(2,2) + sR*da12(2,2) + dcR2*a1(3,2) + cR*da12(3,2)
                 d2a2(3,2) = d2sR*a1(2,2) + dsR1*da12(2,2) + dsR2*da11(2,2) + sR*d2a1(2,2)
                 d2a2(3,2) = d2a2(3,2) + d2cR*a1(3,2) + dcR1*da12(3,2) + dcR2*da11(3,2) + cR*d2a1(3,2)
                 a2(3,3) = sR*a1(2,3) + cR*a1(3,3)
                 da21(3,3) = dsR1*a1(2,3) + sR*da11(2,3) + dcR1*a1(3,3) + cR*da11(3,3)
                 da22(3,3) = dsR2*a1(2,3) + sR*da12(2,3) + dcR2*a1(3,3) + cR*da12(3,3)
                 d2a2(3,3) = d2sR*a1(2,3) + dsR1*da12(2,3) + dsR2*da11(2,3) + sR*d2a1(2,3)
                 d2a2(3,3) = d2a2(3,3) + d2cR*a1(3,3) + dcR1*da12(3,3) + dcR2*da11(3,3) + cR*d2a1(3,3)
 
                 d2aTemp(:,:) = c_0
                 do i1 = 1, 3
                     do i2 = 1, 3
                         do i3 = 1, 3
                             d2aTemp(i1,i2) = d2aTemp(i1,i2) + d2a2(i3,i1)*a1(i3,i2) + da21(i3,i1)*da12(i3,i2)
                             d2aTemp(i1,i2) = d2aTemp(i1,i2) + da22(i3,i1)*da11(i3,i2) + a2(i3,i1)*d2a1(i3,i2)
                         enddo
                     enddo
                 enddo
 
                 aGI(:,:) = c_0
                 do i1 = 1, 3
                     do i2 = 1, 3
                         do i3 = 1, 3
                             aGI(i1,i2) = aGI(i1,i2) + aGL(i1,i3)*d2aTemp(i3,i2)
                         enddo
                     enddo
                 enddo
             endif          		
         endif
 
         return
     end subroutine c_rotateAlpha
 	
 end module ASenD_c_designPropertyFunctions
