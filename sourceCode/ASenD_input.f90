module ASenD_input
    use ASenD_globalData
	use ASenD_constantVals
	
	contains
	
	subroutine findSubstrings(found,string,subStrings,numSub)
	    implicit none
		
		integer, intent(out) :: found
		integer, intent(in) :: numSub
		character(len=256), intent(in) :: string
		character(len=16), intent(in) :: subStrings(numSub)
		
		integer :: i1
		
		found = 0
		do i1 = 1, numSub
		    found = index(string,subStrings(i1))
			if(found .ne. 0) then
			    return
			endif
		enddo
		
	end subroutine findSubstrings
	
	subroutine readModelInput(inFileName)
	    implicit none
		
		character(len=128), intent(in) :: inFileName
		
		character(len=256) :: fileLine
		character(len=16) :: subStrings(16)
		integer :: iosVal
		integer :: i1, i2, i3, i4, i5, i6, found, inElSet, inMaxStress
		integer :: elType, numNds, secNum, matNum, GSet, nuSet
		integer :: layerNum, termNum, mpcNum, loadNum, eSNum, nSNum, elSi, ndSi
		real*8 :: mag
		real*8 :: v1(3), v2(3), v3(3)
        integer :: readInt(9)
        real*8 :: readReal(21), readMat(6,6)
        character(len=64) :: readChar		
		
		open(unit=1, file=inFileName, blank='NULL', action='read')
		
		subStrings(1) = 'nodes:'
		subStrings(2) = 'elements:'
		subStrings(3) = 'sets:'
		subStrings(4) = 'sections:'
		subStrings(5) = 'materials:'
		subStrings(6) = 'constraints:'
		subStrings(7) = 'loads:'
		subStrings(8) = 'initialState:'
		
		fileLine(1:15) = '               '
		numNodes = 0
		numEls = 0
		numNdSets = 0
		ndSetSize = 0
		numElSets = 0
		elSetSize = 0
		inElSet = 0
		numSec = 0
		layupSize = 0
		numMats = 0
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)  !! IOSTAT = 0 perfect, > 0 problem occured, < 0 end of file
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
		    if(i1 .gt. 0) then
			    if(fileLine(i1-5:i1) .eq. 'nodes:') then
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)  !while fileLine does not have a colon
						i2 = index(fileLine,'-')
						if(i2 .gt. 0) then
						    numNodes = numNodes + 1
						endif
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
						i1 = index(fileLine,':')
					enddo
				elseif(fileLine(i1-12:i1) .eq. 'connectivity:') then
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)  !while fileLine does not have a colon
						i2 = index(fileLine,'-')
						if(i2 .ne. 0) then
						    numEls = numEls + 1
						endif
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
						i1 = index(fileLine,':')
					enddo
				elseif(fileLine(i1-4:i1) .eq. 'node:') then
				    inElSet = 0
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,':')
				elseif(fileLine(i1-7:i1) .eq. 'element:') then
				    inElSet = 1
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,':')
				elseif(fileLine(i1-6:i1) .eq. 'labels:') then
				    if(inElSet .eq. 1) then
					    numElSets = numElSets + 1
					    i2 = index(fileLine(i1+1:i1+64),'all')
						if(i2 .gt. 1) then
						    elSetSize = elSetSize + numEls
							read(1,'(A)',iostat=iosVal) fileLine(16:256)
					        i1 = index(fileLine,':')
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					        i1 = index(fileLine,':')
							do while(i1 .eq. 0 .and. iosVal .eq. 0)
							    i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									if(i3 .gt. 0) then
									    readInt(1:3) = (/1,0,1/)
										read(fileLine(i2+1:i3-1),*,end=175) readInt(1:3)
175                                     do i4 = readInt(1), readInt(2), readInt(3)
                                            elSetSize = elSetSize + 1
                                        enddo
									endif
								else
								    i2 = index(fileLine,'-')
									if(i2 .gt. 0) then
									    elSetSize = elSetSize + 1
									endif
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
					            i1 = index(fileLine,':')
							enddo
						endif
					else
					    numNdSets = numNdSets + 1
					    i2 = index(fileLine(i1+1:i1+64),'all')
						if(i2 .gt. 1) then
						    ndSetSize = ndSetSize + numNodes
							read(1,'(A)',iostat=iosVal) fileLine(16:256)
					        i1 = index(fileLine,':')
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					        i1 = index(fileLine,':')
							do while(i1 .eq. 0 .and. iosVal .eq. 0)
							    i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									if(i3 .gt. 0) then
									    readInt(1:3) = (/1,0,1/)
										read(fileLine(i2+1:i3-1),*,end=207) readInt(1:3)
207                                     do i4 = readInt(1), readInt(2), readInt(3)
                                            ndSetSize = ndSetSize + 1
                                        enddo
									endif
								else
								    i2 = index(fileLine,'-')
									if(i2 .gt. 0) then
									    ndSetSize = ndSetSize + 1
									endif
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
					            i1 = index(fileLine,':')
							enddo
						endif					
					endif
				elseif(fileLine(i1-8:i1) .eq. 'sections:') then
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
						i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-4:i1) .eq. 'type:') then
							    numSec = numSec + 1
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-6:i1) .eq. 'layers:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'-')
									if(i2 .gt. 0) then
									    layupSize = layupSize + 1
									endif
								    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								    i1 = index(fileLine,':')
								enddo
							else
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				elseif(fileLine(i1-9:i1) .eq. 'materials:') then
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
						i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-4:i1) .eq. 'name:') then
							    numMats = numMats + 1
							endif
						endif
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
		enddo
		
		if(.not. allocated(nodeList)) then
		    allocate(nodeList(3,numNodes))
			allocate(elementList(8,numEls))
			allocate(elementType(numEls))
			allocate(elementSurfaces(6,numEls))
			allocate(elementSection(numEls))
			allocate(elDepends(numEls))
			allocate(elementSets(elSetSize))
			allocate(elSetName(numElSets))
			allocate(elSetRange(0:numElSets))
			allocate(nodeSets(ndSetSize))
			allocate(ndSetName(numNdSets))
			allocate(ndSetRange(0:numNdSets))
			allocate(elToDRange(0:numEls))
			allocate(ndToDRange(0:numNodes))
			allocate(sectionType(numSec))
			allocate(sectionMatId(numSec))
			allocate(sectionMatName(numSec))
			allocate(sectionOrient(3,3*numSec))
			allocate(secLayupRange(0:numSec))
			allocate(sectionZOffset(numSec))
			allocate(beamProperties(7,numSec))
			allocate(beamStiffness(21,numSec))
			allocate(beamExpLoadCoef(6,numSec))
			allocate(beamMass(21,numSec))
			allocate(beamThermCond(6,numSec))
			allocate(beamSpecHeat(numSec))
			allocate(initialDisp(6,numNodes))
			allocate(initialVel(6,numNodes))
			allocate(initialAcc(6,numNodes))
			allocate(initialTemp(numNodes))
			allocate(initialTdot(numNodes))
			
			nodeList(:,:) = r_0
			elementList(:,:) = r_0
			elementType(:) = r_0
			elementSection(:) = r_0
			elementSets(:) = r_0
			elSetRange(:) = r_0
			nodeSets(:) = r_0
			ndSetRange(:) = r_0
			elToDRange(:) = r_0
			ndToDRange(:) = r_0
			sectionMatId(:) = r_0
			sectionOrient(:,:) = r_0
			secLayupRange(:) = r_0
			sectionZOffset(:) = r_0
			beamProperties(:,:) = r_0
			beamStiffness(:,:) = r_0
			beamExpLoadCoef(:,:) = r_0
			beamMass(:,:) = r_0
			beamThermCond(:,:) = r_0
			beamSpecHeat(:) = r_0
			initialDisp(:,:) = r_0
			initialVel(:,:) = r_0
			initialAcc(:,:) = r_0
			initialTemp(:) = r_0
			initialTdot(:) = r_0
		endif
		
		if(layupSize .gt. 0 .and. .not. allocated(layupMatName)) then
			allocate(layupMatName(layupSize))
			allocate(layupMatID(layupSize))
			allocate(layupThickness(layupSize))
			allocate(layupAngle(layupSize))
			
			layupMatID(:) = r_0
			layupThickness(:) = r_0
			layupAngle(:) = r_0
		endif
		
		if(numMats .gt. 0 .and. .not. allocated(materialName)) then
		    allocate(materialName(numMats))
			allocate(materialDensity(numMats))
			allocate(materialElastic(9,numMats))
			allocate(materialStiffMat(21,numMats))
			allocate(materialThermCond(6,numMats))
			allocate(materialThermExp(6,numMats))
			allocate(materialSpecHeat(numMats))
			allocate(materialMaxStress(9,numMats))
			allocate(materialMaxStrain(9,numMats))
			allocate(materialMaxStrnEngy(numMats))
			
		    materialDensity(:) = r_0
			materialElastic(:,:) = r_0
			materialStiffMat(:,:) = r_0
			materialThermCond(:,:) = r_0
			materialThermExp(:,:) = r_0
			materialSpecHeat(:) = r_0
			materialMaxStress(:,:) = r_0
			materialMaxStrain(:,:) = r_0
			materialMaxStrnEngy(:) = r_0
		endif
		
		rewind(1)
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)  !! IOSTAT = 0 perfect, > 0 problem occured, < 0 end of file
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
		    if(i1 .gt. 0) then
			    if(fileLine(i1-5:i1) .eq. 'nodes:') then
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)  !while fileLine does not have a colon
						i2 = index(fileLine,'-')
						if(i2 .gt. 0) then
						    i3 = index(fileLine,'[')
							i4 = index(fileLine,']')
							read(fileLine(i3+1:i4-1),*) readInt(1), readReal(1:3)
							if(readInt(1) .le. numNodes) then
							    nodeList(:,readInt(1)) = readReal(1:3)
							endif
						endif
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
						i1 = index(fileLine,':')
					enddo
				elseif(fileLine(i1-8:i1) .eq. 'elements:') then
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
					    i1 = index(fileLine,':')
					    if(i1 .gt. 0) then
						    if(fileLine(i1-4:i1) .eq. 'type:') then
							    readChar = fileLine(i1+1:i1+64)
								if(index(readChar,'tet4') .gt. 0) then
								    numNds = 4
									elType = 4
								elseif(index(readChar,'wedge6') .gt. 0) then
								    numNds = 6
									elType = 6
								elseif(index(readChar,'brick8') .gt. 0) then
								    numNds = 8
									elType = 8
								elseif(index(readChar,'brickIM') .gt. 0) then
								    numNds = 8
									elType = 81
								elseif(index(readChar,'shell3') .gt. 0) then
								    numNds = 3
									elType = 3
								elseif(index(readChar,'shell4') .gt. 0) then
								    numNds = 4
									elType = 41
								elseif(index(readChar,'beam2') .gt. 0) then
								    numNds = 2
									elType = 2
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
							elseif(fileLine(i1-12:i1) .eq. 'connectivity:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
										read(fileLine(i2+1:i3-1),*) readInt(1:numNds+1)
										if(readInt(1) .le. numEls) then
										    elementList(1:numNds,readInt(1)) = readInt(2:numNds+1)
											elementType(readInt(1)) = elType
										endif
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
								    i1 = index(fileLine,':')
								enddo
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				elseif(fileLine(i1-4:i1) .eq. 'sets:') then
				    nSNum = 0
					eSNum = 0
				    ndSi = 0
					elSi = 0
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
					    i1 = index(fileLine,':')
						if(i1 .gt. 0) then
							if(fileLine(i1-4:i1) .eq. 'node:') then
								inElSet = 0
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
								call findSubstrings(found,fileLine,subStrings,8)
							elseif(fileLine(i1-7:i1) .eq. 'element:') then
								inElSet = 1
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
								call findSubstrings(found,fileLine,subStrings,8)
							elseif(fileLine(i1-4:i1) .eq. 'name:') then
								if(inElSet .eq. 1) then
									eSNum = eSNum + 1
									read(fileLine(i1+1:i1+64),*) elSetName(esNum)
								else
									nSNum = nSNum + 1
									read(fileLine(i1+1:i1+64),*) ndSetName(nsNum)
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
								call findSubstrings(found,fileLine,subStrings,8)
							elseif(fileLine(i1-6:i1) .eq. 'labels:') then
								if(inElSet .eq. 1) then
									i2 = index(fileLine(i1+1:i1+64),'all')
									if(i2 .gt. 1) then
										do i3 = 1, numEls
											elSi = elSi + 1
											elementSets(elSi) = i3
										enddo
										read(1,'(A)',iostat=iosVal) fileLine(16:256)
										call findSubstrings(found,fileLine,subStrings,8)
									else
										read(1,'(A)',iostat=iosVal) fileLine(16:256)
										i1 = index(fileLine,':')
										do while(i1 .eq. 0 .and. iosVal .eq. 0)
											i2 = index(fileLine,'[')
											if(i2 .gt. 0) then
												i3 = index(fileLine,']')
												if(i3 .gt. 0) then
													readInt(1:3) = (/1,0,1/)
													read(fileLine(i2+1:i3-1),*,end=516) readInt(1:3)
516                                                 do i4 = readInt(1), readInt(2), readInt(3)
														elSi = elSi + 1
														elementSets(elSi) = i4
													enddo
												endif
											else
												i2 = index(fileLine,'-')
												if(i2 .gt. 0) then
													elSi = elSi + 1
													read(fileLine(i2+1:i2+64),*) elementSets(elSi)
												endif
											endif
											read(1,'(A)',iostat=iosVal) fileLine(16:256)
											i1 = index(fileLine,':')
										enddo
									endif
									elSetRange(eSNum) = elSi
									call findSubstrings(found,fileLine,subStrings,8)
								else
									i2 = index(fileLine(i1+1:i1+64),'all')
									if(i2 .gt. 1) then
										do i3 = 1, numNodes
											ndSi = ndSi + 1
											nodeSets(ndSi) = i3
										enddo
										read(1,'(A)',iostat=iosVal) fileLine(16:256)
										call findSubstrings(found,fileLine,subStrings,8)
									else
										read(1,'(A)',iostat=iosVal) fileLine(16:256)
										i1 = index(fileLine,':')
										do while(i1 .eq. 0 .and. iosVal .eq. 0)
											i2 = index(fileLine,'[')
											if(i2 .gt. 0) then
												i3 = index(fileLine,']')
												if(i3 .gt. 0) then
													readInt(1:3) = (/1,0,1/)
													read(fileLine(i2+1:i3-1),*,end=553) readInt(1:3)
553                                                 do i4 = readInt(1), readInt(2), readInt(3)
														ndSi = ndSi + 1
														nodeSets(ndSi) = i4
													enddo
												endif
											else
												i2 = index(fileLine,'-')
												if(i2 .gt. 0) then
													ndSi = ndSi + 1
													read(fileLine(i2+1:i2+64),*) nodeSets(ndSi)
												endif
											endif
											read(1,'(A)',iostat=iosVal) fileLine(16:256)
											i1 = index(fileLine,':')
										enddo
									endif
									ndSetRange(nSNum) = ndSi
									call findSubstrings(found,fileLine,subStrings,8)	
								endif
							else
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
								call findSubstrings(found,fileLine,subStrings,8)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						    call findSubstrings(found,fileLine,subStrings,8)
                        endif							
					enddo
				elseif(fileLine(i1-8:i1) .eq. 'sections:') then
					secNum = 0
					layerNum = 0
					secLayupRange(0) = 0
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
						i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-4:i1) .eq. 'type:') then
                                secNum = secNum + 1
								read(fileLine(i1+1:i1+64),*) sectionType(secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
						    elseif(fileLine(i1-8:i1) .eq. 'material:') then
							    read(fileLine(i1+1:i1+64),*) sectionMatName(secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-11:i1) .eq. 'orientation:') then
							    i2 = index(fileLine,'[')
								i3 = index(fileLine,']')
								read(fileLine(i2+1:i3-1),*) readReal(1:6)
								i4 = secNum*3 - 2
								sectionOrient(:,i4) = readReal(1:3)
								sectionOrient(:,i4+1) = readReal(4:6)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-7:i1) .eq. 'zOffset:') then
							    read(fileLine(i1+1:256),*) sectionZOffset(secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-6:i1) .eq. 'layers:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
										layerNum = layerNum + 1
										read(fileLine(i2+1:i3-1),*) layupMatName(layerNum), readReal(1:2)
										layupThickness(layerNum) = readReal(1)
										layupAngle(layerNum) = readReal(2)
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
								    i1 = index(fileLine,':')
								enddo
								secLayupRange(secNum) = layerNum
							elseif(fileLine(i1-4:i1) .eq. 'area:') then
							    read(fileLine(i1+1:256),*) beamProperties(1,secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-1:i1) .eq. 'I:') then
							    i2 = index(fileLine,'[')
								i3 = index(fileLine,']')
								read(fileLine(i2+1:i3-1),*) beamProperties(2:6,secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-1:i1) .eq. 'J:') then
							    read(fileLine(i1+1:256),*) beamProperties(7,secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-9:i1) .eq. 'stiffness:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							    i1 = index(fileLine,':')
								readMat(:,:) = r_0
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
									    read(fileLine(i2+1:i3-1),*) readInt(1:2), readReal(1)
										readMat(readInt(1),readInt(2)) = readReal(1)
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
							        i1 = index(fileLine,':')
								enddo
								do i4 = 2, 6
								    do i5 = 1, i4-1
									    if(readMat(i4,i5) .ne. r_0) then
										    readMat(i5,i4) = readMat(i4,i5)
										endif
									enddo
								enddo
								i3 = 1
								do i4 = 1, 6
								    do i5 = i4, 6
										beamStiffness(i3,secNum) = readMat(i4,i5)
										i3 = i3 + 1
									enddo
								enddo
							elseif(fileLine(i1-11:i1) .eq. 'expLoadCoef:') then
							    i2 = index(fileLine,'[')
								i3 = index(fileLine,']')
								read(fileLine(i2+1:i3-1),*) beamExpLoadCoef(:,secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-4:i1) .eq. 'mass:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							    i1 = index(fileLine,':')
								readMat(:,:) = r_0
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
									    read(fileLine(i2+1:i3-1),*) readInt(1:2), readReal(1)
										readMat(readInt(1),readInt(2)) = readReal(1)
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
							        i1 = index(fileLine,':')
								enddo
								do i4 = 2, 6
								    do i5 = 1, i4-1
									    if(readMat(i4,i5) .ne. r_0) then
										    readMat(i5,i4) = readMat(i4,i5)
										endif
									enddo
								enddo
								i3 = 1
								do i4 = 1, 6
								    do i5 = i4, 6
										beamMass(i3,secNum) = readMat(i4,i5)
										i3 = i3 + 1
									enddo
								enddo
							elseif(fileLine(i1-12:i1) .eq. 'conductivity:') then
							    i2 = index(fileLine,'[')
								i3 = index(fileLine,']')
								read(fileLine(i2+1:i3-1),*) beamThermCond(:,secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-8:i1) .eq. 'specHeat:') then
							    read(fileLine(i1+1:i1+64),*) beamSpecHeat(secNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-10:i1) .eq. 'elementSet:') then
                                do i2 = 1, numElSets
								    read(fileLine(i1+1:i1+64),*) readChar
									if(readChar .eq. elSetName(i2)) then
									    do i3 = elSetRange(i2-1) + 1, elSetRange(i2)
										    i4 = elementSets(i3)
											elementSection(i4) = secNum
										enddo
									endif
								enddo
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							else
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				elseif(fileLine(i1-9:i1) .eq. 'materials:') then
					matNum = 0
					GSet = 0
					nuSet = 0
					inMaxStress = 0
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
						i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-4:i1) .eq. 'name:') then
							    matNum = matNum + 1
								GSet = 0
					            nuSet = 0
					            inMaxStress = 0
								read(fileLine(i1+1:i1+64),*) materialName(matNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-7:i1) .eq. 'density:') then
							    read(fileLine(i1+1:i1+64),*) materialDensity(matNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-1:i1) .eq. 'E:') then
							    i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									read(fileLine(i2+1:i3-1),*) materialElastic(1:3,matNum)
								else
								    read(fileLine(i1+1:i1+64),*) materialElastic(1,matNum)
									materialElastic(2:3,matNum) = materialElastic(1,matNum)
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-2:i1) .eq. 'nu:') then
							    i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									read(fileLine(i2+1:i3-1),*) materialElastic(4:6,matNum)
								else
								    read(fileLine(i1+1:i1+64),*) materialElastic(4,matNum)
									materialElastic(5:6,matNum) = materialElastic(4,matNum)
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
								nuSet = 1
							elseif(fileLine(i1-1:i1) .eq. 'G:') then
							    i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									read(fileLine(i2+1:i3-1),*) materialElastic(7:9,matNum)
								else
								    read(fileLine(i1+1:i1+64),*) materialElastic(7,matNum)
									materialElastic(8:9,matNum) = materialElastic(7,matNum)
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
								GSet = 1
							elseif(fileLine(i1-9:i1) .eq. 'stiffness:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							    i1 = index(fileLine,':')
								readMat(:,:) = r_0
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
									    read(fileLine(i2+1:i3-1),*) readInt(1:2), readReal(1)
										readMat(readInt(1),readInt(2)) = readReal(1)
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
							        i1 = index(fileLine,':')
								enddo
								do i4 = 2, 6
								    do i5 = 1, i4-1
									    if(readMat(i4,i5) .ne. r_0) then
										    readMat(i5,i4) = readMat(i4,i5)
										endif
									enddo
								enddo
								i3 = 1
								do i4 = 1, 6
								    do i5 = i4, 6
										materialStiffMat(i3,matNum) = readMat(i4,i5)
										i3 = i3 + 1
									enddo
								enddo
							elseif(fileLine(i1-12:i1) .eq. 'conductivity:') then
							    i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									read(fileLine(i2+1:i3-1),*) materialThermCond(1:6,matNum)
								else
								    read(fileLine(i1+1:i1+64),*) materialThermCond(1,matNum)
									materialThermCond(2:3,matNum) = materialThermCond(1,matNum)
									materialThermCond(4:6,matNum) = r_0
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-9:i1) .eq. 'expansion:') then
							    i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									read(fileLine(i2+1:i3-1),*) materialThermExp(1:6,matNum)
								else
								    read(fileLine(i1+1:i1+64),*) materialThermExp(1,matNum)
									materialThermExp(2:3,matNum) = materialThermExp(1,matNum)
									materialThermExp(4:6,matNum) = r_0
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-12:i1) .eq. 'specificHeat:') then
							    read(fileLine(i1+1:i1+64),*) materialSpecHeat(matNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-9:i1) .eq. 'maxStress:') then
							    inMaxStress = 1
							elseif(fileLine(i1-9:i1) .eq. 'maxStrain:') then
							    inMaxStress = 0
							elseif(fileLine(i1-7:i1) .eq. 'tensile:') then
							    i2 = index(fileLine,'[')
								i3 = index(fileLine,']')
								if(inMaxStress .eq. 1) then
							        read(fileLine(i2+1:i3-1),*) materialMaxStress(1:3,matNum)
								else
								    read(fileLine(i2+1:i3-1),*) materialMaxStrain(1:3,matNum)
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-11:i1) .eq. 'compressive:') then
							    i2 = index(fileLine,'[')
								i3 = index(fileLine,']')
							    if(inMaxStress .eq. 1) then
							        read(fileLine(i2+1:i3-1),*) materialMaxStress(4:6,matNum)
								else
								    read(fileLine(i2+1:i3-1),*) materialMaxStrain(4:6,matNum)
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-5:i1) .eq. 'shear:') then
							    i2 = index(fileLine,'[')
								i3 = index(fileLine,']')
							    if(inMaxStress .eq. 1) then
							        read(fileLine(i2+1:i3-1),*) materialMaxStress(7:9,matNum)
								else
								    read(fileLine(i2+1:i3-1),*) materialMaxStrain(7:9,matNum)
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-15:i1) .eq. 'maxStrainEnergy:') then
							    read(fileLine(i1+1:i1+64),*) materialMaxStrnEngy(matNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							else
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
					if(GSet .eq. 0) then
					    materialElastic(7,matNum) = r_p5*materialElastic(1,matNum)/(r_1 + materialElastic(4,matNum))
						materialElastic(8:9,matNum) = materialElastic(7,matNum)
					elseif(nuSet .eq. 0) then
					    materialElastic(4,matNum) = r_p5*materialElastic(1,matNum)/materialElastic(7,matNum) + r_1
						materialElastic(5:6,matNum) = materialElastic(4,matNum)
					endif
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
		enddo
		
		close(1)
		
		!! Populate the section material id numbers using the names read in
		do i1 = 1, numSec
		    do i2 = 1, numMats
			    if(materialName(i2) .eq. sectionMatName(i1)) then
				    sectionMatId(i1) = i2
				endif
			enddo
		enddo
		
		do i1 = 1, layupSize
		    do i2 = 1, numMats
			    if(materialName(i2) .eq. layupMatName(i1)) then
				    layupMatId(i1) = i2
				endif
			enddo
		enddo

		!! Set the direction cosine matrices based on input
		do i1 = 1, numSec
		    i2 = i1*3
			v1 = sectionOrient(:,i2-2)
			v2 = sectionOrient(:,i2-1)
			mag = sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3))
			if(mag .lt. 1e-12) then
			    sectionOrient(:,i2-2) = (/r_1,r_0,r_0/)
				sectionOrient(:,i2-1) = (/r_0,r_1,r_0/)
				sectionOrient(:,i2) = (/r_0,r_0,r_1/)
			else
			    v1 = (r_1/mag)*v1
				v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
				v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
				v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
				mag = sqrt(v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3))
				v3 = (r_1/mag)*v3
				v2(1) = v3(2)*v1(3) - v3(3)*v1(2)
				v2(2) = v3(3)*v1(1) - v3(1)*v1(3)
				v2(3) = v3(1)*v1(2) - v3(2)*v1(1)
				sectionOrient(1,i2-2:i2) = v1
				sectionOrient(2,i2-2:i2) = v2
				sectionOrient(3,i2-2:i2) = v3
			endif
		enddo
	    
	end subroutine readModelInput
	
	subroutine readLoads(inFileName)
	    implicit none
		
		character(len=128), intent(in) :: inFileName
		
		character(len=256) :: fileLine
		character(len=16) :: subStrings(16)
		integer :: iosVal
		integer :: i1, i2, i3, i4, i5, i6, found, inElSet, inMaxStress
		integer :: elType, numNds, secNum, matNum, GSet, nuSet
		integer :: layerNum, termNum, mpcNum, loadNum, eSNum, nSNum, elSi, ndSi
		real*8 :: mag
		real*8 :: v1(3), v2(3), v3(3)
        integer :: readInt(9)
        real*8 :: readReal(21)
        character(len=64) :: readChar		
		
		open(unit=1, file=inFileName, blank='NULL', action='read')
		
		subStrings(1) = 'nodes:'
		subStrings(2) = 'elements:'
		subStrings(3) = 'sets:'
		subStrings(4) = 'sections:'
		subStrings(5) = 'materials:'
		subStrings(6) = 'constraints:'
		subStrings(7) = 'loads:'
		subStrings(8) = 'initialState:'
		
		fileLine(1:15) = '               '
		numLds = 0
		sizeLds = 0
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)  !! IOSTAT = 0 perfect, > 0 problem occured, < 0 end of file
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
		    if(i1 .gt. 0) then
				if(fileLine(i1-5:i1) .eq. 'loads:') then
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
					    i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-4:i1) .eq. 'type:') then
							    numLds = numLds + 1
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-8:i1) .eq. 'setLoads:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'-')
									if(i2 .gt. 0) then
									    sizeLds = sizeLds + 1
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
								    i1 = index(fileLine,':')
								enddo
							else
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
		enddo
		
		if(numLds .gt. 0 .and. .not. allocated(loadNodes)) then
		    allocate(loadNodes(sizeLds))
			allocate(inputLoads(7,sizeLds))
			allocate(loadsRange(0:numLds))
			allocate(loadsActTime(2,numLds))
		    allocate(loadType(numLds))
			
			inputLoads(:,:) = r_0
			loadsRange(:) = r_0
			loadsActTime(:,:) = r_0
		endif
		
		rewind(1)
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)  !! IOSTAT = 0 perfect, > 0 problem occured, < 0 end of file
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
		    if(i1 .gt. 0) then
				if(fileLine(i1-5:i1) .eq. 'loads:') then
					loadsRange(0) = 0
					termNum = 0
					loadNum = 0
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
					    i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-4:i1) .eq. 'type:') then
							    loadNum = loadNum + 1
							    read(fileLine(i1+1:i1+64),*) loadType(loadNum)
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
						    elseif(fileLine(i1-10:i1) .eq. 'activeTime:') then
								i2 = index(fileLine,'[')
								if(i2 .gt. 0) then
								    i3 = index(fileLine,']')
									read(fileLine(i2+1:i3-1),*) loadsActTime(1:2,loadNum)
								else
								    read(fileLine(i1+1:i1+64),*) loadsActTime(1,loadNum)
									loadsActTime(2,loadNum) = r_1*1e+100_8
								endif
								read(1,'(A)',iostat=iosVal) fileLine(16:256)
							elseif(fileLine(i1-8:i1) .eq. 'setLoads:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    termNum = termNum + 1
									    i3 = index(fileLine,']')
										if(loadType(loadNum) .eq. 'gravitational') then
										    read(fileLine(i2+1:i3-1),*) loadNodes(termNum), inputLoads(1:3,termNum)
										elseif(loadType(loadNum) .eq. 'surfacePressure') then
										    read(fileLine(i2+1:i3-1),*) loadNodes(termNum), inputLoads(1:5,termNum)
										elseif(loadType(loadNum) .eq. 'bodyHeatGen') then
										    read(fileLine(i2+1:i3-1),*) loadNodes(termNum), inputLoads(1,termNum)
										elseif(loadType(loadNum) .eq. 'surfaceFlux') then
										    read(fileLine(i2+1:i3-1),*) loadNodes(termNum), inputLoads(1:5,termNum)
										else
										    read(fileLine(i2+1:i3-1),*) loadNodes(termNum), inputLoads(1:7,termNum)
										endif
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
								    i1 = index(fileLine,':')
								enddo
								loadsRange(loadNum) = termNum
							else
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
		enddo
		
		close(1)
	  
	end subroutine readLoads
	
	subroutine readConstraints(inFileName)
	    implicit none
		
		character(len=128), intent(in) :: inFileName
		
		character(len=256) :: fileLine
		character(len=16) :: subStrings(16)
		integer :: iosVal
		integer :: i1, i2, i3, i4, i5, i6, found, inElSet, inMaxStress
		integer :: elType, numNds, secNum, matNum, GSet, nuSet
		integer :: layerNum, termNum, mpcNum, loadNum, eSNum, nSNum, elSi, ndSi
		real*8 :: mag
		real*8 :: v1(3), v2(3), v3(3)
        integer :: readInt(9)
        real*8 :: readReal(21)
        character(len=64) :: readChar		
		
		open(unit=1, file=inFileName, blank='NULL', action='read')
		
		subStrings(1) = 'nodes:'
		subStrings(2) = 'elements:'
		subStrings(3) = 'sets:'
		subStrings(4) = 'sections:'
		subStrings(5) = 'materials:'
		subStrings(6) = 'constraints:'
		subStrings(7) = 'loads:'
		subStrings(8) = 'initialState:'
		
		fileLine(1:15) = '               '
		numMPC = 0
		mpcSize = 0
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)  !! IOSTAT = 0 perfect, > 0 problem occured, < 0 end of file
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
		    if(i1 .gt. 0) then
				if(fileLine(i1-11:i1) .eq. 'constraints:') then
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
						i1 = index(fileLine,':')
						i2 = index(fileLine,'[')
						if(i1 .gt. 0) then
						    if(fileLine(i1-5:i1) .eq. 'terms:') then
								numMPC = numMPC + 1
							endif
						elseif(i2 .gt. 0) then
						    mpcSize = mpcSize + 1
						endif
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
		enddo
		
		if(numMPC .gt. 0 .and. .not. allocated(mpcEqn)) then
		    allocate(mpcEqn(mpcSize))
			allocate(mpcNode(mpcSize))
			allocate(mpcDof(mpcSize))
			allocate(mpcCoef(mpcSize))
			allocate(mpcRHS(numMPC))
			
			mpcEqn(:) = r_0
			mpcDof(:) = r_0
			mpcCoef(:) = r_0
			mpcRHS(:) = r_0
		endif
		
		rewind(1)
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)  !! IOSTAT = 0 perfect, > 0 problem occured, < 0 end of file
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
		    if(i1 .gt. 0) then
				if(fileLine(i1-11:i1) .eq. 'constraints:') then
					termNum = 0
					mpcNum = 0
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
						i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-5:i1) .eq. 'terms:') then
							    mpcNum = mpcNum + 1
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    termNum = termNum + 1
									    i3 = index(fileLine,']')
										read(fileLine(i2+1:i3-1),*) mpcNode(termNum), mpcDof(termNum), mpcCoef(termNum)
										mpcEqn(termNum) = mpcNum
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
									i1 = index(fileLine,':')
								enddo
							elseif(fileLine(i1-3:i1) .eq. 'rhs:') then
							    read(fileLine(i1+1:i1+64),*) mpcRHS(mpcNum)
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
								i1 = index(fileLine,':')
							else
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
		enddo
		
		close(1)
	    
	end subroutine readConstraints
	
	subroutine readInitialState(inFileName)
	    implicit none
		
		character(len=128), intent(in) :: inFileName
		
		character(len=256) :: fileLine
		character(len=16) :: subStrings(16)
		integer :: iosVal
		integer :: i1, i2, i3, i4, i5, i6, found, inElSet, inMaxStress
		integer :: elType, numNds, secNum, matNum, GSet, nuSet
		integer :: layerNum, termNum, mpcNum, loadNum, eSNum, nSNum, elSi, ndSi
		real*8 :: mag
		real*8 :: v1(3), v2(3), v3(3)
        integer :: readInt(9)
        real*8 :: readReal(21)
        character(len=64) :: readChar		
		
		open(unit=1, file=inFileName, blank='NULL', action='read')
		
		subStrings(1) = 'nodes:'
		subStrings(2) = 'elements:'
		subStrings(3) = 'sets:'
		subStrings(4) = 'sections:'
		subStrings(5) = 'materials:'
		subStrings(6) = 'constraints:'
		subStrings(7) = 'loads:'
		subStrings(8) = 'initialState:'
		
		fileLine(1:15) = '               '
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)  !! IOSTAT = 0 perfect, > 0 problem occured, < 0 end of file
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
		    if(i1 .gt. 0) then
				if(fileLine(i1-12:i1) .eq. 'initialState:') then
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
					call findSubstrings(found,fileLine,subStrings,8)
					do while(found .eq. 0 .and. iosVal .eq. 0)
					    i1 = index(fileLine,':')
						if(i1 .gt. 0) then
						    if(fileLine(i1-12:i1) .eq. 'displacement:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					            i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
										read(fileLine(i2+1:i3-1),*) readChar, readReal(1:6)
										read(readChar,*,err=1459) readInt(1)
										initialDisp(:,readInt(1)) = readReal(1:6)
										goto 1467
1459							        do i4 = 1, numNdSets
                                            if(ndSetName(i4) .eq. readChar) then
											    do i5 = ndSetRange(i4-1) + 1, ndSetRange(i4)
												    i6 = nodeSets(i5)
													initialDisp(:,i6) = readReal(1:6)
												enddo
											endif
                                        enddo
1467									i1 = index(fileLine,':')
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
					                i1 = index(fileLine,':')
								enddo
							elseif(fileLine(i1-8:i1) .eq. 'velocity:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					            i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
										read(fileLine(i2+1:i3-1),*) readChar, readReal(1:6)
										read(readChar,*,err=1483) readInt(1)
										initialVel(:,readInt(1)) = readReal(1:6)
										goto 1491
1483							        do i4 = 1, numNdSets
                                            if(ndSetName(i4) .eq. readChar) then
											    do i5 = ndSetRange(i4-1) + 1, ndSetRange(i4)
												    i6 = nodeSets(i5)
													initialVel(:,i6) = readReal(1:6)
												enddo
											endif
                                        enddo
1491									i1 = index(fileLine,':')
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
					                i1 = index(fileLine,':')
								enddo
							elseif(fileLine(i1-12:i1) .eq. 'acceleration:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					            i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
										read(fileLine(i2+1:i3-1),*) readChar, readReal(1:6)
										read(readChar,*,err=1507) readInt(1)
										initialAcc(:,readInt(1)) = readReal(1:6)
										goto 1515
1507							        do i4 = 1, numNdSets
                                            if(ndSetName(i4) .eq. readChar) then
											    do i5 = ndSetRange(i4-1) + 1, ndSetRange(i4)
												    i6 = nodeSets(i5)
													initialAcc(:,i6) = readReal(1:6)
												enddo
											endif
                                        enddo
1515									i1 = index(fileLine,':')
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
					                i1 = index(fileLine,':')
								enddo
						    elseif(fileLine(i1-11:i1) .eq. 'temperature:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					            i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
										read(fileLine(i2+1:i3-1),*) readChar, readReal(1)
										read(readChar,*,err=1531) readInt(1)
										initialTemp(readInt(1)) = readReal(1)
										goto 1539
1531							        do i4 = 1, numNdSets
                                            if(ndSetName(i4) .eq. readChar) then
											    do i5 = ndSetRange(i4-1) + 1, ndSetRange(i4)
												    i6 = nodeSets(i5)
													initialTemp(i6) = readReal(1)
												enddo
											endif
                                        enddo
1539									i1 = index(fileLine,':')
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
					                i1 = index(fileLine,':')
								enddo
							elseif(fileLine(i1-14:i1) .eq. 'tempChangeRate:') then
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					            i1 = index(fileLine,':')
								do while(i1 .eq. 0 .and. iosVal .eq. 0)
								    i2 = index(fileLine,'[')
									if(i2 .gt. 0) then
									    i3 = index(fileLine,']')
										read(fileLine(i2+1:i3-1),*) readChar, readReal(1:6)
										read(readChar,*,err=1555) readInt(1)
										initialTdot(readInt(1)) = readReal(1)
										goto 1563
1555							        do i4 = 1, numNdSets
                                            if(ndSetName(i4) .eq. readChar) then
											    do i5 = ndSetRange(i4-1) + 1, ndSetRange(i4)
												    i6 = nodeSets(i5)
													initialTdot(i6) = readReal(1)
												enddo
											endif
                                        enddo
1563									i1 = index(fileLine,':')
									endif
									read(1,'(A)',iostat=iosVal) fileLine(16:256)
					                i1 = index(fileLine,':')
								enddo
							else
							    read(1,'(A)',iostat=iosVal) fileLine(16:256)
							endif
						else
						    read(1,'(A)',iostat=iosVal) fileLine(16:256)
						endif
						call findSubstrings(found,fileLine,subStrings,8)
					enddo
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
		enddo
		
		close(1)
	    
	end subroutine readInitialState
	
	subroutine readDesignVarInput(inFileName)
	    implicit none
		
		character(len=128), intent(in) :: inFileName
		
		character(len=256) :: fileLine
		character(len=64) :: readChar
		real*8 :: readReal(2)
		integer :: iosVal
		integer :: elListLen, ndListLen, readInt(2)
		integer :: dVarNum, coefNum
		integer :: i1, i2, i3, i4, inserted
		
		open(unit=1, file=inFileName, blank='NULL', action='read')
		
		numDVar = 0
		coefListSize = 0
		
		fileLine(1:15) = '               '
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256)
		
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
			if(i1 .gt. 0) then
			    if(fileLine(i1-8:i1) .eq. 'category:') then
				    numDVar = numDVar + 1
					read(1,'(A)',iostat=iosVal) fileLine(16:256) 
					i1 = index(fileLine,':')			
				elseif(fileLine(i1-12:i1) .eq. 'coefficients:') then
				    read(fileLine(i1+1:i1+64),*,err=1339,end=1339) readReal(1)
					coefListSize = coefListSize + 1
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				    i1 = index(fileLine,':')
					goto 1342
1339                read(1,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)
					     i2 = index(fileLine,'-')
						 if(i2 .gt. 0) then
						     coefListSize = coefListSize + 1
						 endif
						 read(1,'(A)',iostat=iosVal) fileLine(16:256) 
					     i1 = index(fileLine,':')
					enddo
1342                i1 = i1
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,':')
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,':')
			endif
		enddo
		
		
		allocate(r_dVec(numDVar))
		allocate(c_dVec(numDVar))
		allocate(dLdD(numDVar))
		allocate(dCategory(numDVar))
		allocate(dComponent(numDVar))
		allocate(dLayer(numDVar))
		allocate(dActTime(2,numDVar))
		allocate(dElSet(numDVar))
		allocate(dNdSet(numDVar))
		allocate(dCoefList(coefListSize))
		allocate(dCoefListRange(0:numDVar))
		allocate(dsumTermdD(numDVar))
		allocate(dtotVoldD(numDVar))
		
		r_dVec(:) = r_0
		c_dVec(:) = c_0
		dLdD(:) = r_0
		dComponent(:) = r_0
		dLayer(:) = r_0
		dActTime(:,:) = r_0
		dsumTermdD(:) = r_0
		dtotVoldD(:) = r_0
		dCoefList(:) = r_0
		dCoefListRange(:) = 0
		
		rewind(1)
		
		dVarNum = 0
		coefNum = 0
		
		read(1,'(A)',iostat=iosVal) fileLine(16:256) 
		
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
			if(i1 .gt. 0) then
			    if(fileLine(i1-8:i1) .eq. 'category:') then
				    dVarNum = dVarNum + 1
					elListLen = 0
					ndListLen = 0
					read(fileLine(i1+1:i1+64),*) dCategory(dVarNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256) 
					i1 = index(fileLine,':')
				elseif(fileLine(i1-9:i1) .eq. 'component:') then
				    read(fileLine(i1+1:i1+64),*) dComponent(dVarNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,':')
				elseif(fileLine(i1-5:i1) .eq. 'layer:') then
				    read(fileLine(i1+1:i1+64),*) dLayer(dVarNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,':')
				elseif(fileLine(i1-10:i1) .eq. 'activeTime:') then
				    i2 = index(fileLine,'[')
					if(i2 .gt. 0) then
					    i3 = index(fileLine,']')
						read(fileLine(i2+1:i3-1),*) dActTime(1:2,dVarNum)
					else
					    read(fileLine(i1+1:i1+64),*) dActTime(1,dVarNum)
						dActTime(2,dVarNum) = r_1*1e+100_8
					endif
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,':')
				elseif(fileLine(i1-10:i1) .eq. 'elementSet:') then
                    read(fileLine(i1+1:i1+64),*) dElSet(dVarNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,':')
				elseif(fileLine(i1-7:i1) .eq. 'nodeSet:') then
                    read(fileLine(i1+1:i1+64),*) dNdSet(dVarNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,':')
				elseif(fileLine(i1-12:i1) .eq. 'coefficients:') then
                    read(fileLine(i1+1:i1+64),*,err=1440,end=1440) readReal(1)
					coefNum = coefNum + 1
					dCoefList(coefNum) = readReal(1)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				    i1 = index(fileLine,':')
					goto 1450
1440                read(1,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)
					     i2 = index(fileLine,'-')
						 if(i2 .gt. 0) then
						     coefNum = coefNum + 1
							 read(fileLine(i2+1:i2+64),*) dCoefList(coefNum)
						 endif
						 read(1,'(A)',iostat=iosVal) fileLine(16:256) 
					     i1 = index(fileLine,':')
					enddo
1450                dCoefListRange(dVarNum) = coefNum
                else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,':')
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,':')
			endif
		enddo
		
		close(1)
		
		do i1 = 1, numDVar
		    if(dCoefListRange(i1) .eq. 0) then
			    dCoefListRange(i1) = dCoefListRange(i1-1)
			endif
		enddo
		
	end subroutine readDesignVarInput
	
	subroutine readObjectiveInput(inFileName)
	    implicit none
		
		character(len=128), intent(in) :: inFileName
		
		character(len=256) :: fileLine
		real*8 :: readReal
		integer :: iosVal, termNum, objTgti
		integer :: i1, i2, i3
		
		fileLine(1:15) = '               '
		open(unit=1, file=inFileName, blank='NULL', action='read')
		
		numObjTerms = 0
		objTgtSize = 0
		read(1,'(A)',iostat=iosVal) fileLine(16:256)
        do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
			if(i1 .gt. 0) then
			    if(fileLine(i1-8:i1) .eq. 'category:') then
				    numObjTerms = numObjTerms + 1
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-11:i1) .eq. 'targetValue:') then
				    i2 = index(fileLine(i1+1:i1+64),'max')
					if(i2 .gt. 1) then
					    objTgtSize = objTgtSize + 1
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
					else
					    read(fileLine(i1+1:i1+64),*,err=1436,end=1436) readReal
						objTgtSize = objTgtSize + 1
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
						goto 1446
1436					read(1,'(A)',iostat=iosVal) fileLine(16:256)
                        i1 = index(fileLine,':')
						do while(i1 .eq. 0 .and. iosVal .eq. 0)
						    i2 = index(fileLine,'-')
							if(i2 .gt. 0) then
							    objTgtSize = objTgtSize + 1
							endif
							read(1,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,':')
						enddo
1446					i1 = index(fileLine,':')
					endif
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
        enddo
		
		if(numObjTerms .gt. 0) then
		    allocate(objVal(numObjTerms))
		    allocate(objCategory(numObjTerms))
			allocate(objOperator(numObjTerms))
			allocate(objActTime(2,numObjTerms))
			allocate(objComponent(numObjTerms))
			allocate(objLayer(numObjTerms))
			allocate(objCoef(numObjTerms))
			allocate(objExp(numObjTerms))
			allocate(objElSet(numObjTerms))
			allocate(objElSetId(numObjTerms))
			allocate(objNdSet(numObjTerms))
			allocate(objNdSetId(numObjTerms))
			allocate(objTgtTag(numObjTerms))
			allocate(objTgtVal(objTgtSize))
			allocate(objTgtRange(0:numObjTerms))
			
			objVal(:) = r_0
			objActTime(:,:) = r_0
			objComponent(:) = r_0
			objLayer(:) = r_0
			objCoef(:) = r_0
			objExp(:) = r_0
			objElSetId(:) = r_0
			objNdSetId(:) = r_0
			objTgtVal(:) = r_0
			objTgtRange(:) = r_0
		endif
		
		rewind(1)
		
		do i1 = 1, numObjTerms
		    objTgtTag(i1) = 'none'
		enddo
		
		termNum = 0
		objTgti = 0
		read(1,'(A)',iostat=iosVal) fileLine(16:256)
        do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
			if(i1 .gt. 0) then
			    if(fileLine(i1-8:i1) .eq. 'category:') then
				    termNum = termNum + 1
				    read(fileLine(i1+1:i1+64),*) objCategory(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-8:i1) .eq. 'operator:') then
				    read(fileLine(i1+1:i1+64),*) objOperator(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-10:i1) .eq. 'activeTime:') then
				    i2 = index(fileLine,'[')
					if(i2 .gt. 0) then
					    i3 = index(fileLine,']')
					    read(fileLine(i2+1:i3-1),*) objActTime(1:2,termNum)
					else
					    read(fileLine(i1+1:i1+64),*) objActTime(1,termNum)
						objActTime(2,termNum) = r_1*1e+100_8
					endif
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-9:i1) .eq. 'component:') then
				    read(fileLine(i1+1:i1+64),*) objComponent(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-5:i1) .eq. 'layer:') then
				    read(fileLine(i1+1:i1+64),*) objLayer(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-11:i1) .eq. 'coefficient:') then
				    read(fileLine(i1+1:i1+64),*) objCoef(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-8:i1) .eq. 'exponent:') then
				    read(fileLine(i1+1:i1+64),*) objExp(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-10:i1) .eq. 'elementSet:') then
				    read(fileLine(i1+1:i1+64),*) objElSet(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-7:i1) .eq. 'nodeSet:') then
				    read(fileLine(i1+1:i1+64),*) objNdSet(termNum)
					read(1,'(A)',iostat=iosVal) fileLine(16:256)
				elseif(fileLine(i1-11:i1) .eq. 'targetValue:') then
				    i2 = index(fileLine(i1+1:i1+64),'max')
					if(i2 .gt. 1) then
					    objTgti = objTgti + 1
						read(fileLine(i1+1:i1+64),*) objTgtTag(termNum)
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
					else
					    read(fileLine(i1+1:i1+64),*,err=1530,end=1530) readReal
						objTgti = objTgti + 1
						objTgtVal(objTgti) = readReal
						objTgtTag(termNum) = 'single'
						read(1,'(A)',iostat=iosVal) fileLine(16:256)
						goto 1540
1530					read(1,'(A)',iostat=iosVal) fileLine(16:256)
                        i1 = index(fileLine,':')
						do while(i1 .eq. 0 .and. iosVal .eq. 0)
						    i2 = index(fileLine,'-')
							if(i2 .gt. 0) then
							    objTgti = objTgti + 1
								read(fileLine(i2+1:i2+64),*) objTgtVal(objTgti)
							endif
							read(1,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,':')
						enddo
						objTgtTag(termNum) = 'vector'
1540					i1 = index(fileLine,':')
					endif
					objTgtRange(termNum) = objTgti
				else
				    read(1,'(A)',iostat=iosVal) fileLine(16:256)
				endif
			else
			    read(1,'(A)',iostat=iosVal) fileLine(16:256)
			endif
        enddo		
		
		close(1)
		
	end subroutine readObjectiveInput
	
	subroutine readDesignVarValues(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		
		character(len=128) :: fileLine
		real*8 :: readReal
		integer :: i1, i2, i3, i4, readInt, iosVal

        r_dVec(:) = r_0
        fileLine(1:10) = '          '		
		open(unit=1586, file=fileName, blank='NULL', action='read')
		
		read(1586,'(A)',iostat=iosVal) fileLine(11:128)
		do while(iosVal .eq. 0)
		    i2 = index(fileLine,'-')
			if(i2 .gt. 0) then
			    i3 = index(fileLine,'[')
				if(i3 .gt. 0) then
				    i4 = index(fileLine,']')
					read(fileLine(i3+1:i4-1),*) readInt, readReal
					if(readInt .le. numDVar) then
					    r_dVec(readInt) = readReal
					else
					    write(lfUnit,*) 'Error: Design variable label ', readInt, ' in value input file'
						write(lfUnit,*) 'exceeds the size of the design variable vector.'
					endif
				endif
			endif
			read(1586,'(A)',iostat=iosVal) fileLine(11:128)
		enddo
		
		close(1586)
	
	end subroutine readDesignVarValues
	
	subroutine readNodeResults(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		
		character(len=128) :: fileLine
		real*8 :: readReal(6)
		integer :: i1, i2, i3, i4, i5, i6, i7, readInt, iosVal
		
		fileLine(1:15) = '               '
		open(unit=1622, file=fileName, blank='NULL', action='read')
		
		read(1622,'(A)',iostat=iosVal) fileLine(16:128)
		do while(iosVal .eq. 0)
		    i1 = index(fileLine,':')
			if(i1 .gt. 0) then
			    if(fileLine(i1-11:i1) .eq. 'temperature:') then
				    read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)
					    i2 = index(fileLine,'-')
						if(i2 .gt. 0) then
						    i3 = index(fileLine,'[')
							i4 = index(fileLine,']')
							read(fileLine(i3+1:i4-1),*) readInt, readReal(1)
							if(readInt .le. numNodes) then
							    i5 = currentRank(readInt)
								nodeTemp(i5) = readReal(1)
							endif
						endif
						read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					    i1 = index(fileLine,':')
					enddo
				elseif(fileLine(i1-4:i1) .eq. 'tdot:') then
				    read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)
					    i2 = index(fileLine,'-')
						if(i2 .gt. 0) then
						    i3 = index(fileLine,'[')
							i4 = index(fileLine,']')
							read(fileLine(i3+1:i4-1),*) readInt, readReal(1)
							if(readInt .le. numNodes) then
							    i5 = currentRank(readInt)
								nodeTdot(i5) = readReal(1)
							endif
						endif
						read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					    i1 = index(fileLine,':')
					enddo
				elseif(fileLine(i1-12:i1) .eq. 'displacement:') then
				    read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)
					    i2 = index(fileLine,'-')
						if(i2 .gt. 0) then
						    i3 = index(fileLine,'[')
							i4 = index(fileLine,']')
							read(fileLine(i3+1:i4-1),*) readInt, readReal(1:6)
							if(readInt .le. numNodes) then
							    do i5 = 1, 6
								    i6 = nDofIndex(i5,readInt)
								    if(i6 .gt. 0) then
									    nodeDisp(i6) = readReal(i5)
									endif
								enddo
							endif
						endif
						read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					    i1 = index(fileLine,':')
					enddo
				elseif(fileLine(i1-8:i1) .eq. 'velocity:') then
				    read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)
					    i2 = index(fileLine,'-')
						if(i2 .gt. 0) then
						    i3 = index(fileLine,'[')
							i4 = index(fileLine,']')
							read(fileLine(i3+1:i4-1),*) readInt, readReal(1:6)
							if(readInt .le. numNodes) then
							    do i5 = 1, 6
								    i6 = nDofIndex(i5,readInt)
								    if(i6 .gt. 0) then
									    nodeVel(i6) = readReal(i5)
									endif
								enddo
							endif
						endif
						read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					    i1 = index(fileLine,':')
					enddo
				elseif(fileLine(i1-12:i1) .eq. 'acceleration:') then
				    read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					i1 = index(fileLine,':')
					do while(i1 .eq. 0 .and. iosVal .eq. 0)
					    i2 = index(fileLine,'-')
						if(i2 .gt. 0) then
						    i3 = index(fileLine,'[')
							i4 = index(fileLine,']')
							read(fileLine(i3+1:i4-1),*) readInt, readReal(1:6)
							if(readInt .le. numNodes) then
							    do i5 = 1, 6
								    i6 = nDofIndex(i5,readInt)
								    if(i6 .gt. 0) then
									    nodeAcc(i6) = readReal(i5)
									endif
								enddo
							endif
						endif
						read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					    i1 = index(fileLine,':')
					enddo
				else
				    read(1622,'(A)',iostat=iosVal) fileLine(16:128)
					i1 = index(fileLine,':')
				endif
			endif
		enddo
		
		close(1622)
		
	end subroutine readNodeResults
	
!! Binary
	
	subroutine readBinarySolution(stepNum,errFlag)
	    implicit none
		
		integer, intent(in) :: stepNum
		integer, intent(out) :: errFlag
		
		character(len=16) :: stepString
		character(len=32) :: fullFileName
	
        errFlag = 0	
		write(stepString,'(I0)') stepNum
        fullFileName = 'SolutionHistory' // stepString // '.out'
	    open(unit=40, file=fullFileName, form='unformatted', action='read', err=1603)
		
		read(40) prevTemp(:)
		read(40) prevTdot(:)
		read(40) prevDisp(:)
		if(intVecSize .gt. 1) then
		    read(40) prevIntDisp(:)
		endif
		read(40) prevVel(:)
		read(40) prevAcc(:)
		
		close(40)
		return

1603    errFlag = 1
		
	end subroutine readBinarySolution

end module ASenD_input