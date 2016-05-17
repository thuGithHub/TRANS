SUBROUTINE INPUT_Read_GlobalInfo
    USE Global
    IMPLICIT NONE
    CHARACTER (LEN=100), PARAMETER :: GlobalInfoFile ='dimension_global.dat'
    INTEGER:: itemp, i, j, idtemp

    OPEN(UNIT=11,FILE=trim(FilePathPrefix) // GlobalInfoFile, MODE='READ')  !, SHARE='DENYWR')
    READ(11,*)    
    READ(11,*) Num_Block_All !Max_Block
    NAllBlock=Num_Block_All
    allocate(PBlock(Num_Block_All)) 
    allocate( AllBlockName(Num_Block_All) )
    IF (IF_parallel) THEN
        ALLOCATE(ProcNo_All(Num_Block_All))
        ALLOCATE(BlockNo_All(Num_Block_All))
        READ(11,*)
        READ(11,*) itemp
        IF (itemp .ne. NumProcs) THEN
            WRITE(*,*) 'mpi core numbers wrong???'
            STOP
        ENDIF
        READ(11,*)
        DO i=1, NumProcs
            READ(11,*) itemp !the Max_block for that core
            IF (i .eq. MyID+1) THEN
                Max_Block = itemp
                NBlock=itemp
                ALLOCATE(M_blockids(Max_block))
            ENDIF
            DO j=1, itemp
                READ(11,*) idtemp
                ProcNo_All(idtemp) = i-1
                BlockNo_All(idtemp) = j
                PBlock(idtemp)=i
                IF (i .eq. MyID+1) THEN
                    M_blockids(j) = idtemp
                ENDIF
            ENDDO
        ENDDO
    ELSE
        ! only one, sequential
        Max_Block = Num_Block_All
        ALLOCATE(M_blockids(Max_block))
        DO i=1, Max_Block
            M_blockids(i)=i
        ENDDO
    ENDIF
    !added by ydd, read turbomachinery information
    read(11,*)   !P0in, T0in, ma_Income and Pouthub,and blade number
    read(11,*)P0in,T0in,Pouthub,nblade,Npassage
    read(11,*)  !omega,x,y,z, by rpm,should be transformed
    read(11,*)omega(1),omega(2),omega(3)
    read(11,*)      !Original point of rotating, default vaules:(0,0,0)
    read(11,*)Rot_Ori(1),Rot_Ori(2),Rot_Ori(3)
    
    CLOSE(11)   
    
    !Debug
    IF(IF_PARALLEL .and. Myid==0) THEN
        WRITE(*,*) 'Debug Global Info'
        DO i=1,Num_Block_All
            WRITE(*,*) '[', i, ',', ProcNo_All(i), ',', BlockNo_All(i), ']' 
        ENDDO
    ENDIF
    ALLOCATE(AllBlocks(Max_Block))
    
END SUBROUTINE INPUT_Read_GlobalInfo
