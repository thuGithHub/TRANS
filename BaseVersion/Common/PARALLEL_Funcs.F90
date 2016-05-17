SUBROUTINE PARALLEL_Initialize
    use Global
    implicit none

      if(If_parallel == 1)then

         call MPI_INIT( ierr )
         call MPI_COMM_RANK( MPI_COMM_WORLD,myid,ierr )
         call MPI_COMM_SIZE( MPI_COMM_WORLD,numprocs,ierr )
         write(*,*)  "Process ", myid, " of ", numprocs, " is alive!"   ! numprocs���ܼƺ�������myid����ǰ����ĺ��ı�š�
!     	     call CHANGEDIR(myid)
!         write (FilePathPrefix, '(A2,I3.3,A)') 'id', myid, '/'
        MyProc=myid+1
        FilePathPrefix = ''
      else
        FilePathPrefix = ''  !Is there anything wrong? (for length = 6)
      endif
    call MPI_Buffer_Attach(bufferForISend,bufferSize, ierr)

END SUBROUTINE PARALLEL_Initialize




SUBROUTINE PARALLEL_Dimensions
    use Global
    implicit none
    !INTEGER:: Idx_Temp
    INTEGER:: iBlock

      if (IF_parallel == 1) then 

          !ALLOCATE(M_blockids(Max_block))
          !DO iBlock = 1, Max_Block
          !  M_blockids(iBlock) = AllBlocks(iBlock)%Id_Present_Blk
          !ENDDO
          
          NmpiBCpairs_this = 0
          do iBlock = 1, Max_Block
            CALL PARALLEL_define_BC_getids_blk(iBlock)
          enddo

    !	  do idx_temp=1,NmpiBCpairs_this
    !	    write(*,*) myid, (Mpi_bcpairs_this(idx_temp, j),j=1,7)
    !	  enddo

    
          CALL PARALLEL_define_MPI_BCpairs   !BCpairs
      
    !	  CALL define_MPI_round  !TAG sort
      endif

END SUBROUTINE PARALLEL_Dimensions



SUBROUTINE PARALLEL_BCs(sweep)
    use Global
    implicit none
    INTEGER:: sweep  !BC sweeps
    integer:: iBlock, idx_TAGS, Mpitag_now, idx_pairs
    integer:: Mpi_bcnum1, Mpi_id1

         if (If_parallel == 1) then
           DO idx_TAGS = 1,NmpiTAGs
            Mpitag_now = Mpi_TAGS(idx_TAGS)
            DO idx_pairs = 1,Nmpibcpairs
            if ((Mpi_BCpairs(idx_pairs, 1) == Mpitag_now) .and.(myid == Mpi_BCPairs(idx_pairs, 2))) then
!            if ((myid == Mpi_BCPairs(idx_pairs, 2))) then
              DO iBlock = 1,Max_Block
                 if (Mpi_BCpairs(idx_pairs, 3) == M_blockids( iBlock )) then  !call
                    Mpi_bcnum1 = Mpi_BCPairs(idx_pairs, 4)
                    Mpi_id1 = Mpi_BCPairs(idx_pairs, 5)
                    call PARALLEL_BC_sendrecv_blk(iBlock, Mpi_bcnum1,Mpitag_now,Mpi_id1,    sweep)
                 endif
              ENDDO
            endif
            ENDDO

!	write(*,*) myid, idx_TAGS, NmpiTAGs
            CALL MPI_BARRIER( MPI_COMM_WORLD,ierr )
           ENDDO
!!!!!!
            CALL MPI_BARRIER( MPI_COMM_WORLD,ierr )
           endif

END SUBROUTINE PARALLEL_BCs



SUBROUTINE PARALLEL_Finish
    use Global
    implicit none

        if(If_parallel == 1)then
          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
          CALL MPI_FINALIZE(mpi_rc) 
        endif

END SUBROUTINE PARALLEL_Finish







subroutine PARALLEL_define_MPI_BCpairs
    use Global
    implicit none

      integer:: i, idx_temp, idx_tempj, idx_ipairs, idx_ipairs2

      CALL MPI_GATHER(Nmpibcpairs_this, 1, MPI_INTEGER,Nmpibcpairs_all, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr) !���һ�ռ�����������Ϣ������0�ռ�Nmpibcpairs_this��Ϣ;Nmpibcpairs_all��һά���飬ÿ������������ÿ�����̵�parallel P-2-P�߽���

      NmpiBCpairs = 0
      do idx_temp = 1, NumProcs
        NmpiBCpairs = NmpiBCpairs + Nmpibcpairs_all(idx_temp)   !NmpiBCpairs,���е�parallel P-2-P �߽������
      enddo

      CALL MPI_BCAST(NmpiBCpairs,1,MPI_INTEGER, 0,MPI_COMM_WORLD,ierr)! һ�Զ�㲥ͬ����Ϣ������0����NmpiBCpairs��Ϣ

    ! buffer, this
      do idx_temp = 1, NmpiBCpairs_this
      do idx_tempj = 1,7
          Mpi_buffer_this(idx_tempj+7*(idx_temp-1))=Mpi_BCpairs_this(idx_temp, idx_tempj) ! �����̵�buffer��Mpi_buffer_this����Mpi_BCpairs_this��һά��
      enddo
      enddo

      ! displ
      mpi_temp_displ(1)=0
      do idx_temp = 2, NumProcs
         mpi_temp_displ(idx_temp)=mpi_temp_displ(idx_temp-1)+ NmpiBCpairs_all(idx_temp-1) !mpi_temp_displ��n����ǰn-1�����̵�parallel P-2-P�߽�����
        enddo

      do i=1,NumProcs
         mpi_temp_c(i) = Nmpibcpairs_all(i)*7      ! mpi_temp_c(i),��i�����̵�parallel�߽�������7����ʾÿ������Ҫ���͵������С
         mpi_temp_displ(i) = mpi_temp_displ(i)*7   ! mpi_temp_displ(i)��ǰn-1�����̵�parallel P-2-P�߽���������7����������ʱ��ƫ����
      enddo

      CALL MPI_GATHERV(Mpi_buffer_this, NmpiBCpairs_this*7,MPI_INTEGER,Mpi_buffer, mpi_temp_c, mpi_temp_displ, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)!���һ�ռ�����������Ϣ������0�ռ�Mpi_buffer_this��Ϣ;Mpi_buffer��һά���飬����������

      CALL MPI_BCAST(Mpi_buffer,NmpiBCpairs*7,MPI_INTEGER, 0,MPI_COMM_WORLD,ierr) !һ�Զ�㲥ͬ����Ϣ������0����Mpi_buffer��Ϣ
      
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)    !ͬ��

    ! buffer, rev
      do idx_temp = 1, NmpiBCpairs
      do idx_tempj = 1,7
          Mpi_BCpairs(idx_temp, idx_tempj)=Mpi_buffer(idx_tempj+7*(idx_temp-1))  !Mpi_BCpairs������parallel P-2-P �߽����Ϣ
      enddo
      enddo

!! check
!	  write(*,*) 'nmpi'
!	  write(*,*) 'nmpibcpairs=', nmpibcpairs
!	  do idx_temp=1,nmpibcpairs
!	     write(*,*) (mpi_bcpairs(idx_temp,j), j=1,7)
!	  enddo

!	  stop

!	  NMpiTAGs = 2
!	  Mpi_TAGs(1) = 1002001
!	  Mpi_TAGs(2) = 1002002
        NMpiTAGs = 0
      DO idx_ipairs = 1,NmpiBCpairs
      if (Mpi_BCpairs(idx_ipairs, 1) == 0) then !not paired�� ��ʼ������ֵ��Ϊ0
        NMpiTAGs = NMpiTAGs + 1       !���û��paired����NMpiTAGs��ֵ��1
        Mpi_TAGs(NmpiTAGs) = NMpiTAGs + 10000  !dummy ��10000��ʲô���壿
          Mpi_BCpairs(idx_ipairs, 1) = NMpiTAGs + 10000 ! ʲô���壿
        DO idx_ipairs2 = idx_ipairs + 1, NmpiBCpairs 
        if((Mpi_BCpairs(idx_ipairs,2)==Mpi_BCpairs(idx_ipairs2,5)).and.(Mpi_BCpairs(idx_ipairs,3)==Mpi_BCpairs(idx_ipairs2,6)).and.(Mpi_BCpairs(idx_ipairs,4)==Mpi_BCpairs(idx_ipairs2,7))) then
           Mpi_BCpairs(idx_ipairs2, 1) = NMpiTAGs + 10000  !  2=5,3=6,4=7��ͬ���̣�ͬ��ţ�ͬ�߽��?�Լ����Լ������ܳ�����
        endif
        ENDDO
      endif
      ENDDO

end subroutine PARALLEL_define_MPI_BCpairs




subroutine PARALLEL_define_MPI_round
    use Global
    implicit none

	  integer:: nround_all
	  integer:: ita(MaxNumprocs+1)
	  integer:: itb(MaxNumprocs+1)
	  integer:: itx(MaxNumprocs+1)

	  integer:: nround, idx_ipairs, idx_TAGs,idx

	  IF(MOD(numprocs,2)==0)nround_all=numprocs-1
	  If(MOD(numprocs,2)/=0)nround_all=numprocs

	  idx_TAGs=0
	  do nround = 1, nround_all
		CALL PARALLEL_get_MPI_round(numprocs,ita,itb,itx,nround)
	    do idx = 1,numprocs
	      do idx_ipairs=1,NmpiBCpairs
	        if ((ita(idx)<itb(idx)).and. (Mpi_BCpairs(idx_ipairs,2)==ita(idx)).and. (Mpi_BCpairs(idx_ipairs,5)==itb(idx))) then
	            idx_TAGs=idx_TAGs+1
	            Mpi_TAGs(idx_TAGs)=Mpi_BCpairs(idx_ipairs,1)
	        endif
	      enddo
	    enddo
	  enddo

	!check
	    if (idx_TAGs .ne. NmpiTAGs) then
			write(*,*) 'error MPI round sort process'
			stop
		endif
	 
	 end subroutine PARALLEL_define_MPI_round

	

subroutine PARALLEL_get_MPI_round(numprocs,ita,itb,itx,nround)
implicit none
	integer:: numprocs, nround
	integer:: ita(numprocs+1),itb(numprocs+1),itx(numprocs+1)
	
	integer:: i
	
	if(mod(numprocs,2)==0)then
		if(nround==1)then
			do i=1,numprocs
				ita(i)=i-1
			end do
		else
			ita(1)=itx(1)
			ita(2)=itx(numprocs)
			do i=3,numprocs
				ita(i)=itx(i-1)
			end do
		endif
		do i=1,numprocs
			itb(i)=ita(numprocs-i+1)
		end do
		itx=ita
	else
		if(nround==1)then
			do i=1,numprocs+1
				ita(i)=i-1
			end do
		else
			ita(1)=itx(1)
			ita(2)=itx(numprocs+1)
			do i=3,numprocs+1
				ita(i)=itx(i-1)
			end do
		endif
		do i=1,numprocs+1
			itb(i)=ita(numprocs+1-i+1)
		end do
		itx=ita
	endif
	return

	end subroutine PARALLEL_get_MPI_round



subroutine PARALLEL_define_BC_getids_blk(iBlock)   !Ѱ�ұ����Ӧ�ĶԷ����ĺš���š��߽��
    use Global
    implicit none
    integer:: iBlock

      integer:: Nblkbc
      integer:: i_temp_blkex
      integer:: l_foundblk
      
      integer:: Num_Block_BC 
      !Max_Block
      integer:: id_present_blk
      
      INTEGER, POINTER:: IJKBC(:,:)
      
      CALL GLOBAL_SetPointers(iBlock)
      
      Num_Block_BC = ThisBlock%Num_Block_BC
      ! Max_Block
      id_present_blk = ThisBlock%id_present_blk
      IJKBC => ThisBlock%IJKBC

      do Nblkbc = 1, Num_Block_BC
        if (ijkbc(Nblkbc, i_kind) == 9) then
	    i_temp_blkex = ijkbc(Nblkbc, i_blk_ex)
	    l_foundblk = 0
	    !do i_temp_Mblk = 1, Max_Block !N_blocks
	    !  !if (M_blockids(N_blocks) == i_temp_blkex) then ! in this thread
	    !  if (M_blockids(i_temp_Mblk) == i_temp_blkex) then ! in this thread
	    !     l_foundblk = 1
	    !     exit
	    !  endif
	    !enddo
        if (ProcNo_All(i_temp_blkex) == MyID) then
            l_foundblk = 1
        endif

        if (l_foundblk == 0) then !should do MPI
            Nmpibcpairs_this = Nmpibcpairs_this + 1             ! ÿ��һ��parallel P-2-P����Nmpibcpairs_this��1
	        Mpi_BCpairs_this(Nmpibcpairs_this,1)=0  !tag        ����ʱΪ0��tag��
	        Mpi_BCpairs_this(Nmpibcpairs_this,2)=myid  !id0     ������ļ�����ĺţ�id
	        Mpi_BCpairs_this(Nmpibcpairs_this,3)=id_present_blk  !blk0      ��������
	        Mpi_BCpairs_this(Nmpibcpairs_this,4)=Nblkbc  !bc0               ��parallel P-2-P�ı߽��
	        Mpi_BCpairs_this(Nmpibcpairs_this,5)=ProcNo_All(i_temp_blkex) !ijkBC(Nblkbc,i_kindsub) !-1  !id1???  �Է���ļ�����ĺ�
	        Mpi_BCpairs_this(Nmpibcpairs_this,6)=i_temp_blkex !ijkBC(Nblkbc,i_blk_ex)  !blk1   �Է����
	        Mpi_BCpairs_this(Nmpibcpairs_this,7)=ijkBC(Nblkbc,i_blkbd_ex)  !bc1    �Է��߽�ı��
            
            ijkbc(Nblkbc, i_kind) = 10       !���߽����͸�Ϊparallel P-2-P
            ijkbc(Nblkbc, i_kindsub) = ProcNo_All(i_temp_blkex)  ! �����Ϊ�Է�������ĺ�
        write(*,*)"MPI=",ThisBlock%ID_Present_Blk,Nblkbc
	    endif

	  endif
	  enddo
END subroutine PARALLEL_define_BC_getids_blk