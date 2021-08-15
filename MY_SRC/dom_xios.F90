MODULE dom_xios
# if defined key_iomput
   USE dom_oce

   IMPLICIT NONE

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  lon_grid_T, lat_grid_T, area_grid_T
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  lon_grid_U, lat_grid_U, area_grid_U
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  lon_grid_V, lat_grid_V, area_grid_V
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  lon_grid_W, lat_grid_W, area_grid_W

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  bounds_lon_grid_T, bounds_lat_grid_T
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  bounds_lon_grid_U, bounds_lat_grid_U
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  bounds_lon_grid_V, bounds_lat_grid_V
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  bounds_lon_grid_W, bounds_lat_grid_W

   INTEGER, PUBLIC, SAVE :: n_ibegin, n_ni
   INTEGER, PUBLIC, SAVE :: n_jbegin, n_nj
   INTEGER, PUBLIC, SAVE :: n_data_ibegin, n_data_ni
   INTEGER, PUBLIC, SAVE :: n_data_jbegin, n_data_nj
   LOGICAL, PUBLIC, SAVE :: using_xios_coordinates=.FALSE.

CONTAINS


  SUBROUTINE init_dom_xios(iin,ijn,iimppt,ijmppt,ildi,ildj,ilei,ilej)
    INTEGER,INTENT(IN) :: iin(jpnij), ijn(jpnij)
    INTEGER,INTENT(IN) :: iimppt(jpni,jpnj), ijmppt(jpni,jpnj)
    INTEGER,INTENT(IN) :: ildi(jpni,jpnj), ildj(jpni,jpnj)     
    INTEGER,INTENT(IN) :: ilei(jpni,jpnj), ilej(jpni,jpnj)
    INTEGER :: rank(jpni,jpnj)
    LOGICAL :: tag_x(jpni,jpnj)
    LOGICAL :: tag_y(jpni,jpnj)
    INTEGER :: jproc,i,j
    INTEGER :: iend,jend
    
    rank(:,:)=-1
    tag_x(:,:)=.FALSE.
    tag_y(:,:)=.FALSE.
    
    DO jproc = 1, jpnij    
      rank(iin(jproc),ijn(jproc))=jproc
    ENDDO
    

!! Distribute holes to neighbour domains
  
    DO j=2,jpnj
      DO i=1,jpni
        IF (rank(i,j)==-1 .AND. rank(i,j-1)/=-1 .AND. .NOT. tag_x(i,j-1)) THEN 
          rank(i,j)=rank(i,j-1)
          tag_y(i,j)=.TRUE.      
          tag_y(i,j-1)=.TRUE.      
        ENDIF
      ENDDO
    ENDDO

    DO j=jpnj-1,1,-1
      DO i=1,jpni
        IF (rank(i,j)==-1 .AND. rank(i,j+1)/=-1 .AND. .NOT. tag_x(i,j+1)) THEN 
          rank(i,j)=rank(i,j+1)      
          tag_y(i,j)=.TRUE.      
          tag_y(i,j+1)=.TRUE.      
        ENDIF
      ENDDO
    ENDDO

    DO j=1,jpnj
      DO i=2,jpni
        IF (rank(i,j)==-1 .AND. rank(i-1,j)/=-1 .AND. .NOT. tag_y(i-1,j)) THEN 
          rank(i,j)=rank(i-1,j)      
          tag_x(i,j)=.TRUE. 
          tag_x(i-1,j)=.TRUE.
        ENDIF
      ENDDO
    ENDDO

    DO j=1,jpnj
      DO i=jpni-1,1,-1
        IF (rank(i,j)==-1 .AND. rank(i+1,j)/=-1 .AND. .NOT. tag_y(i+1,j)) THEN 
          rank(i,j)=rank(i+1,j)      
          tag_x(i,j)=.TRUE. 
          tag_x(i+1,j)=.TRUE.
        ENDIF
      ENDDO
    ENDDO
    
    

!!!! compute new domain decomposition for xios
    n_ibegin=jpiglo
    n_jbegin=jpjglo
    iend=-1
    jend=-1
    DO j=1,jpnj
      DO i=1,jpni
        IF (rank(i,j)==narea) THEN
          n_ibegin=min(n_ibegin,iimppt(i,j)+ildi(i,j)-1)
          iend=max(iend,iimppt(i,j)+ilei(i,j)-1)
          n_jbegin=min(n_jbegin,ijmppt(i,j)+ildj(i,j)-1)
          jend=max(jend,ijmppt(i,j)+ilej(i,j)-1)
        ENDIF  
      ENDDO
    ENDDO
    
    n_ni=iend-n_ibegin+1
    n_nj=jend-n_jbegin+1
    
    n_data_ibegin=nimpp-n_ibegin
    n_data_ni=jpi
    
    n_data_jbegin=njmpp-n_jbegin
    n_data_nj=jpj

    ALLOCATE(lon_grid_T(n_ni,n_nj), lat_grid_T(n_ni,n_nj),area_grid_T(n_ni,n_nj))
    ALLOCATE(lon_grid_U(n_ni,n_nj), lat_grid_U(n_ni,n_nj),area_grid_U(n_ni,n_nj))
    ALLOCATE(lon_grid_V(n_ni,n_nj), lat_grid_V(n_ni,n_nj),area_grid_V(n_ni,n_nj))
    ALLOCATE(lon_grid_W(n_ni,n_nj), lat_grid_W(n_ni,n_nj),area_grid_W(n_ni,n_nj))
    ALLOCATE(bounds_lon_grid_T(4,n_ni,n_nj), bounds_lat_grid_T(4,n_ni,n_nj))
    ALLOCATE(bounds_lon_grid_U(4,n_ni,n_nj), bounds_lat_grid_U(4,n_ni,n_nj))
    ALLOCATE(bounds_lon_grid_V(4,n_ni,n_nj), bounds_lat_grid_V(4,n_ni,n_nj))
    ALLOCATE(bounds_lon_grid_W(4,n_ni,n_nj), bounds_lat_grid_W(4,n_ni,n_nj))
      
  END SUBROUTINE init_dom_xios

#else

CONTAINS
   USE par_oce
   SUBROUTINE init_dom_xios(iin,ijn,iimppt,ijmppt,ildi,ildj,ilei,ilej)
    INTEGER,INTENT(IN) :: iin(jpnij), ijn(jpnij)
    INTEGER,INTENT(IN) :: iimppt(jpni,jpnj), ijmppt(jpni,jpnj)
    INTEGER,INTENT(IN) :: ildi(jpni,jpnj), ildj(jpni,jpnj)     
    INTEGER,INTENT(IN) :: ilei(jpni,jpnj), ilej(jpni,jpnj)
   END SUBROUTINE init_dom_xios
   
#endif
  
END MODULE dom_xios
