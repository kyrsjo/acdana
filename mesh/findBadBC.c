#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

//Lists to hold temp data
struct tri_listnode {
  int i1;
  int i2;
  int i3;
  int inTet;
  int SS;
  struct tri_listnode* next;
};

struct tri_listnode* lonely_tris;
size_t lonely_tris_num;

FILE* f2;
size_t doubleBC;

void checkTri(int idx1, int idx2, int idx3, int SS, int inTet);
int findInternal(int i1, int i2, int i3);

//Dimensions
size_t ncoords,
  tetexterior,tetexteriorsize,
  tetinterior,tetinteriorsize;

//Data
int tet_ext_id;
nc_type tet_ext_type;
int* tet_ext;

int coords_id;
//nc_type coords_type;
double* coords;

int tet_int_id;
//nc_type tet_int_type;
int* tet_int;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    printf("Usage: findBadBC mesh.ncdf\n");
    exit(1);
  }

  printf("Opening '%s'...\n", argv[1]);

  /* This will be the netCDF ID for the file, dimensions, and data variables. */
  int ncid;
  int ncoords_dim,
    tetexterior_dim,tetexteriorsize_dim,
    tetinterior_dim, tetinteriorsize_dim;

  // (the actual data is global)
  
  /* error handling. */
  int retval;
  
  /* Open the file. NC_NOWRITE tells netCDF we want read-only access
   * to the file.*/
  if ((retval = nc_open(argv[1], NC_NOWRITE, &ncid)))
    ERR(retval);

  if ((retval = nc_inq_dimid(ncid, "ncoords", &ncoords_dim)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, ncoords_dim, &ncoords)))
    ERR(retval);
  printf("ncoords=%zu\n",ncoords);

  if ((retval = nc_inq_dimid(ncid, "tetexterior", &tetexterior_dim)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, tetexterior_dim, &tetexterior)))
    ERR(retval);
  printf("tetexterior=%zu\n",tetexterior);
  
  if ((retval = nc_inq_dimid(ncid, "tetexteriorsize", &tetexteriorsize_dim)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, tetexteriorsize_dim, &tetexteriorsize)))
    ERR(retval);
  printf("tetexteriorsize=%zu\n",tetexteriorsize);
  if (tetexteriorsize != 9) {
    printf("Expected tetexteriorsize==\n9");
    exit(1);
  }

  if ((retval = nc_inq_dimid(ncid, "tetinterior", &tetinterior_dim)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, tetinterior_dim, &tetinterior)))
    ERR(retval);
  printf("tetinterior=%zu\n",tetinterior);
  
  if ((retval = nc_inq_dimid(ncid, "tetinteriorsize", &tetinteriorsize_dim)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncid, tetinteriorsize_dim, &tetinteriorsize)))
    ERR(retval);
  printf("tetinteriorsize=%zu\n",tetinteriorsize);
  if (tetinteriorsize != 5) {
    printf("Expected tetinteriorsize==5\n");
    exit(1);
  }
  
  printf("Loading mesh data...");
  
  if ((retval = nc_inq_varid(ncid, "tetrahedron_exterior", &tet_ext_id)))
    ERR(retval);
  //if ((retval = nc_inq_var(ncid, tet_ext_id, NULL, &tet_ext_type,NULL,NULL,NULL)))
  //  ERR(retval);
  //TODO: Check that type is correct...
  tet_ext = malloc(sizeof(int)*tetexterior*tetexteriorsize);
  if ((retval = nc_get_var_int(ncid, tet_ext_id, tet_ext)))
    ERR(retval);

  if ((retval = nc_inq_varid(ncid, "coords", &coords_id)))
    ERR(retval);
  //TODO: Check that type is correct...
  coords = malloc(sizeof(double)*ncoords*3);
  if ((retval = nc_get_var_double(ncid, coords_id, coords)))
    ERR(retval);

  if ((retval = nc_inq_varid(ncid, "tetrahedron_interior", &tet_int_id)))
    ERR(retval);
  //if ((retval = nc_inq_var(ncid, tet_int_id, NULL, &tet_int_type,NULL,NULL,NULL)))
  //  ERR(retval);
  //TODO: Check that type is correct...
  tet_int = malloc(sizeof(int)*tetinterior*tetinteriorsize);
  if ((retval = nc_get_var_int(ncid, tet_int_id, tet_int)))
    ERR(retval);

  printf(" Done.\n");
  
  //Initialize linked lists
  lonely_tris = NULL;
  
  // Loop over all tris to find pairs
  doubleBC = 0;
  f2 = fopen("findBadBC_doubleBC.csv","w");
  if (f2==NULL) { printf("error opening file findBadBC_doubleBC.csv!\n"); exit(1);}
  fprintf(f2, "x coord,y coord,z coord,scalar\n");

  printf("Looping over exterior tris, calling checkTri() on each.\n This checks if it is seen before (=> it is already in lonely_tris list), if so check that they both are in SS=-1.\n If it matches, delete the old match from lonely_tris.\n It further checks if it matches ANY of the interior tris, if so it should have SS=-1.\n If it is not matching any other tris, then push it to the lonely_tris list.\n");
  for (int i = 0; i < tetexterior; i++) {
    if (i%100 == 0) printf("completed=%6.2f%% len(lonely_tris)=%8zu\n", 100*((double)i)/((double)tetexterior),lonely_tris_num);
    size_t idx_off = i*tetexteriorsize;
    //printf("%i %i %i %i %i %i\n", i, tet_ext[idx_off+1], tet_ext[idx_off+2], tet_ext[idx_off+3], tet_ext[idx_off+4], tet_ext[idx_off+8]);

    checkTri(tet_ext[idx_off+1],tet_ext[idx_off+2],tet_ext[idx_off+3], tet_ext[idx_off+5], i);
    checkTri(tet_ext[idx_off+1],tet_ext[idx_off+2],tet_ext[idx_off+4], tet_ext[idx_off+7], i);
    checkTri(tet_ext[idx_off+1],tet_ext[idx_off+3],tet_ext[idx_off+4], tet_ext[idx_off+6], i);
    checkTri(tet_ext[idx_off+2],tet_ext[idx_off+3],tet_ext[idx_off+4], tet_ext[idx_off+8], i);

  }
  fclose(f2);
  printf(" Done.\n");
  
  //Find the lonely tris with SS=-1:
  printf("Now looking for lonely (surface) tris with SS==-1...\n");
  FILE* f = fopen("findBadBC_noBC.csv","w");
  if (f==NULL) { printf("error opening file findBadBC_noBC.csv!\n"); exit(1);}
  fprintf(f, "x coord,y coord,z coord,scalar\n");
  
  struct tri_listnode* lonely_tris_curr = lonely_tris;
  size_t lonely_noBC = 0;
  while (lonely_tris_curr != NULL) {
    if (lonely_tris_curr->SS == -1) {
      //3 idxs + coordinates of first corner
      //fprintf(f,"%i %i %i %g %g %g\n", lonely_tris_curr->i1, lonely_tris_curr->i2, lonely_tris_curr->i3, coords[lonely_tris_curr->i1*3], coords[lonely_tris_curr->i1*3+1], coords[lonely_tris_curr->i1*3+2]);
      // CSV, can be loaded into paraview
      // http://www.paraview.org/Wiki/ParaView/Data_formats
      fprintf(f,"%.10e,%.10e,%.10e,1\n", coords[lonely_tris_curr->i1*3], coords[lonely_tris_curr->i1*3+1], coords[lonely_tris_curr->i1*3+2]);
      fprintf(f,"%.10e,%.10e,%.10e,1\n", coords[lonely_tris_curr->i2*3], coords[lonely_tris_curr->i2*3+1], coords[lonely_tris_curr->i2*3+2]);
      fprintf(f,"%.10e,%.10e,%.10e,1\n", coords[lonely_tris_curr->i3*3], coords[lonely_tris_curr->i3*3+1], coords[lonely_tris_curr->i3*3+2]);
      lonely_noBC++;
    }
    
    //Advance to the next element
    lonely_tris_curr=lonely_tris_curr->next;
  }
  fclose(f);
  printf(" Done.\n");
  
  printf("lonely_noBC=%zu\n",lonely_noBC);
  printf("doubleBC=%zu\n",doubleBC);
  
  if (lonely_noBC==0 && doubleBC==0)
    printf("Mesh seems good :)\n");
  else
    printf("Bad mesh -- use Paraview to see where the problems are :(\n");
  
  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    ERR(retval);
  
  free(tet_ext);
  tet_ext=NULL;

  free(tet_int);
  tet_int=NULL;

  free(coords);
  coords=NULL;

  lonely_tris_curr = lonely_tris;
  struct tri_listnode* lonely_tris_next = lonely_tris;
  while (lonely_tris_curr != NULL) {
    lonely_tris_next = lonely_tris_curr->next;
    free(lonely_tris_curr);
    lonely_tris_num--;
    lonely_tris_curr=lonely_tris_next;
  }
  lonely_tris=NULL;
  lonely_tris_next=NULL;
}


void checkTri(int idx1, int idx2, int idx3, int SS, int inTet) {
  // Make sure the indices are in increasing order
  int i1, i2, i3;
  i1 = idx1; i2=idx2; i3=idx3;  
  if (!(idx1<idx2 && idx2<idx3)){
    //Mini buble sort
    int itmp;
#define SWAP(x,y) {itmp=x;x=y;y=itmp;}
    if (i1>i2) SWAP(i1,i2);
    if (i2>i3) SWAP(i2,i3);
    if (i1>i2) SWAP(i1,i2);
    //if (i2>i3) SWAP(i2,i3);
#undef SWAP
    //printf("%i %i %i\n",i1,i2,i3);
  }
  /*
  //Debug sorting
  if (!(i1<i2 && i2<i3)){
    printf("ERR\n");
    exit(1);
  }
  */
  
  //Compare the current tri with all known lonely tris:
  struct tri_listnode*  lonely_tris_curr     =  lonely_tris; //List element which we are working on
  struct tri_listnode** lonely_tris_currLink = &lonely_tris; //Pointer to the pointer to the current element (i.e. inside the previous element or the root pointer)
  while (lonely_tris_curr != NULL) {
    if (lonely_tris_curr->i1 == i1 &&
	lonely_tris_curr->i2 == i2 &&
	lonely_tris_curr->i3 == i3    ) {
      
      //Check if there is a non-internal BC set
      if(SS != -1 || lonely_tris_curr->SS != -1) {
	printf("Internal non-lonely %i %i : %i %i %i : %i\n", SS, lonely_tris_curr->SS, i1, i2, i3, lonely_tris_curr->inTet);
	fprintf(f2,"%.10e,%.10e,%.10e,1\n", coords[lonely_tris_curr->i1*3], coords[lonely_tris_curr->i1*3+1], coords[lonely_tris_curr->i1*3+2]);
	fprintf(f2,"%.10e,%.10e,%.10e,1\n", coords[lonely_tris_curr->i2*3], coords[lonely_tris_curr->i2*3+1], coords[lonely_tris_curr->i2*3+2]);
	fprintf(f2,"%.10e,%.10e,%.10e,1\n", coords[lonely_tris_curr->i3*3], coords[lonely_tris_curr->i3*3+1], coords[lonely_tris_curr->i3*3+2]);
	doubleBC++;
      }
      
      // Delete the matched element
      *lonely_tris_currLink=lonely_tris_curr->next;
      free(lonely_tris_curr);
      lonely_tris_num--;
      break;
    }
    //Advance to the next element
    lonely_tris_currLink=&(lonely_tris_curr->next);
    lonely_tris_curr=lonely_tris_curr->next;
  }

  int matchInternal = findInternal(i1,i2,i3);
  if (matchInternal != -1 && SS != -1) {
    printf("Bad tri: SS=%i and it matches an internal tri.\n", SS);
  }
  
  if (lonely_tris_curr == NULL && matchInternal==-1) {
    //There was no match. Add a new tri to the list of lonelies.
    lonely_tris_curr=*lonely_tris_currLink = (struct tri_listnode*) malloc(sizeof(struct tri_listnode));
    lonely_tris_curr->i1=i1;
    lonely_tris_curr->i2=i2;
    lonely_tris_curr->i3=i3;
    lonely_tris_curr->inTet=inTet;
    lonely_tris_curr->SS=SS;
    lonely_tris_curr->next=NULL;
    lonely_tris_num++;
    //printf("%zu\n",lonely_tris_num);
  }
  
}

int matchInternal(int i1,   int i2,   int i3,
		  int idx1, int idx2, int idx3) {
  //Sort idx1,idx2,idx3, then compare it to the (assumed to be sorted)
  // i1,i2,i3. If they match, return 1, otherwise return 0.
  
  // Make sure the indices are in increasing order
  int j1, j2, j3;
  j1 = idx1; j2=idx2; j3=idx3;  
  if (!(j1<j2 && j2<j3)){
    //Mini buble sort
    int jtmp;
#define SWAP(x,y) {jtmp=x;x=y;y=jtmp;}
    if (j1>j2) SWAP(j1,j2);
    if (j2>j3) SWAP(j2,j3);
    if (j1>j2) SWAP(j1,j2);
#undef SWAP
  }
  /*
  //Debug sorting
  if (!(j1<j2 && j2<j3)){
    printf("ERR\n");
    exit(1);
  }
  */

  if (i1==j1 && i2==j2 && i3==j3) {
    //It's a match with an internal TRI
    return 1;
  }
  return 0;
}
int findInternal(int i1, int i2, int i3) {
  //Return the index of an internal tet containing a tri matching (i1,i2,i3),
  // otherwise return -1.
  for(int i = 0; i<tetinterior;i++) {
    size_t idx_off = i*tetinteriorsize;
    if (matchInternal(i1,i2,i3, tet_int[idx_off+1],tet_int[idx_off+2],tet_int[idx_off+3])) return i;
    if (matchInternal(i1,i2,i3, tet_int[idx_off+1],tet_int[idx_off+2],tet_int[idx_off+4])) return i;
    if (matchInternal(i1,i2,i3, tet_int[idx_off+1],tet_int[idx_off+3],tet_int[idx_off+4])) return i;
    if (matchInternal(i1,i2,i3, tet_int[idx_off+2],tet_int[idx_off+3],tet_int[idx_off+4])) return i;
  }
  return -1; //Nothing found
}
