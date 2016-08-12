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
  int SS;
  struct tri_listnode* next;
};

struct tri_listnode* lonely_tris;
size_t lonely_tris_num;

FILE* f2;
size_t doubleBC;

void checkTri(int idx1, int idx2, int idx3, int SS, double* coords);

int main(int argc, char* argv[]) {
  if (argc != 2) {
    printf("Usage: findBadBC mesh.ncdf\n");
    exit(1);
  }

  printf("Opening '%s'...\n", argv[1]);

  /* This will be the netCDF ID for the file, dimensions, and data variables. */
  int ncid;
  int ncoords_dim,tetexterior_dim,tetexteriorsize_dim;

  //Dimensions
  size_t ncoords,tetexterior,tetexteriorsize;

  //Data
  int tet_ext_id;
  nc_type tet_ext_type;
  int* tet_ext;

  int coords_id;
  //nc_type coords_type;
  double* coords;
  
  
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

  /* Get the varid of the data variable, based on its name. */
  if ((retval = nc_inq_varid(ncid, "tetrahedron_exterior", &tet_ext_id)))
    ERR(retval);
  if ((retval = nc_inq_var(ncid, tet_ext_id, NULL, &tet_ext_type,NULL,NULL,NULL)))
    ERR(retval);
  //TODO: Check that type is correct...
  tet_ext = malloc(sizeof(int)*tetexterior*tetexteriorsize);
  if ((retval = nc_get_var_int(ncid, tet_ext_id, tet_ext)))
    ERR(retval);

  if ((retval = nc_inq_varid(ncid, "coords", &coords_id)))
    ERR(retval);
  coords = malloc(sizeof(double)*ncoords*3);
  if ((retval = nc_get_var_double(ncid, coords_id, coords)))
    ERR(retval);

  //Initialize linked lists
  lonely_tris = NULL;
  
  // Loop over all tris to find pairs
  doubleBC = 0;
  f2 = fopen("findBadBC_doubleBC.csv","w");
  if (f2==NULL) { printf("error opening file findBadBC_doubleBC.csv!\n"); exit(1);}
  fprintf(f2, "x coord,y coord,z coord,scalar\n");
  
  for (int i = 0; i < tetexterior; i++) {
    size_t idx_off = i*tetexteriorsize;
    //printf("%zu %i %i %i %i %i\n", idx_off, tet_ext[idx_off+1], tet_ext[idx_off+2], tet_ext[idx_off+3], tet_ext[idx_off+4], tet_ext[idx_off+8]);

    checkTri(tet_ext[idx_off+1],tet_ext[idx_off+2],tet_ext[idx_off+3], tet_ext[idx_off+5], coords);
    checkTri(tet_ext[idx_off+1],tet_ext[idx_off+2],tet_ext[idx_off+4], tet_ext[idx_off+6],coords);
    checkTri(tet_ext[idx_off+1],tet_ext[idx_off+3],tet_ext[idx_off+4], tet_ext[idx_off+7],coords);
    checkTri(tet_ext[idx_off+2],tet_ext[idx_off+3],tet_ext[idx_off+4], tet_ext[idx_off+8], coords);
  }
  fclose(f2);
  
  //Find the lonely tris with SS=-1:
  printf("Writing file for noBC badBCs...\n");
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


  printf("lonely_noBC=%zu\n",lonely_noBC);
  printf("doubleBC=%zu\n",doubleBC);
  
  /* Close the file, freeing all resources. */
  if ((retval = nc_close(ncid)))
    ERR(retval);
  
  free(tet_ext);
  tet_ext=NULL;
}


void checkTri(int idx1, int idx2, int idx3, int SS, double* coords) {
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

    //printf("%i %i %i\n",i1,i2,i3);
  }

  /* 
    //Debug
    if (!(i1<i2 && i2<i3)){
    printf("ERR\n");
    exit(1);
  }
  */

  //Compare the current tri with all known lonely tris:
  struct tri_listnode* lonely_tris_curr = lonely_tris;
  struct tri_listnode** lonely_tris_currLink = &lonely_tris; //Pointer to the pointer to the current element (i.e. inside the previous element or the root pointer)
  while (lonely_tris_curr != NULL) {
    if (lonely_tris_curr->i1 == i1 &&
	lonely_tris_curr->i2 == i2 &&
	lonely_tris_curr->i3 == i3    ) {
      //printf("MATCH!\n");

      //Check if there is a non-internal BC set
      if(SS != -1 || lonely_tris_curr->SS != -1) {
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
  if (lonely_tris_curr == NULL) { //There was no match. Add a new tri.
    lonely_tris_curr=*lonely_tris_currLink = (struct tri_listnode*) malloc(sizeof(struct tri_listnode));
    lonely_tris_curr->i1=i1;
    lonely_tris_curr->i2=i2;
    lonely_tris_curr->i3=i3;
    lonely_tris_curr->SS=SS;
    lonely_tris_curr->next=NULL;
    lonely_tris_num++;
    //printf("%zu\n",lonely_tris_num);
  }
  
}
