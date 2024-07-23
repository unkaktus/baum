/* outputXDMF.c */
/* Anna Neuweiler 3/2022 */

#include "bam.h"
#include "output.h"

#ifdef USEHDF5
#include "hdf5.h"  // this library must be included before Xdmf
#endif

#ifdef XDMF
#include "XdmfAttribute.hpp"
#include "XdmfAttributeType.hpp"
#include "XdmfAttributeCenter.hpp"
#include "XdmfDomain.hpp"
#include "XdmfGeometry.hpp"
#include "XdmfGeometryType.hpp"
#include "XdmfGridCollection.hpp"
#include "XdmfGridCollectionType.hpp"
#include "XdmfHDF5Writer.hpp"
#include "XdmfHDF5Controller.hpp"
#include "XdmfReader.hpp"
#include "XdmfRegularGrid.hpp"
#include "XdmfRectilinearGrid.hpp"
#include "XdmfCurvilinearGrid.hpp"
#include "XdmfTime.hpp"
#include "XdmfTopology.hpp"
#include "XdmfWriter.hpp"
#endif


#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"


/**
 * @brief sanitise_time_steps_curviilinear
 *
 *   Remove all grids with a time value bigger than or equal to the given `time`
 *   from the collection.
 *   This is important when a run is aborted and later continued from a checkpoint
 *   that starts from an earlier time
 *
 * @param collection - collection of curvilinear grids
 * @param time - upper limit for the time value in the collection
 */
#ifdef XDMF
void sanitise_time_steps_curvilinear(XDMFGRIDCOLLECTION *collection,
                                     const double time) {
  // Unlinking data sets might provide new storage space within the hdf5 file,
  // but it will also make the old data inaccessible
  const int unlink_hdf5_data = 0;
  size_t ngrids = XdmfGridCollectionGetNumberCurvilinearGrids(collection);

  if(ngrids>0) {
    XDMFCURVILINEARGRID *last_grid =
        XdmfGridCollectionGetCurvilinearGrid(collection,ngrids-1);

    XDMFTIME *last_time = XdmfCurvilinearGridGetTime(last_grid);

    while(time <= XdmfTimeGetValue(last_time)) {
      //printf("duplicate time %f\n",t);

      XDMFATTRIBUTE *data_attribute_of_last_grid =
          XdmfCurvilinearGridGetAttribute(last_grid,0);

      XDMFHDF5CONTROLLER *heavy_controller =
          (XDMFHDF5CONTROLLER*) XdmfAttributeGetHeavyDataController(
            data_attribute_of_last_grid,0);

      const char *hdf5_file_name = XdmfHDF5ControllerGetFilePath(heavy_controller);
      const char *hdf5_dataset_name = XdmfHDF5ControllerGetDataSetPath(heavy_controller);

      XDMFGEOMETRY *geometry_of_last_grid = XdmfCurvilinearGridGetGeometry(last_grid);

      XDMFHDF5CONTROLLER *geometry_heavy_controller =
          (XDMFHDF5CONTROLLER*) XdmfGeometryGetHeavyDataController(
            geometry_of_last_grid,0);

      const char *geometry_hdf5_file_name =
          XdmfHDF5ControllerGetFilePath(geometry_heavy_controller);
      const char *geometry_hdf5_dataset_name =
          XdmfHDF5ControllerGetDataSetPath(geometry_heavy_controller);


      XdmfGridCollectionRemoveCurvilinearGrid(collection, ngrids-1);
      last_time = NULL; // memory is freed by XdmfGridCollectionRemove...
      last_grid = NULL;
      data_attribute_of_last_grid = NULL;
      heavy_controller = NULL;
      geometry_of_last_grid = NULL;
      geometry_heavy_controller = NULL;


      /*
       * The heavy data persists and there seems to be no way to delete it from the
       * API of the Xdmf library.
       * HDF5 v 1.8 does not have this feature, so it makes sense, that Xdmf cannot
       * do it. Instead H5Ldelete can delete the data links and then one can run
       * the HDF5 repack tool to create a new HDF5 file without the deleted data.
       * HDF5 v 1.10.1+ also removes the data links, but it will fill the space of the
       * unlinked data with new data, if possible.
       * This feature is enabled when H5F_FSPACE_STRATEGY_FSM_AGGR (default) or
       * H5F_FSPACE_STRATEGY_PAGE file space handling strategy is specified
       * at the time of creation of the HDF5 file. It cannot be changed
       * afterwards!
       * You can check the strategy by running `h5stat path/to/file.h5`.
       * Note: v 1.10 uses a different API for file management strategies.
       * More: https://docs.hdfgroup.org/hdf5/rfc/paged_aggregation.pdf
       *       https://docs.hdfgroup.org/hdf5/rfc/RFC-Page_Buffering.pdf
       *
       * At the moment this does not seem to help much. Maybe one has
       * to make the file manager persistent across multiple reopenings of
       * the file.
       */

      if(unlink_hdf5_data) {
#ifdef USEHDF5
        hid_t hdf5_file = H5Fopen(hdf5_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
        if (H5I_INVALID_HID == hdf5_file)
          errorexits("Cannot open file %s",hdf5_file_name);

        if(H5Ldelete(hdf5_file, hdf5_dataset_name, H5P_DEFAULT))
          errorexits("Error when deleting data link in file %s", hdf5_file_name);

        if(H5Fclose(hdf5_file))
          errorexits("Cannot close file %s", hdf5_file_name);

        hid_t geometry_hdf5_file = H5Fopen(geometry_hdf5_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
        if (H5I_INVALID_HID == geometry_hdf5_file)
          errorexits("Cannot open file %s",geometry_hdf5_file_name);

        if(H5Ldelete(geometry_hdf5_file, geometry_hdf5_dataset_name, H5P_DEFAULT))
          errorexits("Error when deleting data link in file %s", geometry_hdf5_file_name);

        if(H5Fclose(geometry_hdf5_file))
          errorexits("Cannot close file %s", geometry_hdf5_file_name);
#else
        char errorstring[1000];
        snprintf(errorstring, 1000, "You must compile with USEHDF5 to unlink the data set "
                                    "%s in %s", hdf5_dataset_name, hdf5_file_name);
#endif
      }

      ngrids = XdmfGridCollectionGetNumberCurvilinearGrids(collection);
      if(ngrids == 0)
        break;
      last_grid =
              XdmfGridCollectionGetCurvilinearGrid(collection,ngrids-1);
      last_time = XdmfCurvilinearGridGetTime(last_grid);
    }
  }
}
#endif

/**
 * @brief sanitise_time_steps_rectilinear
 *
 *   Remove all grids with a time value bigger than or equal to the given `time`
 *   from the collection.
 *   This is important when a run is aborted and later continued from a checkpoint
 *   that starts from an earlier time
 *
 * @param collection - collection of rectilinear grids
 * @param time - upper limit for the time value in the collection
 */
#ifdef XDMF
void sanitise_time_steps_rectilinear(XDMFGRIDCOLLECTION *collection,
                                     const double time) {
  // Unlinking data sets might provide new storage space within the hdf5 file,
  // but it will also make the old data inaccessible
  const int unlink_hdf5_data = 0;
  size_t ngrids = XdmfGridCollectionGetNumberRectilinearGrids(collection);

  if(ngrids>0) {
    XDMFRECTILINEARGRID *last_grid =
        XdmfGridCollectionGetRectilinearGrid(collection,ngrids-1);

    XDMFTIME *last_time = XdmfRectilinearGridGetTime(last_grid);

    while(time <= XdmfTimeGetValue(last_time)) {
      //printf("duplicate time %f\n",t);

      XDMFATTRIBUTE *data_attribute_of_last_grid =
          XdmfRectilinearGridGetAttribute(last_grid,0);

      XDMFHDF5CONTROLLER *heavy_controller =
          (XDMFHDF5CONTROLLER*) XdmfAttributeGetHeavyDataController(
            data_attribute_of_last_grid,0);

      const char *hdf5_file_name = XdmfHDF5ControllerGetFilePath(heavy_controller);
      const char *hdf5_dataset_name = XdmfHDF5ControllerGetDataSetPath(heavy_controller);

      XdmfGridCollectionRemoveRectilinearGrid(collection, ngrids-1);
      last_time = NULL; // memory is freed by XdmfGridCollectionRemove...
      last_grid = NULL;
      data_attribute_of_last_grid = NULL;
      heavy_controller = NULL;

      /*
       * The heavy data persists and there seems to be no way to delete it from the
       * API of the Xdmf library.
       * HDF5 v 1.8 does not have this feature, so it makes sense, that Xdmf cannot
       * do it. Instead H5Ldelete can delete the data links and then one can run
       * the HDF5 repack tool to create a new HDF5 file without the deleted data.
       * HDF5 v 1.10.1+ also removes the data links, but it will fill the space of the
       * unlinked data with new data, if possible.
       * This feature is enabled when H5F_FSPACE_STRATEGY_FSM_AGGR (default) or
       * H5F_FSPACE_STRATEGY_PAGE file space handling strategy is specified
       * at the time of creation of the HDF5 file. It cannot be changed
       * afterwards!
       * You can check the strategy by running `h5stat path/to/file.h5`.
       * Note: v 1.10 uses a different API for file management strategies.
       * More: https://docs.hdfgroup.org/hdf5/rfc/paged_aggregation.pdf
       *       https://docs.hdfgroup.org/hdf5/rfc/RFC-Page_Buffering.pdf
       *
       * At the moment this does not seem to help much. Maybe one has
       * to make the file manager persistent across multiple reopenings of
       * the file.
       */

      if(unlink_hdf5_data) {
#ifdef USEHDF5
        hid_t hdf5_file = H5Fopen(hdf5_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
        if (H5I_INVALID_HID == hdf5_file)
          errorexits("Cannot open file %s",hdf5_file_name);

        if(H5Ldelete(hdf5_file, hdf5_dataset_name, H5P_DEFAULT))
          errorexits("Error when deleting data link in file %s", hdf5_file_name);

        if(H5Fclose(hdf5_file))
          errorexits("Cannot close file %s", hdf5_file_name);
#else
        char errorstring[1000];
        snprintf(errorstring, 1000, "You must compile with USEHDF5 to unlink the data set "
                                    "%s in %s", hdf5_dataset_name, hdf5_file_name);
#endif
      }

      ngrids = XdmfGridCollectionGetNumberRectilinearGrids(collection);
      if(ngrids == 0)
        break;
      last_grid =
              XdmfGridCollectionGetRectilinearGrid(collection,ngrids-1);
      last_time = XdmfRectilinearGridGetTime(last_grid);
    }
  }
}
#endif

void write_xdmf3D_curvilinear(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp)
    mkdir(filename, 0777);
  else
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      XDMFDOMAIN *domain = XdmfDomainNew();
      XDMFGRIDCOLLECTION *collection = XdmfGridCollectionNew();
      XDMFWRITER * writer = XdmfWriterNew(xmfname);
      XDMFHDF5WRITER * heavyWriter = XdmfHDF5WriterNew(h5name, 0);

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);

      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 0);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);

      XdmfDomainFree(domain); domain = NULL;
      XdmfGridCollectionFree(collection); collection = NULL;
      XdmfWriterFree(writer); writer=NULL;
      XdmfHDF5WriterFree(heavyWriter); heavyWriter = NULL;
  }
  else
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE *data_attribute = XdmfAttributeNew();

  /* For curvilinear meshes we are free to choose the order of the points,
   * but Paraview will render the points z running fastest. Hence the
   * data will only look good when ordered in this way.
   */
  for(size_t i=0; i<nx; i++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t k=0; k<nz; k++) {
        if(flt) {
          float float_data = buffer[i + j*nx + k*nx*ny];
          XdmfAttributePushBack(data_attribute, &float_data,
              XDMF_ARRAY_TYPE_FLOAT32, &status);
        }
        else if(dbl) {
          XdmfAttributePushBack(data_attribute, &(buffer[i + j*nx + k*nx*ny]),
            XDMF_ARRAY_TYPE_FLOAT64, &status);
        }
        else
          errorexit("use float or double");
      }
    }
  }

  XdmfAttributeSetType(data_attribute, XDMF_ATTRIBUTE_TYPE_SCALAR, &status);
  XdmfAttributeSetCenter(data_attribute, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(data_attribute, varname, &status);

  /* Create grid */
  XDMFGEOMETRY *geometry = XdmfGeometryNew();
  XdmfGeometrySetType(geometry, XDMF_GEOMETRY_TYPE_XYZ, &status);

  for(size_t i=0; i<nx; i++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t k=0; k<nz; k++) {
        float x = x0 + i*dx;
        float y = y0 + j*dy;
        float z = z0 + k*dz;

        XdmfGeometryPushBack(geometry, &x, XDMF_ARRAY_TYPE_FLOAT32, &status);
        XdmfGeometryPushBack(geometry, &y, XDMF_ARRAY_TYPE_FLOAT32, &status);
        XdmfGeometryPushBack(geometry, &z, XDMF_ARRAY_TYPE_FLOAT32, &status);
      }
    }
  }

  XDMFCURVILINEARGRID *grid = XdmfCurvilinearGridNew3D(nx, ny, nz);
  XdmfCurvilinearGridSetGeometry(grid, geometry, 0);

  /* set time */
  XDMFTIME *time = XdmfTimeNew(t);
  XdmfCurvilinearGridSetTime(grid,time,0);

  /* add data to grid */
  XdmfCurvilinearGridInsertAttribute(grid, data_attribute, 0);

  /* read existing domain from file */
  XDMFREADER *reader = XdmfReaderNew();
  XDMFDOMAIN *domain = (XDMFDOMAIN*)XdmfReaderRead(reader, xmfname, &status);
  XDMFGRIDCOLLECTION *collection = XdmfDomainGetGridCollectionByName(domain, varname);

  sanitise_time_steps_curvilinear(collection, t);

  XdmfGridCollectionInsertCurvilinearGrid(collection, grid, 0);

  /* set xmf writer */
  XDMFWRITER *writer = XdmfWriterNew(xmfname);

  // Write to File
  XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

  XdmfDomainFree(domain); domain = NULL;
  // Segfaults (because it is contained in Domain?):
  //XdmfGridCollectionFree(collection); collection = NULL;

  XdmfAttributeFree(data_attribute); data_attribute = NULL;
  XdmfReaderFree(reader); reader = NULL;
  XdmfTimeFree(time); time = NULL;
  XdmfWriterFree(writer); writer = NULL;

  XdmfGeometryFree(geometry); geometry = NULL;
  XdmfCurvilinearGridFree(grid); grid = NULL;
  #endif
}


void write_xdmf3D_curvilinear_vec(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp)
    mkdir(filename, 0777);
  else
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      XDMFDOMAIN *domain = XdmfDomainNew();
      XDMFGRIDCOLLECTION *collection = XdmfGridCollectionNew();
      XDMFWRITER * writer = XdmfWriterNew(xmfname);
      XDMFHDF5WRITER * heavyWriter = XdmfHDF5WriterNew(h5name, 0);

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);

      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 0);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);

      XdmfDomainFree(domain); domain = NULL;
      XdmfGridCollectionFree(collection); collection = NULL;
      XdmfWriterFree(writer); writer=NULL;
      XdmfHDF5WriterFree(heavyWriter); heavyWriter = NULL;
  }
  else
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE *data_attribute = XdmfAttributeNew();

  /* For curvilinear meshes we are free to choose the order of the points,
   * but Paraview will render the points z running fastest. Hence the
   * data will only look good when ordered in this way.
   */
  for(size_t i=0; i<nx; i++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t k=0; k<nz; k++) {
        for(size_t l=0; l<3; l++) {
          if(flt) {
            float float_data = buffer[i + j*nx + k*nx*ny + l*n];
            XdmfAttributePushBack(data_attribute, &float_data,
              XDMF_ARRAY_TYPE_FLOAT32, &status);
          }
          else if(dbl) {
            XdmfAttributePushBack(data_attribute, &(buffer[i + j*nx + k*nx*ny + l*n]),
              XDMF_ARRAY_TYPE_FLOAT64, &status);
          }
          else
            errorexit("use float or double");
        }
      }
    }
  }

  XdmfAttributeSetType(data_attribute, XDMF_ATTRIBUTE_TYPE_VECTOR, &status);
  XdmfAttributeSetCenter(data_attribute, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(data_attribute, varname, &status);

  /* Create grid */
  XDMFGEOMETRY *geometry = XdmfGeometryNew();
  XdmfGeometrySetType(geometry, XDMF_GEOMETRY_TYPE_XYZ, &status);

  for(size_t i=0; i<nx; i++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t k=0; k<nz; k++) {
        float x = x0 + i*dx;
        float y = y0 + j*dy;
        float z = z0 + k*dz;

        XdmfGeometryPushBack(geometry, &x, XDMF_ARRAY_TYPE_FLOAT32, &status);
        XdmfGeometryPushBack(geometry, &y, XDMF_ARRAY_TYPE_FLOAT32, &status);
        XdmfGeometryPushBack(geometry, &z, XDMF_ARRAY_TYPE_FLOAT32, &status);
      }
    }
  }

  XDMFCURVILINEARGRID *grid = XdmfCurvilinearGridNew3D(nx, ny, nz);
  XdmfCurvilinearGridSetGeometry(grid, geometry, 0);

  /* set time */
  XDMFTIME *time = XdmfTimeNew(t);
  XdmfCurvilinearGridSetTime(grid,time,0);

  /* add data to grid */
  XdmfCurvilinearGridInsertAttribute(grid, data_attribute, 0);

  /* read existing domain from file */
  XDMFREADER *reader = XdmfReaderNew();
  XDMFDOMAIN *domain = (XDMFDOMAIN*)XdmfReaderRead(reader, xmfname, &status);
  XDMFGRIDCOLLECTION *collection = XdmfDomainGetGridCollectionByName(domain, varname);

  sanitise_time_steps_curvilinear(collection, t);

  XdmfGridCollectionInsertCurvilinearGrid(collection, grid, 0);

  /* set xmf writer */
  XDMFWRITER *writer = XdmfWriterNew(xmfname);

  // Write to File
  XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

  XdmfDomainFree(domain); domain = NULL;
  // Segfaults (because it is contained in Domain?):
  //XdmfGridCollectionFree(collection); collection = NULL;

  XdmfAttributeFree(data_attribute); data_attribute = NULL;
  XdmfReaderFree(reader); reader = NULL;
  XdmfTimeFree(time); time = NULL;
  XdmfWriterFree(writer); writer = NULL;

  XdmfGeometryFree(geometry); geometry = NULL;
  XdmfCurvilinearGridFree(grid); grid = NULL;
  #endif
}




void write_xdmf3D_rectilinear(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp)
    mkdir(filename, 0777);
  else
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      XDMFDOMAIN *domain = XdmfDomainNew();
      XDMFGRIDCOLLECTION *collection = XdmfGridCollectionNew();
      XDMFWRITER * writer = XdmfWriterNew(xmfname);
      XDMFHDF5WRITER * heavyWriter = XdmfHDF5WriterNew(h5name, 0);

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);

      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 0);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);

      XdmfDomainFree(domain); domain = NULL;
      XdmfGridCollectionFree(collection); collection = NULL;
      XdmfWriterFree(writer); writer=NULL;
      XdmfHDF5WriterFree(heavyWriter); heavyWriter = NULL;
  }
  else
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE *data_attribute = XdmfAttributeNew();

  // Paraview expects the data to be ordered x running fastest
  for(size_t k=0; k<nz; k++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t i=0; i<nx; i++) {
        if(flt) {
          float float_data = buffer[i + j*nx + k*nx*ny];
          XdmfAttributePushBack(data_attribute, &float_data,
              XDMF_ARRAY_TYPE_FLOAT32, &status);
        }
        else if(dbl) {
          XdmfAttributePushBack(data_attribute, &(buffer[i + j*nx + k*nx*ny]),
            XDMF_ARRAY_TYPE_FLOAT64, &status);
        }
        else
          errorexit("use float or double");
      }
    }
  }

  XdmfAttributeSetType(data_attribute, XDMF_ATTRIBUTE_TYPE_SCALAR, &status);
  XdmfAttributeSetCenter(data_attribute, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(data_attribute, varname, &status);



  /* Create grid */
  XDMFARRAY *xCoordinates = XdmfArrayNew();
  XDMFARRAY *yCoordinates = XdmfArrayNew();
  XDMFARRAY *zCoordinates = XdmfArrayNew();

  for(size_t i=0; i<nx; i++) {
    float x = x0 + i*dx;
    XdmfArrayInsertValue(xCoordinates, i, &x, XDMF_ARRAY_TYPE_FLOAT32, &status);
  }

  for(size_t i=0; i<ny; i++) {
    float y = y0 + i*dy;
    XdmfArrayInsertValue(yCoordinates, i, &y, XDMF_ARRAY_TYPE_FLOAT32, &status);
  }

  for(size_t i=0; i<nz; i++) {
    float z = z0 + i*dz;
    XdmfArrayInsertValue(zCoordinates, i, &z, XDMF_ARRAY_TYPE_FLOAT32, &status);
  }

  XDMFRECTILINEARGRID *grid = XdmfRectilinearGridNew3D(
        xCoordinates, yCoordinates, zCoordinates, 0);

  /* set time */
  XDMFTIME *time = XdmfTimeNew(t);
  XdmfRectilinearGridSetTime(grid,time,0);

  /* add data to grid */
  XdmfRectilinearGridInsertAttribute(grid, data_attribute, 0);

  /* read existing domain from file */
  XDMFREADER *reader = XdmfReaderNew();
  XDMFDOMAIN *domain = (XDMFDOMAIN*)XdmfReaderRead(reader, xmfname, &status);
  XDMFGRIDCOLLECTION *collection = XdmfDomainGetGridCollectionByName(domain, varname);

  sanitise_time_steps_rectilinear(collection, t);

  XdmfGridCollectionInsertRectilinearGrid(collection, grid, 0);

  /* set xmf writer */
  XDMFWRITER *writer = XdmfWriterNew(xmfname);

  // Write to File
  XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

  XdmfDomainFree(domain); domain = NULL;
  // Segfaults (because it is contained in Domain?):
  //XdmfGridCollectionFree(collection); collection = NULL;

  XdmfAttributeFree(data_attribute); data_attribute = NULL;
  XdmfReaderFree(reader); reader = NULL;
  XdmfTimeFree(time); time = NULL;
  XdmfWriterFree(writer); writer = NULL;

  XdmfRectilinearGridFree(grid); grid = NULL;
  XdmfArrayFree(xCoordinates); xCoordinates = NULL;
  XdmfArrayFree(yCoordinates); yCoordinates = NULL;
  XdmfArrayFree(zCoordinates); zCoordinates = NULL;
  #endif
}

void write_xdmf3D_rectilinear_vec(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp)
    mkdir(filename, 0777);
  else
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      XDMFDOMAIN *domain = XdmfDomainNew();
      XDMFGRIDCOLLECTION *collection = XdmfGridCollectionNew();
      XDMFWRITER * writer = XdmfWriterNew(xmfname);
      XDMFHDF5WRITER * heavyWriter = XdmfHDF5WriterNew(h5name, 0);

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);

      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 0);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);

      XdmfDomainFree(domain); domain = NULL;
      XdmfGridCollectionFree(collection); collection = NULL;
      XdmfWriterFree(writer); writer=NULL;
      XdmfHDF5WriterFree(heavyWriter); heavyWriter = NULL;
  }
  else
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE * data_attribute = XdmfAttributeNew();

  // Paraview expects the data to be ordered x running fastest
  for(size_t k=0; k<nz; k++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t i=0; i<nx; i++) {
        for(size_t l=0; l<3; l++) {
          if(flt) {
            float float_data = buffer[i + j*nx + k*nx*ny + l*n];
            XdmfAttributePushBack(data_attribute, &float_data,
              XDMF_ARRAY_TYPE_FLOAT32, &status);
          }
          else if(dbl) {
            XdmfAttributePushBack(data_attribute, &(buffer[i + j*nx + k*nx*ny + l*n]),
              XDMF_ARRAY_TYPE_FLOAT64, &status);
          }
          else
            errorexit("use float or double");
        }
      }
    }
  }

  XdmfAttributeSetType(data_attribute, XDMF_ATTRIBUTE_TYPE_VECTOR, &status);
  XdmfAttributeSetCenter(data_attribute, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(data_attribute, varname, &status);

  /* Create grid */
  XDMFARRAY *xCoordinates = XdmfArrayNew();
  XDMFARRAY *yCoordinates = XdmfArrayNew();
  XDMFARRAY *zCoordinates = XdmfArrayNew();

  for(size_t i=0; i<nx; i++) {
    float x = x0 + i*dx;
    XdmfArrayInsertValue(xCoordinates, i, &x, XDMF_ARRAY_TYPE_FLOAT32, &status);
  }

  for(size_t i=0; i<ny; i++) {
    float y = y0 + i*dy;
    XdmfArrayInsertValue(yCoordinates, i, &y, XDMF_ARRAY_TYPE_FLOAT32, &status);
  }

  for(size_t i=0; i<nz; i++) {
    float z = z0 + i*dz;
    XdmfArrayInsertValue(zCoordinates, i, &z, XDMF_ARRAY_TYPE_FLOAT32, &status);
  }

  XDMFRECTILINEARGRID *grid = XdmfRectilinearGridNew3D(
        xCoordinates, yCoordinates, zCoordinates, 0);

  /* set time */
  XDMFTIME * time = XdmfTimeNew(t);
  XdmfRectilinearGridSetTime(grid,time,0);

  /* add data to grid */
  XdmfRectilinearGridInsertAttribute(grid, data_attribute, 0);

  /* read existing domain from file */
  XDMFREADER *reader = XdmfReaderNew();
  XDMFDOMAIN *domain = (XDMFDOMAIN*)XdmfReaderRead(reader, xmfname, &status);
  XDMFGRIDCOLLECTION *collection = XdmfDomainGetGridCollectionByName(domain, varname);

  sanitise_time_steps_rectilinear(collection, t);

  XdmfGridCollectionInsertRectilinearGrid(collection, grid, 0);

  /* set xmf writer */
  XDMFWRITER *writer = XdmfWriterNew(xmfname);

  // Write to File
  XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

  XdmfDomainFree(domain); domain = NULL;
  // Segfaults (because it is contained in Domain?):
  //XdmfGridCollectionFree(collection); collection = NULL;

  XdmfAttributeFree(data_attribute); data_attribute = NULL;
  XdmfReaderFree(reader); reader = NULL;
  XdmfTimeFree(time); time = NULL;
  XdmfWriterFree(writer); writer = NULL;

  XdmfRectilinearGridFree(grid); grid = NULL;
  XdmfArrayFree(xCoordinates); xCoordinates = NULL;
  XdmfArrayFree(yCoordinates); yCoordinates = NULL;
  XdmfArrayFree(zCoordinates); zCoordinates = NULL;
  #endif
}



void write_xdmf3D(char *varname, char *dirsuffix, char *suffix, int level, double t, 
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp) 
    mkdir(filename, 0777);  
  else 
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf", 
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5", 
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      XDMFDOMAIN *domain = XdmfDomainNew();
      XDMFGRIDCOLLECTION *collection = XdmfGridCollectionNew();
      XDMFWRITER * writer = XdmfWriterNew(xmfname);
      XDMFHDF5WRITER * heavyWriter = XdmfHDF5WriterNew(h5name, 0);

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);

      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 0);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);

      XdmfDomainFree(domain); domain = NULL;
      XdmfGridCollectionFree(collection); collection = NULL;
      XdmfWriterFree(writer); writer=NULL;
      XdmfHDF5WriterFree(heavyWriter); heavyWriter = NULL;
  }
  else 
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE * array = XdmfAttributeNew();

  // Paraview and Visit expect the data to be ordered x running fastest
  for(size_t k=0; k<nz; k++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t i=0; i<nx; i++) {
        if(flt) {
          float float_data = buffer[i + j*nx + k*nx*ny];
          XdmfAttributePushBack(array, &float_data,
              XDMF_ARRAY_TYPE_FLOAT32, &status);
        }
        else if(dbl) {
          XdmfAttributePushBack(array, &(buffer[i + j*nx + k*nx*ny]),
            XDMF_ARRAY_TYPE_FLOAT64, &status);
        }
        else
          errorexit("use float or double");
      }
    }
  }

  XdmfAttributeSetType(array, XDMF_ATTRIBUTE_TYPE_SCALAR, &status);
  XdmfAttributeSetCenter(array, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(array, varname, &status);


  /* write to file */
  // There is a bug in Paraview/VTK. Therefore the values for x and z are switched in this call
  // That means the XDMF that is generated will actually be wrong!
  XDMFREGULARGRID * grid = XdmfRegularGridNew3D(dz, dy, dx, nz, ny, nx, z0, y0, x0);
  // Correct would be:
  // XDMFREGULARGRID * grid = XdmfRegularGridNew3D(dx, dy, dz, nx, ny, nz, x0, y0, z0);

  XDMFTIME * time = XdmfTimeNew(t);
  XdmfRegularGridSetTime(grid,time,0);

  XdmfRegularGridInsertAttribute(grid, array, 0);

  XDMFREADER * reader = XdmfReaderNew();
  //XDMFDOMAIN * readDomain = XdmfDomainNew();
  //XDMFGRIDCOLLECTION * readCollection = XdmfGridCollectionNew();
  XDMFDOMAIN *readDomain = (XDMFDOMAIN*)XdmfReaderRead(reader, xmfname, &status);
  XDMFGRIDCOLLECTION * readCollection = XdmfDomainGetGridCollectionByName(readDomain, varname);
  XdmfGridCollectionInsertRegularGrid(readCollection, grid, 0);

  /* set xmf writer */
  XDMFWRITER *writer = XdmfWriterNew(xmfname);

  // Write to File
  XdmfDomainAccept(readDomain, (XDMFVISITOR *)writer, &status);

  XdmfDomainFree(readDomain); readDomain = NULL;
  // Segfaults (because it is contained in Domain?):
  // XdmfGridCollectionFree(readCollection); readCollection = NULL;

  XdmfAttributeFree(array); array = NULL;
  XdmfReaderFree(reader); reader = NULL;
  XdmfTimeFree(time); time = NULL;
  XdmfWriterFree(writer); writer = NULL;

  XdmfRegularGridFree(grid); grid = NULL;
  #endif

}

void write_xdmf3D_vec(char *varname, char *dirsuffix, char *suffix, int level, double t, 
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp) 
    mkdir(filename, 0777);  
  else 
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf", 
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5", 
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* set xmf writer */
  XDMFWRITER * writer = XdmfWriterNew(xmfname);
  void * heavyWriter = XdmfHDF5WriterNew(h5name, 0);
  int currentindex = XdmfHeavyDataWriterGetFileIndex(heavyWriter);
  XdmfHeavyDataWriterSetFileIndex(heavyWriter, currentindex);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      void * domain = XdmfDomainNew();
      void * collection = XdmfGridCollectionNew();

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);
      
      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 1);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);
  }
  else 
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE * array = XdmfAttributeNew();
  int i,j;

  // Paraview and Visit expect the data to be ordered x running fastest
  for(size_t k=0; k<nz; k++) {
    for(size_t j=0; j<ny; j++) {
      for(size_t i=0; i<nx; i++) {
        for(size_t l=0; l<3; l++) {
          if(flt) {
            float float_data = buffer[i + j*nx + k*nx*ny + l*n];
            XdmfAttributePushBack(array, &float_data,
              XDMF_ARRAY_TYPE_FLOAT32, &status);
          }
          else if(dbl) {
            XdmfAttributePushBack(array, &(buffer[i + j*nx + k*nx*ny + l*n]),
              XDMF_ARRAY_TYPE_FLOAT64, &status);
          }
          else
            errorexit("use float or double");
        }
      }
    }
  }

  XdmfAttributeSetType(array, XDMF_ATTRIBUTE_TYPE_VECTOR, &status);
  XdmfAttributeSetCenter(array, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(array, varname, &status);

  // Write to File
  XDMFWRITER * writerWithHeavy = XdmfWriterNewSpecifyHeavyDataWriter(xmfname, heavyWriter);

  /* write to file */
  void * grid = XdmfRegularGridNew3D(dz, dy, dx, nz, ny, nx, z0, y0, x0);
  void * time = XdmfTimeNew(t);
  XdmfRegularGridSetTime(grid,time,1);

  XdmfRegularGridInsertAttribute(grid, array, 1);

  void * reader = XdmfReaderNew();
  void * readDomain = XdmfDomainNew();
  void * readCollection = XdmfGridCollectionNew();
  readDomain = XdmfReaderRead(reader, xmfname, &status);
  readCollection = XdmfDomainGetGridCollectionByName(readDomain, varname);
  XdmfGridCollectionInsertRegularGrid(readCollection, grid, 1);

  XdmfDomainAccept(readDomain, (XDMFVISITOR *)writerWithHeavy, &status);

  XdmfAttributeFree(array);
  XdmfWriterFree(writer);
  free(reader);
  free(time);
  free(writerWithHeavy);
  #endif
}


void write_xdmf2D(char *varname, char *dirsuffix, char *suffix, int level, double t, 
        int n, double *buffer, int nv, int j, int dbl, int flt, int text, int binary,
        int nx, int ny, double x0, double y0, double dx, double dy)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp)
    mkdir(filename, 0777);
  else
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5",
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* set xmf writer */
  XDMFWRITER * writer = XdmfWriterNew(xmfname);
  void * heavyWriter = XdmfHDF5WriterNew(h5name, 0);
  int currentindex = XdmfHeavyDataWriterGetFileIndex(heavyWriter);
  XdmfHeavyDataWriterSetFileIndex(heavyWriter, currentindex);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      void * domain = XdmfDomainNew();
      void * collection = XdmfGridCollectionNew();

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);

      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 1);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);
  }
  else
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE * array = XdmfAttributeNew();
  int i;

  /* maybe also should account for big endian, text/binary, ... ? */
  if (dbl){
    double xdouble;
    for (i = 0; i < n; i++){
      xdouble = buffer[nv*i+j];
      XdmfAttributePushBack(array, &xdouble , XDMF_ARRAY_TYPE_FLOAT64, &status);
    }
  }
  if (flt){
    float xfloat;
    for (i = 0; i < n; i++){
      xfloat = (float) buffer[nv*i+j];
      XdmfAttributePushBack(array, &xfloat, XDMF_ARRAY_TYPE_FLOAT32, &status);
    }
  }
  else
    errorexit("write_raw(): use float or double");

  XdmfAttributeSetType(array, XDMF_ATTRIBUTE_TYPE_SCALAR, &status);
  XdmfAttributeSetCenter(array, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(array, varname, &status);

  // Write to File
  XDMFWRITER * writerWithHeavy = XdmfWriterNewSpecifyHeavyDataWriter(xmfname, heavyWriter);

  /* write to file */
  void * grid = XdmfRegularGridNew2D(dy, dx, ny, nx, x0, y0);
  void * time = XdmfTimeNew(t);
  XdmfRegularGridSetTime(grid,time,1);

  XdmfRegularGridInsertAttribute(grid, array, 1);

  void * reader = XdmfReaderNew();
  void * readDomain = XdmfDomainNew();
  void * readCollection = XdmfGridCollectionNew();
  readDomain = XdmfReaderRead(reader, xmfname, &status);
  readCollection = XdmfDomainGetGridCollectionByName(readDomain, varname);
  XdmfGridCollectionInsertRegularGrid(readCollection, grid, 1);

  XdmfDomainAccept(readDomain, (XDMFVISITOR *)writerWithHeavy, &status);

  XdmfAttributeFree(array);
  XdmfWriterFree(writer);
  free(reader);
  free(time);
  free(writerWithHeavy);
  #endif
}



void write_xdmf2D_vec(char *varname, char *dirsuffix, char *suffix, int level, double t, 
        int n, double *buffer, int nv, int j, int dbl, int flt, int text, int binary,
        int nx, int ny, double x0, double y0, double dx, double dy, int k, int l)
{
  #ifdef XDMF
  char filename[1000], str[1000], xmfname[1000], h5name[1000];
  sprintf(str,"outdir_%s",dirsuffix);
  char *outdir = Gets(str);
  FILE *fp,*fp2;
  int status = 0;


  /* make sure subdirectory exists */
  snprintf(filename, 1000, "%s/%s.%s%d_xdmf", outdir, varname, suffix, level);
  fp = fopen(filename, "r");
  if (!fp) 
    mkdir(filename, 0777);  
  else 
    fclose(fp);

  /* get file name*/
  snprintf(xmfname, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.xmf", 
           outdir, varname, suffix, level, varname, suffix, level, boxstring);
  snprintf(h5name, 1000, "%s/%s.%s%d_xdmf/%s.%s%d%s.h5", 
           outdir, varname, suffix, level, varname, suffix, level, boxstring);

  /* set xmf writer */
  XDMFWRITER * writer = XdmfWriterNew(xmfname);
  void * heavyWriter = XdmfHDF5WriterNew(h5name, 0);
  int currentindex = XdmfHeavyDataWriterGetFileIndex(heavyWriter);
  XdmfHeavyDataWriterSetFileIndex(heavyWriter, currentindex);

  /* check if file already exists or it is first call */
  fp2 = fopen(xmfname, "r");
  if (!fp2){
      /* create xmf file with important info */
      void * domain = XdmfDomainNew();
      void * collection = XdmfGridCollectionNew();

      // set grid collection type and name
      XdmfGridCollectionSetType(collection, XDMF_GRID_COLLECTION_TYPE_TEMPORAL, &status);
      XdmfGridCollectionSetName(collection, varname, &status);
      
      // write xmf file
      XdmfDomainInsertGridCollection(domain, collection, 1);
      XdmfDomainAccept(domain, (XDMFVISITOR *)writer, &status);

      // write hdf5
      XdmfHDF5WriterOpenFile(heavyWriter, &status);
      XdmfDomainAccept(domain, (XDMFVISITOR *)heavyWriter, &status);
      XdmfHDF5WriterCloseFile(heavyWriter, &status);
  }
  else 
    fclose(fp2);

  /* write raw data as in 'write_raw_vtk' */
  XDMFATTRIBUTE * array = XdmfAttributeNew();
  int i;

  /* maybe also should account for big endian, text/binary, ... ? */
  if (dbl){
    double xdouble;
    for (i = 0; i < n; i++){
      xdouble = buffer[nv*i+j+k];
      XdmfAttributePushBack(array, &xdouble, XDMF_ARRAY_TYPE_FLOAT64, &status);
      xdouble = buffer[nv*i+j+l];
      XdmfAttributePushBack(array, &xdouble, XDMF_ARRAY_TYPE_FLOAT64, &status);
    }
  }
  if (flt){
    float xfloat;
    for (i = 0; i < n; i++){
      xfloat = (float) buffer[nv*i+j+k];
      XdmfAttributePushBack(array, &xfloat, XDMF_ARRAY_TYPE_FLOAT32, &status);
      xfloat = (float) buffer[nv*i+j+l];
      XdmfAttributePushBack(array, &xfloat, XDMF_ARRAY_TYPE_FLOAT32, &status);
    }
  }
  else
    errorexit("write_raw(): use float or double");

  XdmfAttributeSetType(array, XDMF_ATTRIBUTE_TYPE_VECTOR, &status);
  XdmfAttributeSetCenter(array, XDMF_ATTRIBUTE_CENTER_NODE, &status);
  XdmfAttributeSetName(array, varname, &status);

  // Write to File
  XDMFWRITER * writerWithHeavy = XdmfWriterNewSpecifyHeavyDataWriter(xmfname, heavyWriter);

  /* write to file */
  void * grid = XdmfRegularGridNew2D(dy, dx, ny, nx, x0, y0);
  void * time = XdmfTimeNew(t);
  XdmfRegularGridSetTime(grid,time,1);

  XdmfRegularGridInsertAttribute(grid, array, 1);

  void * reader = XdmfReaderNew();
  void * readDomain = XdmfDomainNew();
  void * readCollection = XdmfGridCollectionNew();
  readDomain = XdmfReaderRead(reader, xmfname, &status);
  readCollection = XdmfDomainGetGridCollectionByName(readDomain, varname);
  XdmfGridCollectionInsertRegularGrid(readCollection, grid, 1);

  XdmfDomainAccept(readDomain, (XDMFVISITOR *)writerWithHeavy, &status);

  XdmfAttributeFree(array);
  XdmfWriterFree(writer);
  free(reader);
  free(time);
  free(writerWithHeavy);
  #endif
}
