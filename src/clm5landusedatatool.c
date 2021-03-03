/* Generate Global CLM5 Surface Data from LUH format time series and MODIS and MIRCA2000 current day reference data */
/* Author Peter Lawrence - Terrestrial Sciences Section - National Center for Atmospheric Research */
/* Contact lawrence@ucar.edu 303 - 497 1727 */

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define MAXCLMPIX 1440
#define MAXCLMLIN 720

#define CLMPIXSIZE 0.25
#define CLMLLX -180.0
#define CLMLLY -90.0

#define PI 4.0*atan(1.0)
#define EarthCir 40075.017

#define MAXPFT 15
#define MAXCFTRAW 31
#define MAXCFT 64

#define firsttreepft 1
#define lasttreepft 8

#define RANK_natpft 1
#define RANK_cft 1
#define RANK_EDGEN 0
#define RANK_EDGEE 0
#define RANK_EDGES 0
#define RANK_EDGEW 0
#define RANK_LAT 1
#define RANK_LATIXY 2
#define RANK_LON 1
#define RANK_LONGXY 2
#define RANK_LANDMASK 2
#define RANK_LANDFRAC 2
#define RANK_AREA 2
#define RANK_PCT_GLACIER 2
#define RANK_PCT_LAKE 2
#define RANK_PCT_WETLAND 2
#define RANK_PCT_URBAN 2
#define RANK_PCT_NATVEG 2
#define RANK_PCT_CROP 2
#define RANK_PCT_NAT_PFT 3
#define RANK_PCT_CFT 3
#define RANK_FERTNITRO_CFT 3
#define RANK_HARVEST_VH1 2
#define RANK_HARVEST_VH2 2
#define RANK_HARVEST_SH1 2
#define RANK_HARVEST_SH2 2
#define RANK_HARVEST_SH3 2
#define RANK_GRAZING 2
#define RANK_UNREPRESENTED_PFT_LULCC 3
#define RANK_UNREPRESENTED_CFT_LULCC 3

long MAXOUTPIX = MAXCLMPIX;
long MAXOUTLIN = MAXCLMLIN;
long OUTLONOFFSET = 0;
long OUTLATOFFSET = 0;
float OUTPIXSIZE = CLMPIXSIZE;
float OUTLLX = CLMLLX;
float OUTLLY = CLMLLY;

long OUTDATASIZE = sizeof(float) * MAXCLMPIX * MAXCLMLIN;
long OUTDBLDATASIZE = sizeof(double) * MAXCLMPIX * MAXCLMLIN;

/* Namelist Variables */

char regionfilename[1024];
char outputdir[1024];
char outputseries[1024];
int firstrefyear;
int refyear;
int firstyear;
int startyear;
int endyear;
char clmcurrentsurfdb[1024];
char clmLUHforestdb[1024];
char clmLUHpasturedb[1024];
char clmLUHotherdb[1024];
char clmLUHc3anndb[1024];
char clmLUHc4anndb[1024];
char clmLUHc3perdb[1024];
char clmLUHc4perdb[1024];
char clmLUHc3nfxdb[1024];
char refstatesdb[1024];
char refmanagementdb[1024];
char reftransitionsdb[1024];
char luhstatesdb[1024];
char luhmanagementdb[1024];
char luhtransitionsdb[1024];
char pftparamfile[1024];
char cftrawparamfile[1024];
char cftparamfile[1024];
int flipLUHgrids;
int includeOcean;

char PFTluhtype[MAXPFT][256];
char CFTRAWluhtype[MAXCFTRAW][256];
char CFTluhtype[MAXCFT][256];

float *tempGrid;
float *tempflipGrid;
float *tempoutGrid;
float *translossGrid;

int *innatpft;
int *incft;
float inEDGEN;
float inEDGEE;
float inEDGES;
float inEDGEW;
float *inLAT;
float *inLATIXY;
float *inLON;
float *inLONGXY;

float *inLANDMASKGrid;
float *inLANDFRACGrid;
float *inAREAGrid;
float *inPCTGLACIERGrid;
float *inPCTLAKEGrid;
float *inPCTWETLANDGrid;
float *inPCTURBANGrid;
float *inPCTNATVEGGrid;
float *inPCTCROPGrid;
float *inCURRENTPCTPFTGrid[MAXPFT];
float *inCURRENTPCTCFTGrid[MAXCFT];

float *inFORESTPCTPFTGrid[MAXPFT];
float *inPASTUREPCTPFTGrid[MAXPFT];
float *inOTHERPCTPFTGrid[MAXPFT];

float *inC3ANNPCTCFTGrid[MAXCFTRAW];
float *inC4ANNPCTCFTGrid[MAXCFTRAW];
float *inC3PERPCTCFTGrid[MAXCFTRAW];
float *inC4PERPCTCFTGrid[MAXCFTRAW];
float *inC3NFXPCTCFTGrid[MAXCFTRAW];

float *inBASEPRIMFGrid;
float *inBASEPRIMNGrid;
float *inBASESECDFGrid;
float *inBASESECDNGrid;
float *inBASEPASTRGrid;
float *inBASERANGEGrid;
float *inBASEC3ANNGrid;
float *inBASEC4ANNGrid;
float *inBASEC3PERGrid;
float *inBASEC4PERGrid;
float *inBASEC3NFXGrid;
float *inBASEURBANGrid;

float *inCURRPRIMFGrid;
float *inCURRPRIMNGrid;
float *inCURRSECDFGrid;
float *inCURRSECDNGrid;
float *inCURRPASTRGrid;
float *inCURRRANGEGrid;
float *inCURRC3ANNGrid;
float *inCURRC4ANNGrid;
float *inCURRC3PERGrid;
float *inCURRC4PERGrid;
float *inCURRC3NFXGrid;
float *inCURRURBANGrid;

float *inPREVSECDFGrid;
float *inPREVSECDNGrid;
float *inPREVPASTRGrid;
float *inPREVRANGEGrid;
float *inPREVC3ANNGrid;
float *inPREVC4ANNGrid;
float *inPREVC3PERGrid;
float *inPREVC4PERGrid;
float *inPREVC3NFXGrid;

float *inHARVESTVH1Grid;
float *inHARVESTVH2Grid;
float *inHARVESTSH1Grid;
float *inHARVESTSH2Grid;
float *inHARVESTSH3Grid;

float *inBIOHVH1Grid;
float *inBIOHVH2Grid;
float *inBIOHSH1Grid;
float *inBIOHSH2Grid;
float *inBIOHSH3Grid;

float *inUNREPSECDFGrid;
float *inUNREPSECDNGrid;
float *inUNREPPASTRGrid;
float *inUNREPRANGEGrid;
float *inUNREPC3ANNGrid;
float *inUNREPC4ANNGrid;
float *inUNREPC3PERGrid;
float *inUNREPC4PERGrid;
float *inUNREPC3NFXGrid;

float *inFERTC3ANNGrid;
float *inFERTC4ANNGrid;
float *inFERTC3PERGrid;
float *inFERTC4PERGrid;
float *inFERTC3NFXGrid;

float *inIRRIGC3ANNGrid;
float *inIRRIGC4ANNGrid;
float *inIRRIGC3PERGrid;
float *inIRRIGC4PERGrid;
float *inIRRIGC3NFXGrid;

float *inBASEFORESTTOTALGrid;
float *inBASENONFORESTTOTALGrid;
float *inBASECROPTOTALGrid;
float *inBASEMISSINGGrid;
float *inBASEOTHERGrid;
float *inBASENATVEGGrid;

float *inCURRFORESTTOTALGrid;
float *inCURRNONFORESTTOTALGrid;
float *inCURRCROPTOTALGrid;
float *inCURRMISSINGGrid;
float *inCURROTHERGrid;
float *inCURRNATVEGGrid;

float *inUNREPFORESTGrid;
float *inUNREPOTHERGrid;

float *outPCTNATVEGGrid;
float *outPCTCROPGrid;
float *outPCTPFTGrid[MAXPFT];
float *outPCTCFTGrid[MAXCFT];

float *outFERTNITROGrid[MAXCFT];

float *outUNREPPFTGrid[MAXPFT];
float *outUNREPCFTGrid[MAXCFT];

float *outHARVESTVH1Grid;
float *outHARVESTVH2Grid;
float *outHARVESTSH1Grid;
float *outHARVESTSH2Grid;
float *outHARVESTSH3Grid;

float *outBIOHVH1Grid;
float *outBIOHVH2Grid;
float *outBIOHSH1Grid;
float *outBIOHSH2Grid;
float *outBIOHSH3Grid;

double *outLANDFRACdblGrid;
double *outAREAdblGrid;
double *outPCTGLACIERdblGrid;
double *outPCTLAKEdblGrid;
double *outPCTWETLANDdblGrid;
double *outPCTURBANdblGrid;

double *outPCTNATVEGdblGrid;
double *outPCTCROPdblGrid;
double *outPCTPFTdblGrid[MAXPFT];
double *outPCTCFTdblGrid[MAXCFT];

double *outFERTNITROdblGrid[MAXCFT];

double *outUNREPPFTdblGrid[MAXPFT];
double *outUNREPCFTdblGrid[MAXCFT];

double *outBIOHVH1dblGrid;
double *outBIOHVH2dblGrid;
double *outBIOHSH1dblGrid;
double *outBIOHSH2dblGrid;
double *outBIOHSH3dblGrid;

/* Out Surface Data NetCDF variables */
int  stat;  /* return status */
int  ncid;  /* netCDF id */

/* dimension ids */
int natpft_dim;
int cft_dim;
int lon_dim;
int lat_dim;
int nchar_dim;

/* dimension lengths */
size_t natpft_len = 15;
size_t cft_len = 64;
size_t lon_len = MAXCLMPIX;
size_t lat_len = MAXCLMLIN;
size_t nchar_len = 128;

/* variable ids */
int natpft_id;
int cft_id;
int EDGEN_id;
int EDGEE_id;
int EDGES_id;
int EDGEW_id;
int LAT_id;
int LATIXY_id;
int LON_id;
int LONGXY_id;
int LANDMASK_id;
int LANDFRAC_id;
int AREA_id;
int PCT_GLACIER_id;
int PCT_LAKE_id;
int PCT_WETLAND_id;
int PCT_URBAN_id;
int PCT_NATVEG_id;
int PCT_CROP_id;
int PCT_NAT_PFT_id;
int PCT_CFT_id;
int FERTNITRO_CFT_id;
int HARVEST_VH1_id;
int HARVEST_VH2_id;
int HARVEST_SH1_id;
int HARVEST_SH2_id;
int HARVEST_SH3_id;
int GRAZING_id;
int UNREPRESENTED_PFT_LULCC_id;
int UNREPRESENTED_CFT_LULCC_id;

/* variable shapes */
int natpft_dims[RANK_natpft];
int cft_dims[RANK_cft];
int LAT_dims[RANK_LAT];
int LATIXY_dims[RANK_LATIXY];
int LON_dims[RANK_LON];
int LONGXY_dims[RANK_LONGXY];
int LANDMASK_dims[RANK_LANDMASK];
int LANDFRAC_dims[RANK_LANDFRAC];
int AREA_dims[RANK_AREA];
int PCT_GLACIER_dims[RANK_PCT_GLACIER];
int PCT_LAKE_dims[RANK_PCT_LAKE];
int PCT_WETLAND_dims[RANK_PCT_WETLAND];
int PCT_URBAN_dims[RANK_PCT_URBAN];
int PCT_NATVEG_dims[RANK_PCT_NATVEG];
int PCT_CROP_dims[RANK_PCT_CROP];
int PCT_NAT_PFT_dims[RANK_PCT_NAT_PFT];
int PCT_CFT_dims[RANK_PCT_CFT];
int FERTNITRO_CFT_dims[RANK_FERTNITRO_CFT];
int HARVEST_VH1_dims[RANK_HARVEST_VH1];
int HARVEST_VH2_dims[RANK_HARVEST_VH2];
int HARVEST_SH1_dims[RANK_HARVEST_SH1];
int HARVEST_SH2_dims[RANK_HARVEST_SH2];
int HARVEST_SH3_dims[RANK_HARVEST_SH3];
int GRAZING_dims[RANK_GRAZING];
int UNREPRESENTED_PFT_LULCC_dims[RANK_UNREPRESENTED_PFT_LULCC];
int UNREPRESENTED_CFT_LULCC_dims[RANK_UNREPRESENTED_CFT_LULCC];

int readnamelist(char *namelist) {

  FILE *namelistfile;
  char fieldname[256];

  printf("Reading Namelist: %s\n",namelist);
  namelistfile = fopen(namelist,"r");
  
  fscanf(namelistfile,"%s %s",fieldname,regionfilename);
  fscanf(namelistfile,"%s %s",fieldname,outputdir);
  fscanf(namelistfile,"%s %s",fieldname,outputseries);
  fscanf(namelistfile,"%s %d",fieldname,&firstrefyear);
  fscanf(namelistfile,"%s %d",fieldname,&refyear);
  fscanf(namelistfile,"%s %d",fieldname,&firstyear);
  fscanf(namelistfile,"%s %d",fieldname,&startyear);
  fscanf(namelistfile,"%s %d",fieldname,&endyear);
  fscanf(namelistfile,"%s %s",fieldname,clmcurrentsurfdb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHforestdb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHpasturedb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHotherdb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHc3anndb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHc4anndb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHc3perdb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHc4perdb);
  fscanf(namelistfile,"%s %s",fieldname,clmLUHc3nfxdb);
  fscanf(namelistfile,"%s %s",fieldname,refstatesdb);
  fscanf(namelistfile,"%s %s",fieldname,luhstatesdb);
  fscanf(namelistfile,"%s %s",fieldname,luhmanagementdb);
  fscanf(namelistfile,"%s %s",fieldname,luhtransitionsdb);
  fscanf(namelistfile,"%s %s",fieldname,pftparamfile);
  fscanf(namelistfile,"%s %s",fieldname,cftrawparamfile);
  fscanf(namelistfile,"%s %s",fieldname,cftparamfile);  
  fscanf(namelistfile,"%s %d",fieldname,&flipLUHgrids);
  fscanf(namelistfile,"%s %d",fieldname,&includeOcean);

  return 0;

}

int setregionoptions() {

  FILE *pftparamfile;
  float lllon, lllat, urlon, urlat, pixsize;
  
  printf("Reading %s\n",regionfilename);
  pftparamfile = fopen(regionfilename,"r");

  fscanf(pftparamfile,"%f",&lllon);
  fscanf(pftparamfile,"%f",&lllat);
  fscanf(pftparamfile,"%f",&urlon);
  fscanf(pftparamfile,"%f",&urlat);
  fscanf(pftparamfile,"%f",&pixsize);
  
  OUTPIXSIZE = pixsize;
  MAXOUTPIX = (long) (urlon - lllon) / OUTPIXSIZE;
  MAXOUTLIN = (long) (urlat - lllat) / OUTPIXSIZE;

  lon_len = MAXOUTPIX;
  lat_len = MAXOUTLIN;
  
  OUTLLX = lllon;
  OUTLLY = lllat;
  
  OUTLATOFFSET = (long) (90.0 - urlat) / OUTPIXSIZE;
  OUTLONOFFSET = (long) (lllon + 180.0) / OUTPIXSIZE;

  OUTDATASIZE = MAXOUTPIX * MAXOUTLIN * sizeof(float);
  OUTDBLDATASIZE = MAXOUTPIX * MAXOUTLIN * sizeof(double);
  
  return 0;

}

int createallgrids() {

  int pftid, cftid;

  tempGrid = (float *) malloc(OUTDATASIZE);
  tempoutGrid = (float *) malloc(OUTDATASIZE);
  tempflipGrid = (float *) malloc(OUTDATASIZE);
  translossGrid = (float *) malloc(OUTDATASIZE);

  innatpft = (int *) malloc(MAXPFT * sizeof(int));
  incft = (int *) malloc(MAXCFT * sizeof(int));
  inLAT = (float *) malloc(MAXOUTLIN * sizeof(float));
  inLATIXY = (float *) malloc(OUTDATASIZE);
  inLON = (float *) malloc(MAXOUTPIX * sizeof(float));
  inLONGXY = (float *) malloc(OUTDATASIZE);

  inLANDMASKGrid = (float *) malloc(OUTDATASIZE);
  inLANDFRACGrid = (float *) malloc(OUTDATASIZE);
  inAREAGrid = (float *) malloc(OUTDATASIZE);
  inPCTGLACIERGrid = (float *) malloc(OUTDATASIZE);
  inPCTLAKEGrid = (float *) malloc(OUTDATASIZE);
  inPCTWETLANDGrid = (float *) malloc(OUTDATASIZE);
  inPCTURBANGrid = (float *) malloc(OUTDATASIZE);
  inPCTNATVEGGrid = (float *) malloc(OUTDATASIZE);
  inPCTCROPGrid = (float *) malloc(OUTDATASIZE);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      inCURRENTPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {  
      inCURRENTPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      inFORESTPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
      inPASTUREPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
      inOTHERPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      inC3ANNPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC4ANNPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC3PERPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC4PERPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
      inC3NFXPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }

  inBASEPRIMFGrid = (float *) malloc(OUTDATASIZE);
  inBASEPRIMNGrid = (float *) malloc(OUTDATASIZE);
  inBASESECDFGrid = (float *) malloc(OUTDATASIZE);
  inBASESECDNGrid = (float *) malloc(OUTDATASIZE);
  inBASEPASTRGrid = (float *) malloc(OUTDATASIZE);
  inBASERANGEGrid = (float *) malloc(OUTDATASIZE);
  inBASEC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inBASEC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inBASEC3PERGrid = (float *) malloc(OUTDATASIZE);
  inBASEC4PERGrid = (float *) malloc(OUTDATASIZE);
  inBASEC3NFXGrid = (float *) malloc(OUTDATASIZE);
  inBASEURBANGrid = (float *) malloc(OUTDATASIZE);

  inCURRPRIMFGrid = (float *) malloc(OUTDATASIZE);
  inCURRPRIMNGrid = (float *) malloc(OUTDATASIZE);
  inCURRSECDFGrid = (float *) malloc(OUTDATASIZE);
  inCURRSECDNGrid = (float *) malloc(OUTDATASIZE);
  inCURRPASTRGrid = (float *) malloc(OUTDATASIZE);
  inCURRRANGEGrid = (float *) malloc(OUTDATASIZE);
  inCURRC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inCURRC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inCURRC3PERGrid = (float *) malloc(OUTDATASIZE);
  inCURRC4PERGrid = (float *) malloc(OUTDATASIZE);
  inCURRC3NFXGrid = (float *) malloc(OUTDATASIZE);
  inCURRURBANGrid = (float *) malloc(OUTDATASIZE);

  inPREVSECDFGrid = (float *) malloc(OUTDATASIZE);
  inPREVSECDNGrid = (float *) malloc(OUTDATASIZE);
  inPREVPASTRGrid = (float *) malloc(OUTDATASIZE);
  inPREVRANGEGrid = (float *) malloc(OUTDATASIZE);
  inPREVC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inPREVC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inPREVC3PERGrid = (float *) malloc(OUTDATASIZE);
  inPREVC4PERGrid = (float *) malloc(OUTDATASIZE);
  inPREVC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inHARVESTVH1Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTVH2Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTSH1Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTSH2Grid = (float *) malloc(OUTDATASIZE);
  inHARVESTSH3Grid = (float *) malloc(OUTDATASIZE);

  inBIOHVH1Grid = (float *) malloc(OUTDATASIZE);
  inBIOHVH2Grid = (float *) malloc(OUTDATASIZE);
  inBIOHSH1Grid = (float *) malloc(OUTDATASIZE);
  inBIOHSH2Grid = (float *) malloc(OUTDATASIZE);
  inBIOHSH3Grid = (float *) malloc(OUTDATASIZE);

  inUNREPSECDFGrid = (float *) malloc(OUTDATASIZE);
  inUNREPSECDNGrid = (float *) malloc(OUTDATASIZE);
  inUNREPPASTRGrid = (float *) malloc(OUTDATASIZE);
  inUNREPRANGEGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC3PERGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC4PERGrid = (float *) malloc(OUTDATASIZE);
  inUNREPC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inFERTC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inFERTC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inFERTC3PERGrid = (float *) malloc(OUTDATASIZE);
  inFERTC4PERGrid = (float *) malloc(OUTDATASIZE);
  inFERTC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inIRRIGC3ANNGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC4ANNGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC3PERGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC4PERGrid = (float *) malloc(OUTDATASIZE);
  inIRRIGC3NFXGrid = (float *) malloc(OUTDATASIZE);

  inBASEFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inBASENONFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inBASECROPTOTALGrid = (float *) malloc(OUTDATASIZE);
  inBASEMISSINGGrid = (float *) malloc(OUTDATASIZE);
  inBASEOTHERGrid = (float *) malloc(OUTDATASIZE);
  inBASENATVEGGrid = (float *) malloc(OUTDATASIZE);

  inCURRFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inCURRNONFORESTTOTALGrid = (float *) malloc(OUTDATASIZE);
  inCURRCROPTOTALGrid = (float *) malloc(OUTDATASIZE);
  inCURRMISSINGGrid = (float *) malloc(OUTDATASIZE);
  inCURROTHERGrid = (float *) malloc(OUTDATASIZE);
  inCURRNATVEGGrid = (float *) malloc(OUTDATASIZE);

  inUNREPFORESTGrid = (float *) malloc(OUTDATASIZE);
  inUNREPOTHERGrid = (float *) malloc(OUTDATASIZE);

  outPCTNATVEGGrid = (float *) malloc(OUTDATASIZE);
  outPCTCROPGrid = (float *) malloc(OUTDATASIZE);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outPCTPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outPCTCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outUNREPPFTGrid[pftid] = (float *) malloc(OUTDATASIZE);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outUNREPCFTGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outFERTNITROGrid[cftid] = (float *) malloc(OUTDATASIZE);
  }
  
  outHARVESTVH1Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTVH2Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTSH1Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTSH2Grid = (float *) malloc(OUTDATASIZE);
  outHARVESTSH3Grid = (float *) malloc(OUTDATASIZE);

  outBIOHVH1Grid = (float *) malloc(OUTDATASIZE);
  outBIOHVH2Grid = (float *) malloc(OUTDATASIZE);
  outBIOHSH1Grid = (float *) malloc(OUTDATASIZE);
  outBIOHSH2Grid = (float *) malloc(OUTDATASIZE);
  outBIOHSH3Grid = (float *) malloc(OUTDATASIZE);

  outLANDFRACdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outAREAdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTGLACIERdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTLAKEdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTWETLANDdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTURBANdblGrid = (double *) malloc(OUTDBLDATASIZE);

  outPCTNATVEGdblGrid = (double *) malloc(OUTDBLDATASIZE);
  outPCTCROPdblGrid = (double *) malloc(OUTDBLDATASIZE);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outPCTPFTdblGrid[pftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outPCTCFTdblGrid[cftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      outUNREPPFTdblGrid[pftid] = (double *) malloc(OUTDBLDATASIZE);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outUNREPCFTdblGrid[cftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      outFERTNITROdblGrid[cftid] = (double *) malloc(OUTDBLDATASIZE);
  }
  
  outBIOHVH1dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outBIOHVH2dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outBIOHSH1dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outBIOHSH2dblGrid = (double *) malloc(OUTDBLDATASIZE);
  outBIOHSH3dblGrid = (double *) malloc(OUTDBLDATASIZE);

  return 0;

}


int
readpftparamfile() {

  FILE *pftparaminfile;
  int inpft, inpftid;
  char inPFTluhtype[256], inPFTname[256];
  
  printf("Reading %s\n",pftparamfile);
  pftparaminfile = fopen(pftparamfile,"r");

  for (inpft = 0; inpft < MAXPFT; inpft++) {
      fscanf(pftparaminfile,"%d%s%s",&inpftid,inPFTluhtype,inPFTname);
      sprintf(PFTluhtype[inpft],"%s",inPFTluhtype);
  }  
  
  return 0;

}

int
readcftrawparamfile() {

  FILE *cftrawparaminfile;
  int incftraw, incftrawid;
  char inCFTRAWluhtype[256], inCFTRAWname[256];
  
  printf("Reading %s\n",cftrawparamfile);
  cftrawparaminfile = fopen(cftrawparamfile,"r");

  for (incftraw = 0; incftraw < MAXCFTRAW; incftraw++) {
      fscanf(cftrawparaminfile,"%d%s%s",&incftrawid,inCFTRAWluhtype,inCFTRAWname);
      sprintf(CFTRAWluhtype[incftraw],"%s",inCFTRAWluhtype);
  }  
  
  return 0;

}

int
readcftparamfile() {

  FILE *cftparaminfile;
  int incft, incftid;
  char inCFTluhtype[256], inCFTname[256];
  
  printf("Reading %s\n",cftparamfile);
  cftparaminfile = fopen(cftparamfile,"r");

  for (incft = 0; incft < MAXCFT; incft++) {
      fscanf(cftparaminfile,"%d%s%s",&incftid,inCFTluhtype,inCFTname);
      sprintf(CFTluhtype[incft],"%s",inCFTluhtype);
  }  
  
  return 0;

}


int initializeGrids() {

  long clmlin, clmpix;
  long pftid, cftid;
  float clmlat, clmlatdistance, clmlondistance, clmarea;
  
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          outPCTNATVEGGrid[MAXOUTPIX * MAXOUTLIN] = 0.0;
          outPCTCROPGrid[MAXOUTPIX * MAXOUTLIN] = 0.0;
          for (pftid = 0; pftid < MAXPFT; pftid++) {
              outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              outUNREPPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
          }

          outHARVESTVH1Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outHARVESTVH2Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outHARVESTSH1Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outHARVESTSH2Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outHARVESTSH3Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;

          outBIOHVH1Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outBIOHVH2Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outBIOHSH1Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outBIOHSH2Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          outBIOHSH3Grid[clmlin * MAXOUTPIX + clmpix] = 0.0;

          for (cftid = 0; cftid < MAXCFT; cftid++) {
              outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              outFERTNITROGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              outUNREPCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
          }
      }
  }

  return 0;
                
}

void
check_err(const int stat, const int line, const char *file) {

    if (stat != NC_NOERR) {
        (void)fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
        fflush(stderr);
        exit(1);
    }
}

int
openncinputfile(char *netcdffilename) {

    printf("Opening NetCDF File: %s\n",netcdffilename); 
    stat = nc_open(netcdffilename, NC_NOWRITE, &ncid);
    check_err(stat,__LINE__,__FILE__);

    return 0;

}

int
openncoutputfile(char *netcdffilename) {

    printf("Opening NetCDF File: %s\n",netcdffilename); 
    stat = nc_open(netcdffilename, NC_WRITE, &ncid);
    check_err(stat,__LINE__,__FILE__);

    return 0;

}

int
createncoutputfile(char *netcdffilename) {

    printf("Creating NetCDF File: %s\n",netcdffilename); 

    /* enter define mode */
    stat = nc_create(netcdffilename, NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL, &ncid);
    check_err(stat,__LINE__,__FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "natpft", natpft_len, &natpft_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "cft", cft_len, &cft_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "lon", lon_len, &lon_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "lat", lat_len, &lat_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "nchar", nchar_len, &nchar_dim);
    check_err(stat,__LINE__,__FILE__);

    /* define variables */

    natpft_dims[0] = natpft_dim;
    stat = nc_def_var(ncid, "natpft", NC_INT, RANK_natpft, natpft_dims, &natpft_id);
    check_err(stat,__LINE__,__FILE__);

    cft_dims[0] = cft_dim;
    stat = nc_def_var(ncid, "cft", NC_INT, RANK_cft, cft_dims, &cft_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGEN", NC_FLOAT, RANK_EDGEN, 0, &EDGEN_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGEE", NC_FLOAT, RANK_EDGEE, 0, &EDGEE_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGES", NC_FLOAT, RANK_EDGES, 0, &EDGES_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "EDGEW", NC_FLOAT, RANK_EDGEW, 0, &EDGEW_id);
    check_err(stat,__LINE__,__FILE__);

    LAT_dims[0] = lat_dim;
    stat = nc_def_var(ncid, "LAT", NC_FLOAT, RANK_LAT, LAT_dims, &LAT_id);
    check_err(stat,__LINE__,__FILE__);

    LATIXY_dims[0] = lat_dim;
    LATIXY_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LATIXY", NC_FLOAT, RANK_LATIXY, LATIXY_dims, &LATIXY_id);
    check_err(stat,__LINE__,__FILE__);

    LON_dims[0] = lon_dim;
    stat = nc_def_var(ncid, "LON", NC_FLOAT, RANK_LON, LON_dims, &LON_id);
    check_err(stat,__LINE__,__FILE__);

    LONGXY_dims[0] = lat_dim;
    LONGXY_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LONGXY", NC_FLOAT, RANK_LONGXY, LONGXY_dims, &LONGXY_id);
    check_err(stat,__LINE__,__FILE__);

    LANDMASK_dims[0] = lat_dim;
    LANDMASK_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LANDMASK", NC_FLOAT, RANK_LANDMASK, LANDMASK_dims, &LANDMASK_id);
    check_err(stat,__LINE__,__FILE__);

    LANDFRAC_dims[0] = lat_dim;
    LANDFRAC_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "LANDFRAC", NC_DOUBLE, RANK_LANDFRAC, LANDFRAC_dims, &LANDFRAC_id);
    check_err(stat,__LINE__,__FILE__);

    AREA_dims[0] = lat_dim;
    AREA_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "AREA", NC_DOUBLE, RANK_AREA, AREA_dims, &AREA_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_GLACIER_dims[0] = lat_dim;
    PCT_GLACIER_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_GLACIER", NC_DOUBLE, RANK_PCT_GLACIER, PCT_GLACIER_dims, &PCT_GLACIER_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_LAKE_dims[0] = lat_dim;
    PCT_LAKE_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_LAKE", NC_DOUBLE, RANK_PCT_LAKE, PCT_LAKE_dims, &PCT_LAKE_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_WETLAND_dims[0] = lat_dim;
    PCT_WETLAND_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_WETLAND", NC_DOUBLE, RANK_PCT_WETLAND, PCT_WETLAND_dims, &PCT_WETLAND_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_URBAN_dims[0] = lat_dim;
    PCT_URBAN_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_URBAN", NC_DOUBLE, RANK_PCT_URBAN, PCT_URBAN_dims, &PCT_URBAN_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_NATVEG_dims[0] = lat_dim;
    PCT_NATVEG_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_NATVEG", NC_DOUBLE, RANK_PCT_NATVEG, PCT_NATVEG_dims, &PCT_NATVEG_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_CROP_dims[0] = lat_dim;
    PCT_CROP_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "PCT_CROP", NC_DOUBLE, RANK_PCT_CROP, PCT_CROP_dims, &PCT_CROP_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_NAT_PFT_dims[0] = natpft_dim;
    PCT_NAT_PFT_dims[1] = lat_dim;
    PCT_NAT_PFT_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "PCT_NAT_PFT", NC_DOUBLE, RANK_PCT_NAT_PFT, PCT_NAT_PFT_dims, &PCT_NAT_PFT_id);
    check_err(stat,__LINE__,__FILE__);

    PCT_CFT_dims[0] = cft_dim;
    PCT_CFT_dims[1] = lat_dim;
    PCT_CFT_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "PCT_CFT", NC_DOUBLE, RANK_PCT_CFT, PCT_CFT_dims, &PCT_CFT_id);
    check_err(stat,__LINE__,__FILE__);

    FERTNITRO_CFT_dims[0] = cft_dim;
    FERTNITRO_CFT_dims[1] = lat_dim;
    FERTNITRO_CFT_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "FERTNITRO_CFT", NC_DOUBLE, RANK_FERTNITRO_CFT, FERTNITRO_CFT_dims, &FERTNITRO_CFT_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_VH1_dims[0] = lat_dim;
    HARVEST_VH1_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_VH1", NC_DOUBLE, RANK_HARVEST_VH1, HARVEST_VH1_dims, &HARVEST_VH1_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_VH2_dims[0] = lat_dim;
    HARVEST_VH2_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_VH2", NC_DOUBLE, RANK_HARVEST_VH2, HARVEST_VH2_dims, &HARVEST_VH2_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_SH1_dims[0] = lat_dim;
    HARVEST_SH1_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_SH1", NC_DOUBLE, RANK_HARVEST_SH1, HARVEST_SH1_dims, &HARVEST_SH1_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_SH2_dims[0] = lat_dim;
    HARVEST_SH2_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_SH2", NC_DOUBLE, RANK_HARVEST_SH2, HARVEST_SH2_dims, &HARVEST_SH2_id);
    check_err(stat,__LINE__,__FILE__);

    HARVEST_SH3_dims[0] = lat_dim;
    HARVEST_SH3_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "HARVEST_SH3", NC_DOUBLE, RANK_HARVEST_SH3, HARVEST_SH3_dims, &HARVEST_SH3_id);
    check_err(stat,__LINE__,__FILE__);

    GRAZING_dims[0] = lat_dim;
    GRAZING_dims[1] = lon_dim;
    stat = nc_def_var(ncid, "GRAZING", NC_DOUBLE, RANK_GRAZING, GRAZING_dims, &GRAZING_id);
    check_err(stat,__LINE__,__FILE__);

    UNREPRESENTED_PFT_LULCC_dims[0] = natpft_dim;
    UNREPRESENTED_PFT_LULCC_dims[1] = lat_dim;
    UNREPRESENTED_PFT_LULCC_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "UNREPRESENTED_PFT_LULCC", NC_DOUBLE, RANK_UNREPRESENTED_PFT_LULCC, UNREPRESENTED_PFT_LULCC_dims, &UNREPRESENTED_PFT_LULCC_id);
    check_err(stat,__LINE__,__FILE__);

    UNREPRESENTED_CFT_LULCC_dims[0] = cft_dim;
    UNREPRESENTED_CFT_LULCC_dims[1] = lat_dim;
    UNREPRESENTED_CFT_LULCC_dims[2] = lon_dim;
    stat = nc_def_var(ncid, "UNREPRESENTED_CFT_LULCC", NC_DOUBLE, RANK_UNREPRESENTED_CFT_LULCC, UNREPRESENTED_CFT_LULCC_dims, &UNREPRESENTED_CFT_LULCC_id);
    check_err(stat,__LINE__,__FILE__);

    /* assign global attributes */

    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "source", 20, "Peter Lawrence, NCAR");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "creation_date", 28, "Tue Jun 13 16:42:45 MDT 2017");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "title", 18, "mksrf_file.nc");
    check_err(stat,__LINE__,__FILE__);
    }


    /* assign per-variable attributes */

    {
    stat = nc_put_att_text(ncid, natpft_id, "long_name", 23, "indices of natural PFTs");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, natpft_id, "units", 5, "index");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, cft_id, "long_name", 15, "indices of CFTs");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, cft_id, "units", 5, "index");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEN_id, "long_name", 29, "northern edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEN_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEE_id, "long_name", 28, "eastern edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEE_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGES_id, "long_name", 29, "southern edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGES_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEW_id, "long_name", 28, "western edge of surface grid");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, EDGEW_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LAT_id, "long_name", 3, "lat");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LAT_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const float mksrf_file__FillValue_att[1] = {((float)9.96921e+36)} ;
    stat = nc_put_att_float(ncid, LATIXY_id, "_FillValue", NC_FLOAT, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LATIXY_id, "long_name", 11, "latitude-2d");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LATIXY_id, "units", 13, "degrees north");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LON_id, "long_name", 3, "lon");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LON_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const float mksrf_file__FillValue_att[1] = {((float)9.96921e+36)} ;
    stat = nc_put_att_float(ncid, LONGXY_id, "_FillValue", NC_FLOAT, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LONGXY_id, "long_name", 12, "longitude-2d");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LONGXY_id, "units", 12, "degrees east");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDMASK_id, "long_name", 9, "land mask");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDMASK_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDFRAC_id, "long_name", 25, "land fraction of gridcell");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, LANDFRAC_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, AREA_id, "long_name", 16, "area of gridcell");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, AREA_id, "units", 4, "km^2");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_GLACIER_id, "long_name", 30, "total percent glacier landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_GLACIER_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_GLACIER_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_LAKE_id, "long_name", 27, "total percent lake landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_LAKE_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_LAKE_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_WETLAND_id, "long_name", 30, "total percent wetland landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_WETLAND_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_WETLAND_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_URBAN_id, "long_name", 28, "total percent urban landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_URBAN_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_URBAN_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NATVEG_id, "long_name", 41, "total percent natural vegetation landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NATVEG_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_NATVEG_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CROP_id, "long_name", 27, "total percent crop landunit");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CROP_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_CROP_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NAT_PFT_id, "long_name", 73, "percent plant functional type on the natural veg landunit (% of landunit)");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_NAT_PFT_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_NAT_PFT_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CFT_id, "long_name", 65, "percent crop functional type on the crop landunit (% of landunit)");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, PCT_CFT_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, PCT_CFT_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, FERTNITRO_CFT_id, "long_name", 33, "nitrogen fertilizer for each crop");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, FERTNITRO_CFT_id, "units", 8, "gN/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, FERTNITRO_CFT_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH1_id, "long_name", 27, "harvest from primary forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH1_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_VH1_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH2_id, "long_name", 31, "harvest from primary non-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_VH2_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_VH2_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH1_id, "long_name", 36, "harvest from secondary mature-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH1_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_SH1_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH2_id, "long_name", 35, "harvest from secondary young-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH2_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_SH2_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH3_id, "long_name", 33, "harvest from secondary non-forest");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, HARVEST_SH3_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, HARVEST_SH3_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, GRAZING_id, "long_name", 25, "grazing of herbacous pfts");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, GRAZING_id, "units", 8, "gC/m2/yr");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, GRAZING_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_PFT_LULCC_id, "long_name", 41, "unrepresented PFT gross LULCC transitions");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_PFT_LULCC_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, UNREPRESENTED_PFT_LULCC_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_CFT_LULCC_id, "long_name", 42, "unrepresented crop gross LULCC transitions");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, UNREPRESENTED_CFT_LULCC_id, "units", 8, "unitless");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double mksrf_file__FillValue_att[1] = {((double)-9999)} ;
    stat = nc_put_att_double(ncid, UNREPRESENTED_CFT_LULCC_id, "_FillValue", NC_DOUBLE, 1, mksrf_file__FillValue_att);
    check_err(stat,__LINE__,__FILE__);
    }


    /* leave define mode */
    stat = nc_enddef (ncid);
    check_err(stat,__LINE__,__FILE__);

    /* assign variable data */

    stat = nc_close(ncid);
    check_err(stat,__LINE__,__FILE__);
    return 0;
}

int
closencfile() {

    stat = nc_close(ncid);
    check_err(stat,__LINE__,__FILE__);

    return 0;

}

int readnc0dfield(char *FieldName, float *targetvalue) {

    int varid;
        
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_get_var_float(ncid, varid, targetvalue);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int readnc1dfield(char *FieldName, float *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_get_var_float(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int readnc1dintfield(char *FieldName, int *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_get_var_int(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int readnc2dfield(char *FieldName, float *targetgrid, int flipgrid) {

    int varid;
    long clmlin, clmpix, fliplin;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    if (flipgrid == 0) {
        stat =  nc_get_var_float(ncid, varid, targetgrid);
        check_err(stat,__LINE__,__FILE__);
    }
    else {
        stat =  nc_get_var_float(ncid, varid, tempflipGrid);
        check_err(stat,__LINE__,__FILE__);
        for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
            fliplin = MAXOUTLIN - clmlin - 1;
            for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
                targetgrid[clmlin * MAXOUTPIX + clmpix] = tempflipGrid[fliplin * MAXOUTPIX + clmpix];
            }
	}
    }
    
    return 0;

}


int readnc3dfield(char *FieldName, int index3d, float *targetgrid, int flipgrid) {

    int varid;
    long clmlin, clmpix, fliplin;
    size_t start[3], count[3];
    
    count[0] = 1;
    count[1] = MAXOUTLIN;
    count[2] = MAXOUTPIX;
    start[0] = index3d;
    start[1] = 0;
    start[2] = 0;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    if (flipgrid == 0) {
        stat =  nc_get_vara_float(ncid, varid, start, count, targetgrid);
        check_err(stat,__LINE__,__FILE__);
    }
    else {
        stat =  nc_get_vara_float(ncid, varid, start, count, tempflipGrid);
        check_err(stat,__LINE__,__FILE__);
        for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
            fliplin = MAXOUTLIN - clmlin - 1;
            for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
                targetgrid[clmlin * MAXOUTPIX + clmpix] = tempflipGrid[fliplin * MAXOUTPIX + clmpix];
            }
	}
    }
        
    return 0;
    
}

int writenc0dfield(char *FieldName, float *targetvalue) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_float(ncid, varid, targetvalue);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc1dfield(char *FieldName, float *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_float(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc1dintfield(char *FieldName, int *targetarray) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_int(ncid, varid, targetarray);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc2dfield(char *FieldName, float *targetgrid) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_float(ncid, varid, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc3dfield(char *FieldName, int index3d, float *targetgrid) {

    int varid;
    size_t start[3], count[3];
    
    count[0] = 1;
    count[1] = MAXOUTLIN;
    count[2] = MAXOUTPIX;
    start[0] = index3d;
    start[1] = 0;
    start[2] = 0;
        
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_vara_float(ncid, varid, start, count, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;
    
}

int writenc2ddblfield(char *FieldName, double *targetgrid) {

    int varid;
    
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_var_double(ncid, varid, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;

}

int writenc3ddblfield(char *FieldName, int index3d, double *targetgrid) {

    int varid;
    size_t start[3], count[3];
    
    count[0] = 1;
    count[1] = MAXOUTLIN;
    count[2] = MAXOUTPIX;
    start[0] = index3d;
    start[1] = 0;
    start[2] = 0;
        
    stat =  nc_inq_varid(ncid, FieldName, &varid);
    check_err(stat,__LINE__,__FILE__);

    stat =  nc_put_vara_double(ncid, varid, start, count, targetgrid);
    check_err(stat,__LINE__,__FILE__);
    
    return 0;
    
}

int readclmcurrentGrids() {

  int pftid, cftid;

  openncinputfile(clmcurrentsurfdb);
  
  readnc1dintfield("natpft",innatpft);
  readnc1dintfield("cft",incft);
  readnc0dfield("EDGEN",&inEDGEN);
  readnc0dfield("EDGEE",&inEDGEE);
  readnc0dfield("EDGES",&inEDGES);
  readnc0dfield("EDGEW",&inEDGEW);
  readnc1dfield("LAT",inLAT);
  readnc2dfield("LATIXY",inLATIXY,0);
  readnc1dfield("LON",inLON);
  readnc2dfield("LONGXY",inLONGXY,0);
  readnc2dfield("LANDMASK",inLANDMASKGrid,0);
  readnc2dfield("LANDFRAC",inLANDFRACGrid,0);
  readnc2dfield("AREA",inAREAGrid,0);
  readnc2dfield("PCT_GLACIER",inPCTGLACIERGrid,0);
  readnc2dfield("PCT_LAKE",inPCTLAKEGrid,0);
  readnc2dfield("PCT_WETLAND",inPCTWETLANDGrid,0);
  readnc2dfield("PCT_URBAN",inPCTURBANGrid,0);
  readnc2dfield("PCT_NATVEG",inPCTNATVEGGrid,0);
  readnc2dfield("PCT_CROP",inPCTCROPGrid,0);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc3dfield("PCT_NAT_PFT",pftid,inCURRENTPCTPFTGrid[pftid],0);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      readnc3dfield("PCT_CFT",cftid,inCURRENTPCTCFTGrid[cftid],0);
  }

  closencfile();

  return 0;
  
}


int readclmLUHforestGrids() {

  int pftid, cftid;

  openncinputfile(clmLUHforestdb);  

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc3dfield("PCT_NAT_PFT",pftid,inFORESTPCTPFTGrid[pftid],0);
  }
  
  closencfile();
  
  return 0;
  
}
  
  
int readclmLUHpastureGrids() {

  int pftid, cftid;

  openncinputfile(clmLUHpasturedb);  

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc3dfield("PCT_NAT_PFT",pftid,inPASTUREPCTPFTGrid[pftid],0);
  }
  
  closencfile();
  
  return 0;
  
}


int readclmLUHotherGrids() {

  int pftid, cftid;

  openncinputfile(clmLUHotherdb);  

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      readnc3dfield("PCT_NAT_PFT",pftid,inOTHERPCTPFTGrid[pftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readclmLUHc3annGrids() {

  int pftid, cftid;

  openncinputfile(clmLUHc3anndb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc3dfield("PCT_CFT",cftid,inC3ANNPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readclmLUHc4annGrids() {

  int pftid, cftid;

  openncinputfile(clmLUHc4anndb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc3dfield("PCT_CFT",cftid,inC4ANNPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readclmLUHc3perGrids() {

  int pftid, cftid;

  openncinputfile(clmLUHc3perdb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc3dfield("PCT_CFT",cftid,inC3PERPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readclmLUHc4perGrids() {

  int pftid, cftid;
  long clmlin, clmpix;

  openncinputfile(clmLUHc4perdb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc3dfield("PCT_CFT",cftid,inC4PERPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}

int readclmLUHc3nfxGrids() {

  int pftid, cftid;
  long clmlin, clmpix;

  openncinputfile(clmLUHc3nfxdb);  

  for (cftid = 0; cftid < MAXCFTRAW; cftid++) {
      readnc3dfield("PCT_CFT",cftid,inC3NFXPCTCFTGrid[cftid],0);
  }
  
  closencfile();
  
  return 0;
  
}


int readLUHbasestateGrids() {

  int yearindex;
  
  yearindex = refyear - firstrefyear;
  
  openncinputfile(refstatesdb); 

  readnc3dfield("primf",yearindex,inBASEPRIMFGrid,flipLUHgrids);
  readnc3dfield("primn",yearindex,inBASEPRIMNGrid,flipLUHgrids);
  readnc3dfield("secdf",yearindex,inBASESECDFGrid,flipLUHgrids);
  readnc3dfield("secdn",yearindex,inBASESECDNGrid,flipLUHgrids);
  readnc3dfield("pastr",yearindex,inBASEPASTRGrid,flipLUHgrids);
  readnc3dfield("range",yearindex,inBASERANGEGrid,flipLUHgrids);
  readnc3dfield("c3ann",yearindex,inBASEC3ANNGrid,flipLUHgrids);
  readnc3dfield("c4ann",yearindex,inBASEC4ANNGrid,flipLUHgrids);
  readnc3dfield("c3per",yearindex,inBASEC3PERGrid,flipLUHgrids);
  readnc3dfield("c4per",yearindex,inBASEC4PERGrid,flipLUHgrids);
  readnc3dfield("c3nfx",yearindex,inBASEC3NFXGrid,flipLUHgrids);
  readnc3dfield("urban",yearindex,inBASEURBANGrid,flipLUHgrids);
  
  closencfile();

  return 0;
  
}


int readLUHcurrstateGrids(int currentyear) {

  int yearindex;
  
  yearindex = currentyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  openncinputfile(luhstatesdb); 

  readnc3dfield("primf",yearindex,inCURRPRIMFGrid,flipLUHgrids);
  readnc3dfield("primn",yearindex,inCURRPRIMNGrid,flipLUHgrids);
  readnc3dfield("secdf",yearindex,inCURRSECDFGrid,flipLUHgrids);
  readnc3dfield("secdn",yearindex,inCURRSECDNGrid,flipLUHgrids);
  readnc3dfield("pastr",yearindex,inCURRPASTRGrid,flipLUHgrids);
  readnc3dfield("range",yearindex,inCURRRANGEGrid,flipLUHgrids);
  readnc3dfield("c3ann",yearindex,inCURRC3ANNGrid,flipLUHgrids);
  readnc3dfield("c4ann",yearindex,inCURRC4ANNGrid,flipLUHgrids);
  readnc3dfield("c3per",yearindex,inCURRC3PERGrid,flipLUHgrids);
  readnc3dfield("c4per",yearindex,inCURRC4PERGrid,flipLUHgrids);
  readnc3dfield("c3nfx",yearindex,inCURRC3NFXGrid,flipLUHgrids);
  readnc3dfield("urban",yearindex,inCURRURBANGrid,flipLUHgrids);
  
  closencfile();

  return 0;
  
}

int readLUHprevstateGrids(int prevyear) {

  int yearindex;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhstatesdb); 

  readnc3dfield("secdf",yearindex,inPREVSECDFGrid,flipLUHgrids);
  readnc3dfield("secdn",yearindex,inPREVSECDNGrid,flipLUHgrids);
  readnc3dfield("pastr",yearindex,inPREVPASTRGrid,flipLUHgrids);
  readnc3dfield("range",yearindex,inPREVRANGEGrid,flipLUHgrids);
  readnc3dfield("c3ann",yearindex,inPREVC3ANNGrid,flipLUHgrids);
  readnc3dfield("c4ann",yearindex,inPREVC4ANNGrid,flipLUHgrids);
  readnc3dfield("c3per",yearindex,inPREVC3PERGrid,flipLUHgrids);
  readnc3dfield("c4per",yearindex,inPREVC4PERGrid,flipLUHgrids);
  readnc3dfield("c3nfx",yearindex,inPREVC3NFXGrid,flipLUHgrids);
  
  closencfile();

  return 0;
  
}


int readLUHwoodharvestGrids(int prevyear) {

  int yearindex;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("primf_harv",yearindex,inHARVESTVH1Grid,flipLUHgrids);
  readnc3dfield("primn_harv",yearindex,inHARVESTVH2Grid,flipLUHgrids);
  readnc3dfield("secmf_harv",yearindex,inHARVESTSH1Grid,flipLUHgrids);
  readnc3dfield("secyf_harv",yearindex,inHARVESTSH2Grid,flipLUHgrids);
  readnc3dfield("secnf_harv",yearindex,inHARVESTSH3Grid,flipLUHgrids);
  readnc3dfield("primf_bioh",yearindex,inBIOHVH1Grid,flipLUHgrids);
  readnc3dfield("primn_bioh",yearindex,inBIOHVH2Grid,flipLUHgrids);
  readnc3dfield("secmf_bioh",yearindex,inBIOHSH1Grid,flipLUHgrids);
  readnc3dfield("secyf_bioh",yearindex,inBIOHSH2Grid,flipLUHgrids);
  readnc3dfield("secnf_bioh",yearindex,inBIOHSH3Grid,flipLUHgrids);
  
  closencfile();
  
  return 0;
  
}


int readUNREPSECDFGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("secdf_to_secdn",yearindex,translossGrid,flipLUHgrids);
  
  readnc3dfield("secdf_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdf_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdf_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdf_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdf_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdf_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdf_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdf_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVSECDFGrid[clmlin * MAXOUTPIX + clmpix] - inCURRSECDFGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPSECDFGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }
  
  closencfile();

  return 0;
  
}


int readUNREPSECDNGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("secdn_to_secdf",yearindex,translossGrid,flipLUHgrids);

  readnc3dfield("secdn_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdn_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdn_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdn_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdn_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdn_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdn_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("secdn_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVSECDNGrid[clmlin * MAXOUTPIX + clmpix] - inCURRSECDNGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPSECDNGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readUNREPC3ANNGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("c3ann_to_secdn",yearindex,translossGrid,flipLUHgrids);
  
  readnc3dfield("c3ann_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3ann_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3ann_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3ann_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3ann_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3ann_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3ann_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3ann_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVC3ANNGrid[clmlin * MAXOUTPIX + clmpix] - inCURRC3ANNGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPC3ANNGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readUNREPC4ANNGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("c4ann_to_secdn",yearindex,translossGrid,flipLUHgrids);

  readnc3dfield("c4ann_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4ann_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4ann_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4ann_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4ann_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4ann_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4ann_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4ann_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVC4ANNGrid[clmlin * MAXOUTPIX + clmpix] - inCURRC4ANNGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPC4ANNGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readUNREPC3PERGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("c3per_to_secdn",yearindex,translossGrid,flipLUHgrids);

  readnc3dfield("c3per_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3per_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3per_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3per_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3per_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3per_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3per_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3per_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVC3PERGrid[clmlin * MAXOUTPIX + clmpix] - inCURRC3PERGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPC3PERGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readUNREPC4PERGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("c4per_to_secdn",yearindex,translossGrid,flipLUHgrids);

  readnc3dfield("c4per_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4per_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4per_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4per_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4per_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4per_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4per_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c4per_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVC4PERGrid[clmlin * MAXOUTPIX + clmpix] - inCURRC4PERGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPC4PERGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readUNREPC3NFXGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("c3nfx_to_secdn",yearindex,translossGrid,flipLUHgrids);

  readnc3dfield("c3nfx_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3nfx_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3nfx_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3nfx_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3nfx_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3nfx_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3nfx_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("c3nfx_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVC3NFXGrid[clmlin * MAXOUTPIX + clmpix] - inCURRC3NFXGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPC3NFXGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readUNREPPASTRGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("pastr_to_secdn",yearindex,translossGrid,flipLUHgrids);

  readnc3dfield("pastr_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("pastr_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("pastr_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("pastr_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("pastr_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("pastr_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("pastr_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("pastr_to_range",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVPASTRGrid[clmlin * MAXOUTPIX + clmpix] - inCURRPASTRGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPPASTRGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readUNREPRANGEGrids(int prevyear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = prevyear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhtransitionsdb); 
  
  readnc3dfield("range_to_secdn",yearindex,translossGrid,flipLUHgrids);

  readnc3dfield("range_to_urban",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("range_to_c3ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("range_to_c4ann",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("range_to_c3per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("range_to_c4per",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("range_to_c3nfx",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("range_to_pastr",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          translossGrid[clmlin * MAXOUTPIX + clmpix] = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
      }
  }

  readnc3dfield("range_to_secdf",yearindex,tempGrid,flipLUHgrids);
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          totaltransloss = translossGrid[clmlin * MAXOUTPIX + clmpix] + tempGrid[clmlin * MAXOUTPIX + clmpix];
          unreploss = totaltransloss - (inPREVRANGEGrid[clmlin * MAXOUTPIX + clmpix] - inCURRRANGEGrid[clmlin * MAXOUTPIX + clmpix]);
          if (unreploss < 0.0) {
              unreploss = 0.0;
          }
          inUNREPRANGEGrid[clmlin * MAXOUTPIX + clmpix] = unreploss;
      }
  }

  closencfile();
  
  return 0;
  
}


int readLUHcropmanagementGrids(int curryear) {

  int yearindex;
  long clmlin, clmpix;
  float totaltransloss;
  float unreploss;
  
  yearindex = curryear - firstyear;
  if (yearindex < 0) {
      yearindex = 0;
  }
  
  openncinputfile(luhmanagementdb); 
  
  readnc3dfield("fertl_c3ann",yearindex,inFERTC3ANNGrid,flipLUHgrids);
  readnc3dfield("fertl_c4ann",yearindex,inFERTC4ANNGrid,flipLUHgrids);
  readnc3dfield("fertl_c3per",yearindex,inFERTC3PERGrid,flipLUHgrids);
  readnc3dfield("fertl_c4per",yearindex,inFERTC4PERGrid,flipLUHgrids);
  readnc3dfield("fertl_c3nfx",yearindex,inFERTC3NFXGrid,flipLUHgrids);
  readnc3dfield("irrig_c3ann",yearindex,inIRRIGC3ANNGrid,flipLUHgrids);
  readnc3dfield("irrig_c4ann",yearindex,inIRRIGC4ANNGrid,flipLUHgrids);
  readnc3dfield("irrig_c3per",yearindex,inIRRIGC3PERGrid,flipLUHgrids);
  readnc3dfield("irrig_c4per",yearindex,inIRRIGC4PERGrid,flipLUHgrids);
  readnc3dfield("irrig_c3nfx",yearindex,inIRRIGC3NFXGrid,flipLUHgrids);

  closencfile();
  
  return 0;
  
}

int generateLUHcollectionGrids() {

  long clmlin, clmpix;
  
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
      
          inBASEFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] = inBASEPRIMFGrid[clmlin * MAXOUTPIX + clmpix] + inBASESECDFGrid[clmlin * MAXOUTPIX + clmpix];
          inBASENONFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] = inBASEPRIMNGrid[clmlin * MAXOUTPIX + clmpix] + inBASESECDNGrid[clmlin * MAXOUTPIX + clmpix];
          inBASECROPTOTALGrid[clmlin * MAXOUTPIX + clmpix] = inBASEC3ANNGrid[clmlin * MAXOUTPIX + clmpix] + inBASEC4ANNGrid[clmlin * MAXOUTPIX + clmpix] + inBASEC3PERGrid[clmlin * MAXOUTPIX + clmpix] + inBASEC4PERGrid[clmlin * MAXOUTPIX + clmpix] + inBASEC3NFXGrid[clmlin * MAXOUTPIX + clmpix];
           inBASEMISSINGGrid[clmlin * MAXOUTPIX + clmpix] = 1.0 - inBASEFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] - inBASENONFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] - inBASEPASTRGrid[clmlin * MAXOUTPIX + clmpix] - inBASERANGEGrid[clmlin * MAXOUTPIX + clmpix] - inBASECROPTOTALGrid[clmlin * MAXOUTPIX + clmpix];
          if (inBASEMISSINGGrid[clmlin * MAXOUTPIX + clmpix] < 0.0) {
              inBASEMISSINGGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          }
          if (inBASEMISSINGGrid[clmlin * MAXOUTPIX + clmpix] > 1.0) {
              inBASEMISSINGGrid[clmlin * MAXOUTPIX + clmpix] = 1.0;
          }
          inBASEOTHERGrid[clmlin * MAXOUTPIX + clmpix] = inBASEPRIMNGrid[clmlin * MAXOUTPIX + clmpix] + inBASESECDNGrid[clmlin * MAXOUTPIX + clmpix] + inBASERANGEGrid[clmlin * MAXOUTPIX + clmpix] + inBASEMISSINGGrid[clmlin * MAXOUTPIX + clmpix];
          inBASENATVEGGrid[clmlin * MAXOUTPIX + clmpix] = inBASEFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] + inBASEPASTRGrid[clmlin * MAXOUTPIX + clmpix] + inBASEOTHERGrid[clmlin * MAXOUTPIX + clmpix];

          inCURRFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] = inCURRPRIMFGrid[clmlin * MAXOUTPIX + clmpix] + inCURRSECDFGrid[clmlin * MAXOUTPIX + clmpix];
          inCURRNONFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] = inCURRPRIMNGrid[clmlin * MAXOUTPIX + clmpix] + inCURRSECDNGrid[clmlin * MAXOUTPIX + clmpix];
          inCURRCROPTOTALGrid[clmlin * MAXOUTPIX + clmpix] = inCURRC3ANNGrid[clmlin * MAXOUTPIX + clmpix] + inCURRC4ANNGrid[clmlin * MAXOUTPIX + clmpix] + inCURRC3PERGrid[clmlin * MAXOUTPIX + clmpix] + inCURRC4PERGrid[clmlin * MAXOUTPIX + clmpix] + inCURRC3NFXGrid[clmlin * MAXOUTPIX + clmpix];
          inCURRMISSINGGrid[clmlin * MAXOUTPIX + clmpix] = 1.0 - inCURRFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] - inCURRNONFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] - inCURRPASTRGrid[clmlin * MAXOUTPIX + clmpix] - inCURRRANGEGrid[clmlin * MAXOUTPIX + clmpix] - inCURRCROPTOTALGrid[clmlin * MAXOUTPIX + clmpix];
          if (inCURRMISSINGGrid[clmlin * MAXOUTPIX + clmpix] < 0.0) {
              inCURRMISSINGGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          }
          if (inCURRMISSINGGrid[clmlin * MAXOUTPIX + clmpix] > 1.0) {
              inCURRMISSINGGrid[clmlin * MAXOUTPIX + clmpix] = 1.0;
          }
          inCURROTHERGrid[clmlin * MAXOUTPIX + clmpix] = inCURRPRIMNGrid[clmlin * MAXOUTPIX + clmpix] + inCURRSECDNGrid[clmlin * MAXOUTPIX + clmpix] + inCURRRANGEGrid[clmlin * MAXOUTPIX + clmpix] + inCURRMISSINGGrid[clmlin * MAXOUTPIX + clmpix];
          inCURRNATVEGGrid[clmlin * MAXOUTPIX + clmpix] = inCURRFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] + inCURRPASTRGrid[clmlin * MAXOUTPIX + clmpix] + inCURROTHERGrid[clmlin * MAXOUTPIX + clmpix];

          inUNREPFORESTGrid[clmlin * MAXOUTPIX + clmpix] = inUNREPSECDFGrid[clmlin * MAXOUTPIX + clmpix] - inHARVESTSH1Grid[clmlin * MAXOUTPIX + clmpix] - inHARVESTSH2Grid[clmlin * MAXOUTPIX + clmpix];
          if (inUNREPFORESTGrid[clmlin * MAXOUTPIX + clmpix] < 0.0) {
              inUNREPFORESTGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          }
          if (inUNREPFORESTGrid[clmlin * MAXOUTPIX + clmpix] > 1.0) {
              inUNREPFORESTGrid[clmlin * MAXOUTPIX + clmpix] = 1.0;
          }
          inUNREPOTHERGrid[clmlin * MAXOUTPIX + clmpix] = inUNREPSECDNGrid[clmlin * MAXOUTPIX + clmpix] - inHARVESTSH3Grid[clmlin * MAXOUTPIX + clmpix];
          if (inUNREPOTHERGrid[clmlin * MAXOUTPIX + clmpix] < 0.0) {
              inUNREPOTHERGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          }
          if (inUNREPOTHERGrid[clmlin * MAXOUTPIX + clmpix] > 1.0) {
              inUNREPOTHERGrid[clmlin * MAXOUTPIX + clmpix] = 1.0;
          }
      }
  }
            
  return 0;
  
}


int generateclmPFTGrids() {

  long clmlin, clmpix;
  int pftid;
  float pctnatvegval, pctnatvegbase, forestunrepval, pastureunrepval, otherunrepval;
  float foresttotalbaseval, foresttotalfracval, foresttotalfracdelta, foresttotalcurrentval;
  float pasturebaseval, pasturefracval, pasturecurrentval, pasturefracdelta;
  float otherbaseval, otherfracval, othercurrentval, otherfracdelta;
  float currentpctforestpft, deltapctforestpft;
  float currentpctpasturepft, deltapctpasturepft;
  float currentpctotherpft, deltapctotherpft;
  float unrepforestfrac, unrepotherfrac;
  float newpctpft, unreppctpft, newpctpfttotal;
  
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          if (inLANDMASKGrid[clmlin * MAXOUTPIX + clmpix] == 1) {
              pctnatvegval = inCURRNATVEGGrid[clmlin * MAXOUTPIX + clmpix] * 100.0;
              if (pctnatvegval > 0.0) {
                  outPCTNATVEGGrid[clmlin * MAXOUTPIX + clmpix] = pctnatvegval;
                  pctnatvegbase = inBASENATVEGGrid[clmlin * MAXOUTPIX + clmpix] * 100.0;
                  forestunrepval = inUNREPFORESTGrid[clmlin * MAXOUTPIX + clmpix];
                  otherunrepval = inUNREPOTHERGrid[clmlin * MAXOUTPIX + clmpix];
                  if (pctnatvegbase > 0.0) {
                      foresttotalbaseval = inBASEFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegbase * 100.0;
                      foresttotalfracval = inCURRFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegval * 100.0;
                      foresttotalfracdelta = foresttotalfracval - foresttotalbaseval;
                      if (foresttotalfracdelta >= 0.0) {
                          foresttotalcurrentval = foresttotalbaseval;
                      }
                      else {
                          foresttotalcurrentval = foresttotalbaseval + foresttotalfracdelta;
                          foresttotalfracdelta = 0.0;
                      }
                      pasturebaseval = inBASEPASTRGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegbase * 100.0;
                      pasturefracval = inCURRPASTRGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegval * 100.0;
                      pasturecurrentval = 0.0;
                      pasturefracdelta = pasturefracval;
                      otherbaseval = inBASEOTHERGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegbase * 100.0;
                      otherfracval = inCURROTHERGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegval * 100.0;
                      otherfracdelta = otherfracval - otherbaseval;
                      if (otherfracdelta >= 0.0) {
                          othercurrentval = otherbaseval;
                      }
                      else {
                          othercurrentval = otherbaseval + otherfracdelta;
                          otherfracdelta = 0.0;
                      }
                  }
                  else {
                      foresttotalcurrentval = 0.0;
                      foresttotalfracdelta = inCURRFORESTTOTALGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegval * 100.0;
                      pasturecurrentval = 0.0;
                      pasturefracdelta = inCURRPASTRGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegval * 100.0;
                      othercurrentval = 0.0;
                      otherfracdelta = inCURROTHERGrid[clmlin * MAXOUTPIX + clmpix] / pctnatvegval * 100.0;
                  }
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      currentpctforestpft = foresttotalcurrentval * inCURRENTPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                      deltapctforestpft = foresttotalfracdelta * inFORESTPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                      unrepforestfrac = forestunrepval * (currentpctforestpft + deltapctforestpft) / 100.0;
                      currentpctpasturepft = pasturecurrentval * inCURRENTPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                      deltapctpasturepft = pasturefracdelta * inPASTUREPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                      currentpctotherpft = othercurrentval * inCURRENTPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                      deltapctotherpft = otherfracdelta * inOTHERPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                      newpctpft = currentpctforestpft + deltapctforestpft + currentpctpasturepft + deltapctpasturepft + currentpctotherpft + deltapctotherpft;
                      outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = newpctpft;
                      outUNREPPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = unrepforestfrac;
                  }
                  newpctpfttotal = 0.0;
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      newpctpfttotal = newpctpfttotal + outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                  }
                  if (newpctpfttotal > 0.0) {
                      for (pftid = 0; pftid < MAXPFT; pftid++) {
                          newpctpft = outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                          if (newpctpft > 0.0) {
                              newpctpft = newpctpft / newpctpfttotal * 100.0;
                              unreppctpft = unreppctpft / newpctpfttotal * 100.0;
                              if (unreppctpft > newpctpft) {
                                  unreppctpft = newpctpft;
                              }
                              outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = newpctpft;
                          }
                          else {
                              outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                          }
                      }
                  }
              }
              else {
                  outPCTNATVEGGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
                  outPCTPFTGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
                  outUNREPPFTGrid[0][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                      outUNREPPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  }
              }
          }
      }
  }

  return 0;
  
}


int generateclmCFTGrids() {

  long clmlin, clmpix;
  int cftid, rawcftid, rainfedcftid, irrigcftid;
  float pctcropval, c3annunrepval, c4annunrepval, c3perunrepval, c4perunrepval, c3nfxunrepval;
  float newpctrainfedcft, newpctirrigcft, newunreprainfedval, newunrepirrigval;
  float newpctcroptotal, newpctcft;

  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          if (inLANDMASKGrid[clmlin * MAXOUTPIX + clmpix] == 1) {
              pctcropval = inCURRCROPTOTALGrid[clmlin * MAXOUTPIX + clmpix] * 100.0;
              if (pctcropval > 0.0 && pctcropval <= 100.0) {
                  outPCTCROPGrid[clmlin * MAXOUTPIX + clmpix] = pctcropval;
                  c3annunrepval = inUNREPC3ANNGrid[clmlin * MAXOUTPIX + clmpix];
                  c4annunrepval = inUNREPC4ANNGrid[clmlin * MAXOUTPIX + clmpix];
                  c3perunrepval = inUNREPC3PERGrid[clmlin * MAXOUTPIX + clmpix];
                  c4perunrepval = inUNREPC4PERGrid[clmlin * MAXOUTPIX + clmpix];
                  c3nfxunrepval = inUNREPC3NFXGrid[clmlin * MAXOUTPIX + clmpix];
                  for (rawcftid = 0; rawcftid < MAXCFTRAW; rawcftid++) {
                      rainfedcftid = 2 * (rawcftid + 1);
                      irrigcftid = 2 * (rawcftid + 1) + 1;
                      newpctrainfedcft = inCURRC3ANNGrid[clmlin * MAXOUTPIX + clmpix] * (1.0 - inIRRIGC3ANNGrid[clmlin * MAXOUTPIX + clmpix]) * inC3ANNPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newpctirrigcft = inCURRC3ANNGrid[clmlin * MAXOUTPIX + clmpix] * (inIRRIGC3ANNGrid[clmlin * MAXOUTPIX + clmpix]) * inC3ANNPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newunreprainfedval = c3annunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c3annunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunreprainfedval;
                          outFERTNITROGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC3ANNGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] + newpctirrigcft;
                          outUNREPCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] + newunrepirrigval;
                          outFERTNITROGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC3ANNGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      newpctrainfedcft = inCURRC4ANNGrid[clmlin * MAXOUTPIX + clmpix] * (1.0 - inIRRIGC4ANNGrid[clmlin * MAXOUTPIX + clmpix]) * inC4ANNPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newpctirrigcft = inCURRC4ANNGrid[clmlin * MAXOUTPIX + clmpix] * (inIRRIGC4ANNGrid[clmlin * MAXOUTPIX + clmpix]) * inC4ANNPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newunreprainfedval = c4annunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c4annunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunreprainfedval;
                          outFERTNITROGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC4ANNGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunrepirrigval;
                          outFERTNITROGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC4ANNGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      newpctrainfedcft = inCURRC3PERGrid[clmlin * MAXOUTPIX + clmpix] * (1.0 - inIRRIGC3PERGrid[clmlin * MAXOUTPIX + clmpix]) * inC3PERPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newpctirrigcft = inCURRC3PERGrid[clmlin * MAXOUTPIX + clmpix] * (inIRRIGC3PERGrid[clmlin * MAXOUTPIX + clmpix]) * inC3PERPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newunreprainfedval = c3perunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c3perunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunreprainfedval;
                          outFERTNITROGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC3PERGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunrepirrigval;
                          outFERTNITROGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC3PERGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      newpctrainfedcft = inCURRC4PERGrid[clmlin * MAXOUTPIX + clmpix] * (1.0 - inIRRIGC4PERGrid[clmlin * MAXOUTPIX + clmpix]) * inC4PERPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newpctirrigcft = inCURRC4PERGrid[clmlin * MAXOUTPIX + clmpix] * (inIRRIGC4PERGrid[clmlin * MAXOUTPIX + clmpix]) * inC4PERPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newunreprainfedval = c4perunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c4perunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunreprainfedval;
                          outFERTNITROGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC4PERGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunrepirrigval;
                          outFERTNITROGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC4PERGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      newpctrainfedcft = inCURRC3NFXGrid[clmlin * MAXOUTPIX + clmpix] * (1.0 - inIRRIGC3NFXGrid[clmlin * MAXOUTPIX + clmpix]) * inC3NFXPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newpctirrigcft = inCURRC3NFXGrid[clmlin * MAXOUTPIX + clmpix] * (inIRRIGC3NFXGrid[clmlin * MAXOUTPIX + clmpix]) * inC3NFXPCTCFTGrid[rawcftid][clmlin * MAXOUTPIX + clmpix];
                      newunreprainfedval = c3nfxunrepval * newpctrainfedcft / 100.0;
                      newunrepirrigval = c3nfxunrepval * newpctirrigcft / 100.0;
                      if (newpctrainfedcft > 0.0) {
                          outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newpctrainfedcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunreprainfedval;
                          outFERTNITROGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC3NFXGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                      if (newpctirrigcft > 0.0) {
                          outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = outPCTCFTGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] + newpctirrigcft;
                          outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[rainfedcftid][clmlin * MAXOUTPIX + clmpix] + newunrepirrigval;
                          outFERTNITROGrid[irrigcftid][clmlin * MAXOUTPIX + clmpix] = inFERTC3NFXGrid[clmlin * MAXOUTPIX + clmpix] / 10.0;
                      }
                  }
                  newpctcroptotal = 0.0;
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      newpctcroptotal = newpctcroptotal + outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix];
                  }
                  if (newpctcroptotal > 0.0) {
                      for (cftid = 0; cftid < MAXCFT; cftid++) {
                          newpctcft = outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix];
                          if (newpctcft > 0.0) {
                              newpctcft = newpctcft / newpctcroptotal * 100.0;
                              outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix] = newpctcft;
                          }
                          else {
                              outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                          }
                      }
                  }
              }
              else {
                  outPCTCROPGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
                  outPCTCFTGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
                  outUNREPCFTGrid[0][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  for (cftid = 1; cftid < MAXCFT; cftid++) {
                      outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                      outUNREPCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  }
              }
          }
      }
  }

  return 0;
  
}


int generateclmwoodharvestGrids() {

  long clmlin, clmpix;
  int pftid;
  float TreePFTArea, TreeFrac, TreeScale, PFTArea;
  float newharvestvh1, newharvestvh2, newharvestsh1, newharvestsh2, newharvestsh3;
  float newbiohvh1, newbiohvh2, newbiohsh1, newbiohsh2, newbiohsh3;

  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          if (inLANDMASKGrid[clmlin * MAXOUTPIX + clmpix] == 1.0) {
              TreePFTArea = 0.0;
              TreeFrac = 0.0;
              PFTArea = inAREAGrid[clmlin * MAXOUTPIX + clmpix] * inLANDFRACGrid[clmlin * MAXOUTPIX + clmpix] * outPCTNATVEGGrid[clmlin * MAXOUTPIX + clmpix] / 100.0 * 1.0e6;              
              for (pftid = firsttreepft; pftid <= lasttreepft; pftid++) {
                  TreePFTArea = TreePFTArea + PFTArea * outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] / 100.0;
                  TreeFrac = TreeFrac + outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix] / 100.0;
              }
              TreeScale = 1.0;
              if (TreePFTArea > 1.0e6) {
                   newharvestvh1 = inHARVESTVH1Grid[clmlin * MAXOUTPIX + clmpix] * TreeScale;
                  if (newharvestvh1 < 0.0 || newharvestvh1 > 9.0e4) {
                      newharvestvh1 = 0.0;
                  }
                  if (newharvestvh1 > 0.98) {
                      newharvestvh1 = 0.98;
                  }
                  outHARVESTVH1Grid[clmlin * MAXOUTPIX + clmpix] = newharvestvh1;
                   newbiohvh1 = inBIOHVH1Grid[clmlin * MAXOUTPIX + clmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohvh1 < 0.0) {
                      newbiohvh1 = 0.0;
                  }
                  if (newbiohvh1 > 10000.0) {
                      newbiohvh1 = 10000.0;
                  }
                  outBIOHVH1Grid[clmlin * MAXOUTPIX + clmpix] = newbiohvh1;
                   newharvestvh2 = inHARVESTVH2Grid[clmlin * MAXOUTPIX + clmpix] * TreeScale;
                  if (newharvestvh2 < 0.0 || newharvestvh2 > 9.0e4) {
                      newharvestvh2 = 0.0;
                  }
                  if (newharvestvh2 > 0.98) {
                      newharvestvh2 = 0.98;
                  }
                  outHARVESTVH2Grid[clmlin * MAXOUTPIX + clmpix] = newharvestvh2;
                  newbiohvh2 = inBIOHVH2Grid[clmlin * MAXOUTPIX + clmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohvh2 < 0.0) {
                      newbiohvh2 = 0.0;
                  }
                  if (newbiohvh2 > 10000.0) {
                      newbiohvh2 = 10000.0;
                  }
                  outBIOHVH2Grid[clmlin * MAXOUTPIX + clmpix] = newbiohvh2;
                   newharvestsh1 = inHARVESTSH1Grid[clmlin * MAXOUTPIX + clmpix] * TreeScale;
                  if (newharvestsh1 < 0.0 || newharvestsh1 > 9.0e4) {
                      newharvestsh1 = 0.0;
                  }
                  if (newharvestsh1 > 0.98) {
                      newharvestsh1 = 0.98;
                  }
                  outHARVESTSH1Grid[clmlin * MAXOUTPIX + clmpix] = newharvestsh1;
                  newbiohsh1 = inBIOHSH1Grid[clmlin * MAXOUTPIX + clmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohsh1 < 0.0) {
                      newbiohsh1 = 0.0;
                  }
                  if (newbiohsh1 > 10000.0) {
                      newbiohsh1 = 10000.0;
                  }
                  outBIOHSH1Grid[clmlin * MAXOUTPIX + clmpix] = newbiohsh1;
                   newharvestsh2 = inHARVESTSH2Grid[clmlin * MAXOUTPIX + clmpix] * TreeScale;
                  if (newharvestsh2 < 0.0 || newharvestsh2 > 9.0e4) {
                      newharvestsh2 = 0.0;
                  }
                  if (newharvestsh2 > 0.98) {
                      newharvestsh2 = 0.98;
                  }
                  outHARVESTSH2Grid[clmlin * MAXOUTPIX + clmpix] = newharvestsh2;
                  newbiohsh2 = inBIOHSH2Grid[clmlin * MAXOUTPIX + clmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohsh2 < 0.0) {
                      newbiohsh2 = 0.0;
                  }
                  if (newbiohsh2 > 10000.0) {
                      newbiohsh2 = 10000.0;
                  }
                  outBIOHSH2Grid[clmlin * MAXOUTPIX + clmpix] = newbiohsh2;
                   newharvestsh3 = inHARVESTSH3Grid[clmlin * MAXOUTPIX + clmpix] * TreeScale;
                  if (newharvestsh3 < 0.0 || newharvestsh3 > 9.0e4) {
                      newharvestsh3 = 0.0;
                  }
                  if (newharvestsh3 > 0.98) {
                      newharvestsh3 = 0.98;
                  }
                  outHARVESTSH3Grid[clmlin * MAXOUTPIX + clmpix] = newharvestsh3;
                  newbiohsh3 = inBIOHSH3Grid[clmlin * MAXOUTPIX + clmpix] * 1000.0 / TreePFTArea * TreeScale;
                  if (newbiohsh3 < 0.0) {
                      newbiohsh3 = 0.0;
                  }
                  if (newbiohsh3 > 10000.0) {
                      newbiohsh3 = 10000.0;
                  }
                  outBIOHSH3Grid[clmlin * MAXOUTPIX + clmpix] = newbiohsh3;
              }
          }
      }
  }

  return 0;
  
}

int generatedblGrids() {

  double AllFrac, OtherFrac, AllPFTs, AllCFTs, tempdblPCT;
  long clmlin, clmpix;
  int pftid, cftid;
  
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          outAREAdblGrid[clmlin * MAXOUTPIX + clmpix] = inAREAGrid[clmlin * MAXOUTPIX + clmpix];
          if (inLANDMASKGrid[clmlin * MAXOUTPIX + clmpix] == 1.0) {
              outLANDFRACdblGrid[clmlin * MAXOUTPIX + clmpix] = inLANDFRACGrid[clmlin * MAXOUTPIX + clmpix];
              outPCTGLACIERdblGrid[clmlin * MAXOUTPIX + clmpix] = inPCTGLACIERGrid[clmlin * MAXOUTPIX + clmpix];
              outPCTLAKEdblGrid[clmlin * MAXOUTPIX + clmpix] = inPCTLAKEGrid[clmlin * MAXOUTPIX + clmpix];
              outPCTWETLANDdblGrid[clmlin * MAXOUTPIX + clmpix] = inPCTWETLANDGrid[clmlin * MAXOUTPIX + clmpix];
              outPCTURBANdblGrid[clmlin * MAXOUTPIX + clmpix] = inPCTURBANGrid[clmlin * MAXOUTPIX + clmpix];
              outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix] = outPCTCROPGrid[clmlin * MAXOUTPIX + clmpix];
              outPCTNATVEGdblGrid[clmlin * MAXOUTPIX + clmpix] = outPCTNATVEGGrid[clmlin * MAXOUTPIX + clmpix];	      
              OtherFrac = outPCTGLACIERdblGrid[clmlin * MAXOUTPIX + clmpix] + outPCTLAKEdblGrid[clmlin * MAXOUTPIX + clmpix] + outPCTWETLANDdblGrid[clmlin * MAXOUTPIX + clmpix] + outPCTURBANdblGrid[clmlin * MAXOUTPIX + clmpix];
              AllFrac = outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix] + outPCTNATVEGdblGrid[clmlin * MAXOUTPIX + clmpix];
              AllPFTs = 0.0;
              for (pftid = 0; pftid < MAXPFT; pftid++) {
                  AllPFTs = AllPFTs + outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
              }
              AllCFTs = 0.0;
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  AllCFTs = AllCFTs + outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix];
              }
              if (AllFrac == 0.0) {
                  outPCTNATVEGdblGrid[clmlin * MAXOUTPIX + clmpix] = 100.0;
                  outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
                  outPCTPFTdblGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
                  for (pftid = 1; pftid < MAXPFT; pftid++) {
                      outPCTPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  }
                  outPCTCFTdblGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
                  for (cftid = 1; cftid < MAXCFT; cftid++) {
                      outPCTCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outFERTNITROdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  }
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      outUNREPPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outUNREPCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                  }
                  outBIOHVH1dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
                  outBIOHVH2dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
                  outBIOHSH1dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
                  outBIOHSH2dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
                  outBIOHSH3dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              else {
                  outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix] = outPCTCROPGrid[clmlin * MAXOUTPIX + clmpix];
                  outPCTNATVEGdblGrid[clmlin * MAXOUTPIX + clmpix] = 100.0 - outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix];
                  if (AllPFTs == 0.0) {
                      outPCTPFTdblGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
                      for (pftid = 1; pftid < MAXPFT; pftid++) {
                          outPCTPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                      }
                  }
                  else {
                      for (pftid = 0; pftid < MAXPFT; pftid++) {
                          tempdblPCT = outPCTPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                          outPCTPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = tempdblPCT * 100.0 / AllPFTs;
                      }
                  }
                  if (AllCFTs == 0.0) {
                      outPCTCFTdblGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
                      for (cftid = 1; cftid < MAXCFT; cftid++) {
                          outPCTCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
                      }
                  }
                  else {
                      for (cftid = 0; cftid < MAXCFT; cftid++) {
                          tempdblPCT = outPCTCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix];
                          outPCTCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = tempdblPCT * 100.0 / AllCFTs;
                      }
                  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outFERTNITROdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = outFERTNITROGrid[cftid][clmlin * MAXOUTPIX + clmpix];
                  }
                  for (pftid = 0; pftid < MAXPFT; pftid++) {
                      outUNREPPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = outUNREPPFTGrid[pftid][clmlin * MAXOUTPIX + clmpix];
                  }
                  for (cftid = 0; cftid < MAXCFT; cftid++) {
                      outUNREPCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = outUNREPCFTGrid[cftid][clmlin * MAXOUTPIX + clmpix];
                  }
                  outBIOHVH1dblGrid[clmlin * MAXOUTPIX + clmpix] = outBIOHVH1Grid[clmlin * MAXOUTPIX + clmpix];
                  outBIOHVH2dblGrid[clmlin * MAXOUTPIX + clmpix] = outBIOHVH2Grid[clmlin * MAXOUTPIX + clmpix];
                  outBIOHSH1dblGrid[clmlin * MAXOUTPIX + clmpix] = outBIOHSH1Grid[clmlin * MAXOUTPIX + clmpix];
                  outBIOHSH2dblGrid[clmlin * MAXOUTPIX + clmpix] = outBIOHSH2Grid[clmlin * MAXOUTPIX + clmpix];
                  outBIOHSH3dblGrid[clmlin * MAXOUTPIX + clmpix] = outBIOHSH3Grid[clmlin * MAXOUTPIX + clmpix];
              }        
          } 
          else {
              outLANDFRACdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outPCTGLACIERdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outPCTLAKEdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outPCTWETLANDdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outPCTURBANdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outPCTNATVEGdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              for (pftid = 0; pftid < MAXPFT; pftid++) {
                  outPCTPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outPCTCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outFERTNITROdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              for (pftid = 0; pftid < MAXPFT; pftid++) {
                  outUNREPPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outUNREPCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              outBIOHVH1dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outBIOHVH2dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outBIOHSH1dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outBIOHSH2dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outBIOHSH3dblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
          }
      }
  }
  
  return 0;

}

int swapoceanGrids() {

  double scalelandunits;
  long clmlin, clmpix;
  int pftid, cftid;
  
  for (clmlin = 0; clmlin < MAXOUTLIN; clmlin++) {
      for (clmpix = 0; clmpix < MAXOUTPIX; clmpix++) {
          if (inLANDMASKGrid[clmlin * MAXOUTPIX + clmpix] == 0.0) {
              inLANDMASKGrid[clmlin * MAXOUTPIX + clmpix] = 1.0;
	      outLANDFRACdblGrid[clmlin * MAXOUTPIX + clmpix] = 1.0;
	      outPCTLAKEdblGrid[clmlin * MAXOUTPIX + clmpix] = 100.0;
	      outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix] = 0.0;
              outPCTPFTdblGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
              for (pftid = 1; pftid < MAXPFT; pftid++) {
                  outPCTPFTdblGrid[pftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              outPCTCFTdblGrid[0][clmlin * MAXOUTPIX + clmpix] = 100.0;
              for (cftid = 1; cftid < MAXCFT; cftid++) {
                  outPCTCFTdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }
              for (cftid = 0; cftid < MAXCFT; cftid++) {
                  outFERTNITROdblGrid[cftid][clmlin * MAXOUTPIX + clmpix] = 0.0;
              }	      
	  }
	  else {
	      scalelandunits = outLANDFRACdblGrid[clmlin * MAXOUTPIX + clmpix];
	      outLANDFRACdblGrid[clmlin * MAXOUTPIX + clmpix] = 1.0;
	      outPCTGLACIERdblGrid[clmlin * MAXOUTPIX + clmpix] = scalelandunits * outPCTGLACIERdblGrid[clmlin * MAXOUTPIX + clmpix];
	      outPCTLAKEdblGrid[clmlin * MAXOUTPIX + clmpix] = scalelandunits * outPCTLAKEdblGrid[clmlin * MAXOUTPIX + clmpix] + (1.0 - scalelandunits) * 100.0;
	      outPCTWETLANDdblGrid[clmlin * MAXOUTPIX + clmpix] = scalelandunits * outPCTWETLANDdblGrid[clmlin * MAXOUTPIX + clmpix];
	      outPCTURBANdblGrid[clmlin * MAXOUTPIX + clmpix] = scalelandunits * outPCTURBANdblGrid[clmlin * MAXOUTPIX + clmpix];
	      outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix] = scalelandunits * outPCTCROPdblGrid[clmlin * MAXOUTPIX + clmpix];
	      outPCTNATVEGdblGrid[clmlin * MAXOUTPIX + clmpix] = scalelandunits * outPCTNATVEGdblGrid[clmlin * MAXOUTPIX + clmpix];
	 }
      }
  }
	      
  return 0;
  
}


int writegrids(int currentyear) {

  char outncfilename[1024];
  long clmlin, clmpix;
  int pftid, cftid;
  char pftidstr[256];
  char cftidstr[256];

  sprintf(outncfilename,"%s/%s_%d.nc",outputdir,outputseries,currentyear);
  createncoutputfile(outncfilename);
  openncoutputfile(outncfilename);
  
  writenc1dintfield("natpft",innatpft);
  writenc1dintfield("cft",incft);
  writenc0dfield("EDGEN",&inEDGEN);
  writenc0dfield("EDGEE",&inEDGEE);
  writenc0dfield("EDGES",&inEDGES);
  writenc0dfield("EDGEW",&inEDGEW);
  writenc1dfield("LAT",inLAT);
  writenc2dfield("LATIXY",inLATIXY);
  writenc1dfield("LON",inLON);
  writenc2dfield("LONGXY",inLONGXY);
  writenc2dfield("LANDMASK",inLANDMASKGrid);
  writenc2ddblfield("LANDFRAC",outLANDFRACdblGrid);
  writenc2ddblfield("AREA",outAREAdblGrid);
  writenc2ddblfield("PCT_GLACIER",outPCTGLACIERdblGrid);
  writenc2ddblfield("PCT_LAKE",outPCTLAKEdblGrid);
  writenc2ddblfield("PCT_WETLAND",outPCTWETLANDdblGrid);
  writenc2ddblfield("PCT_URBAN",outPCTURBANdblGrid);
  writenc2ddblfield("PCT_NATVEG",outPCTNATVEGdblGrid);
  writenc2ddblfield("PCT_CROP",outPCTCROPdblGrid);
  
  for (pftid = 0; pftid < MAXPFT; pftid++) {
      writenc3ddblfield("PCT_NAT_PFT",pftid,outPCTPFTdblGrid[pftid]);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      writenc3ddblfield("PCT_CFT",cftid,outPCTCFTdblGrid[cftid]);
  }

  for (cftid = 0; cftid < MAXCFT; cftid++) {
      writenc3ddblfield("FERTNITRO_CFT",cftid,outFERTNITROdblGrid[cftid]);
  }

  for (pftid = 0; pftid < MAXPFT; pftid++) {
      writenc3ddblfield("UNREPRESENTED_PFT_LULCC",pftid,outUNREPPFTdblGrid[pftid]);
  }
  
  for (cftid = 0; cftid < MAXCFT; cftid++) {
      writenc3ddblfield("UNREPRESENTED_CFT_LULCC",cftid,outUNREPCFTdblGrid[cftid]);
  }

  writenc2ddblfield("HARVEST_VH1",outBIOHVH1dblGrid);
  writenc2ddblfield("HARVEST_VH2",outBIOHVH2dblGrid);
  writenc2ddblfield("HARVEST_SH1",outBIOHSH1dblGrid);
  writenc2ddblfield("HARVEST_SH2",outBIOHSH2dblGrid);
  writenc2ddblfield("HARVEST_SH3",outBIOHSH3dblGrid);

  closencfile();
  
  return 0;

}


main(long narg, char **argv) {

  int yearnumber;
    
  if(narg != 2){
        printf("Usage clm5landdatatool namelistfile\n");
        return 0;
  }
  
  readnamelist(argv[1]);
  setregionoptions();
  readpftparamfile();
  readcftrawparamfile();
  readcftparamfile();

  createallgrids();

  readclmcurrentGrids();
  readclmLUHforestGrids();
  readclmLUHpastureGrids();
  readclmLUHotherGrids();
  readclmLUHc3annGrids();
  readclmLUHc4annGrids();
  readclmLUHc3perGrids();
  readclmLUHc4perGrids();
  readclmLUHc3nfxGrids();

  readLUHbasestateGrids();
  
  for (yearnumber = startyear; yearnumber <= endyear; yearnumber++) {
  
      initializeGrids();
      
      readLUHcurrstateGrids(yearnumber);
      readLUHprevstateGrids(yearnumber-1);
  
      readLUHwoodharvestGrids(yearnumber-1);
  
/*      readUNREPSECDFGrids(yearnumber-1);
      readUNREPSECDNGrids(yearnumber-1);
      readUNREPPASTRGrids(yearnumber-1);
      readUNREPRANGEGrids(yearnumber-1);
      readUNREPC3ANNGrids(yearnumber-1);
      readUNREPC4ANNGrids(yearnumber-1);
      readUNREPC3PERGrids(yearnumber-1);
      readUNREPC4PERGrids(yearnumber-1);
      readUNREPC3NFXGrids(yearnumber-1); */

      readLUHcropmanagementGrids(yearnumber);

      generateLUHcollectionGrids();
      generateclmPFTGrids();
      generateclmCFTGrids();
      generateclmwoodharvestGrids();
      generatedblGrids();
      
      if (includeOcean == 0) {
          swapoceanGrids();
      }
      
      writegrids(yearnumber);

  }
  
  return 1;
  
}
