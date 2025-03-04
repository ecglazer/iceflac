/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Role:
**    Converts FLAC's binary outputs to VTK format
**
** Author:
**
**    Eunseo Choi
**    Lamont-Doherty Earth Observatory (echoi@ldeo.columbia.edu)
**    P.O. Box 1000
**    61 Rt. 9W
**    Palisades, NY 10964, USA,
**
**    Eh Tan
**    Institute for Geophysics (tan2@mail.utexas.edu)
**    University of Texas, Austin
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <limits.h>
#ifndef PATH_MAX
	#define PATH_MAX 1024
#endif

#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2

//#define DEBUG


void ConvertTimeStep( int rank, unsigned int dumpIteration, double time, unsigned int gnode, unsigned int gelement, unsigned int nodeNumber[2], unsigned int elementNumber[2] );
void writeNodeVectorFloat(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile);
void writeNodeScalarDouble(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile);
void writeNodeScalarFloat(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile);
void writeElementScalarDouble(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile);
void writeElementScalarFloat(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile);
void writeElementScalarInt(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile);

char		path[PATH_MAX];
FILE*		coordIn;
FILE*		velIn;
FILE*		phaseIn;
FILE*		apsIn;
FILE*		eIIn;
FILE*		eIIIn;
FILE*		strainRateIn;
FILE*		stressIn;
FILE*		sxxIn;
FILE*		szzIn;
FILE*		sxzIn;
FILE*		pressureIn;
FILE*		tempIn;
FILE*		densIn;
FILE*		viscIn;
FILE*		srcIn;
FILE*		dissIn;

unsigned int	nodeNumber[2];
unsigned int	elementNumber[2];
unsigned int	globalNodeNumber;
unsigned int	globalElementNumber;
unsigned int	globalBoundaryNumber;
int 		doTemp = 1;
int 		doForce = 1;
int 		doAps = 1;
int 		doHPr = 1;
int 		doVisc = 1;
const unsigned	numVelocityVectorComponent = 2; /* 2 components in stress vector */
const unsigned	numStressVectorComponent = 6; /* 6 components in stress vector */
const unsigned	numStressComponentsPerElement = 6; /* 6 averaged components per element */

int main( int argc, char* argv[]) 
{
    char		tmpBuf[PATH_MAX];
    FILE*		nxnzIn;
    FILE*		timeStepIn;
    unsigned int	rank;
    unsigned int	dumpIteration, junk;
    float		time;
    unsigned int	nrec=0, step=0, stepMin=0, stepMax=0;
	
    if( argc < 2 || argc > 4 ) {
		fprintf(stderr,"usage: flac2vtk path [step_min [step_max]]\n\n"
                        "Processing flac data output to VTK format.\n"
                        "If step_max is not given, processing to latest steps\n"
                        "If both step_min and step_max are not given, processing all steps\n");
		exit(1);
    }

    /*
     * Set the default input/output path and range of time steps to process.
     */
	sprintf( path, "%s", argv[1] );

        /* Read the max steps */
	sprintf( tmpBuf, "%s/_contents.0", path );
	if( (timeStepIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	while(!feof(timeStepIn)) {
            fscanf( timeStepIn, "%d %d %g", &nrec, &step, &time );
        }
	rewind( timeStepIn );

        /* Set the range of steps to process. */
        if( argc == 3 ) {
            /* Read steps from arguments */
            stepMin = (unsigned int) atoi(argv[2]);
            stepMax = nrec;
        } else if( argc == 4 ) {
            /* Read steps from arguments */
            stepMin = (unsigned int) atoi(argv[2]);
            stepMax = (unsigned int) atoi(argv[3]);
            if(stepMax > nrec) stepMax = nrec;
        } else {
            /* all steps */
            stepMin = 1;
            stepMax = nrec;
        }

        /* Get mesh size */
	sprintf( tmpBuf, "%s/nxnz.0", path );
	if( (nxnzIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}

	fscanf( nxnzIn, "%d %d", &elementNumber[0], &elementNumber[1] );
	fclose( nxnzIn );

	nodeNumber[0] = elementNumber[0] + 1;
	nodeNumber[1] = elementNumber[1] + 1;
	globalNodeNumber = nodeNumber[0]*nodeNumber[1];
	globalElementNumber = elementNumber[0]*elementNumber[1];

	/* Print out some information */
	fprintf(stderr, "Time step range:  %u <-> %u, # of steps=%d\n", stepMin, stepMax, stepMax-stepMin+1 );
	fprintf(stderr, "nnode:%u nelem:%u nx:%u nz:%u ex:%u ez:%u\n",
			globalNodeNumber, globalElementNumber, nodeNumber[0], nodeNumber[1],
			elementNumber[0], elementNumber[1]);

	/* Start processing for each rank */
	rank = 0;
			
	sprintf( tmpBuf, "%s/mesh.%u", path, rank );
	if( (coordIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/vel.%u", path, rank );
	if( (velIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/temperature.%u", path, rank );
	if( (tempIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/srII.%u", path, rank );
	if( (strainRateIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/eI.%u", path, rank );
	if( (eIIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/eII.%u", path, rank );
	if( (eIIIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/phase.%u", path, rank );
	if( (phaseIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/aps.%u", path, rank );
	if( (apsIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/sII.%u", path, rank );
	if( (stressIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/sxx.%u", path, rank );
	if( (sxxIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/szz.%u", path, rank );
	if( (szzIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/sxz.%u", path, rank );
	if( (sxzIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/pres.%u", path, rank );
	if( ( pressureIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		doHPr = 0;
	}
	sprintf( tmpBuf, "%s/visc.%u", path, rank );
	if( (viscIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/density.%u", path, rank );
	if( (densIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/src.%u", path, rank );
	if( (srcIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	sprintf( tmpBuf, "%s/diss.%u", path, rank );
	if( (dissIn = fopen( tmpBuf, "r" )) == NULL ) {
		fprintf(stderr, "\"%s\" not found\n", tmpBuf );
		exit(1);
	}
	/*
	 * Read in loop information and write VTK files for wanted time steps.
	 */
	for( dumpIteration = 1; dumpIteration < stepMin; dumpIteration++ ) {
            fscanf( timeStepIn, "%d %d %g", &junk, &step, &time );
        }
	for( dumpIteration = stepMin-1; dumpIteration < stepMax; dumpIteration++ ) {
            fscanf( timeStepIn, "%d %d %g", &junk, &step, &time );
            fprintf(stderr,"trying %d-th conversion: model time=%.2e Myr)\n", dumpIteration+1, time);
            ConvertTimeStep( rank, dumpIteration, time, globalNodeNumber, globalElementNumber,
                             nodeNumber, elementNumber );
	}
		
	/*
	 * Close the input files 
	 */
        fclose( timeStepIn );
	fclose( coordIn );
	fclose( velIn );
	fclose( apsIn );
	fclose( phaseIn );
	fclose( eIIn );
	fclose( eIIIn );
	fclose( strainRateIn );
	fclose( pressureIn );
	fclose( stressIn );
	fclose( sxxIn );
	fclose( szzIn );
	fclose( sxzIn );
	fclose( tempIn );
	fclose( viscIn );
	fclose( densIn );
	fclose( srcIn );
	fclose( dissIn );
	
    return 0;
}

void ConvertTimeStep( 
					 int rank, 
					 unsigned int dumpIteration, 
					 double time,
					 unsigned int globalNodeNumber,
					 unsigned int globalElementNumber,
					 unsigned int nodeNumber[2], 
					 unsigned int elementNumber[2]
					  ) 
{
    char		tmpBuf[PATH_MAX];
    FILE*		vtkOut;
    FILE*		topoOut;
    unsigned int	node_gI,node_gJ,dimI;
	
    /*
     * Open the output file 
     */
    sprintf( tmpBuf, "%s/flac.%06u.vts", path, dumpIteration+1 );
    fprintf(stderr, "%s\n", tmpBuf);
    if( (vtkOut = fopen( tmpBuf, "w+" )) == NULL ) {
		fprintf(stderr, "Cannot open \"%s\" for writing\n", tmpBuf );
		exit(1);
    }
    sprintf( tmpBuf, "%s/topo.%06u.dat", path, dumpIteration+1 );
    if( (topoOut = fopen( tmpBuf, "w+" )) == NULL ) {
		fprintf(stderr, "Cannot open \"%s\" for writing\n", tmpBuf );
		exit(1);
    }
	
    /*
     * Write out simulation information 
     */
    fprintf( vtkOut, "<?xml version=\"1.0\"?>\n" );
    fprintf( vtkOut, "<VTKFile type=\"StructuredGrid\"  version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
    fprintf( vtkOut, "  <StructuredGrid WholeExtent=\"%d %d %d %d 0 0\">\n",
			 0,elementNumber[0],0,elementNumber[1]);
    fprintf( vtkOut, "    <Piece Extent=\"%d %d %d %d 0 0\">\n",
			 0,elementNumber[0],0,elementNumber[1]);
	
    /*
     *
     *  			Start the ---node--- section 
     *
     */
    fprintf( vtkOut, "      <PointData Vectors=\"Velocity\">\n");

    /*
     * Write out the velocity information 
     */
	writeNodeVectorFloat(dumpIteration, "Velocity", velIn, vtkOut );

    /*
     * Write out the temperature information
     */
      	writeNodeScalarFloat(dumpIteration, "Temperature", tempIn, vtkOut );

    fprintf( vtkOut, "      </PointData>\n");

    /* 
     *
     *  			Start the ---element--- section 
     *
     */
    fprintf( vtkOut, "      <CellData Scalars=\"Phase\">\n");

    /*
     * Write out the phase information 
     */
        writeElementScalarInt(dumpIteration, "Phase", phaseIn, vtkOut);

    /*
     * Write out the plastic strain information 
     */
      	writeElementScalarFloat(dumpIteration, "Plastic strain", apsIn, vtkOut);

	/*
     * Write out the total strain information 
     */	
      	writeElementScalarFloat(dumpIteration, "eI", eIIn, vtkOut);
      	writeElementScalarFloat(dumpIteration, "eII", eIIIn, vtkOut);

    /*
     * Write out the strain rate information 
     */
      	writeElementScalarFloat(dumpIteration, "Strain rate", strainRateIn, vtkOut);
	

    /*
     * Write out the stress (second inv. of dev. stress) information 
     */
      	writeElementScalarFloat(dumpIteration, "Stress", stressIn, vtkOut);

    /*
     * Write out the xx,yy,xy components of stress
     */
      	writeElementScalarFloat(dumpIteration, "Sxx", sxxIn, vtkOut);
      	writeElementScalarFloat(dumpIteration, "Szz", szzIn, vtkOut);
      	writeElementScalarFloat(dumpIteration, "Sxz", sxzIn, vtkOut);

    /*
     * Write out the pressure information 
     */
      	writeElementScalarFloat(dumpIteration, "Pressure", pressureIn, vtkOut);

    /*
     * Write out the thermal information 
     */
      	writeElementScalarFloat(dumpIteration, "Viscosity", viscIn, vtkOut);
        writeElementScalarFloat(dumpIteration, "Density", densIn, vtkOut );
      	writeElementScalarFloat(dumpIteration, "HeatSource", srcIn, vtkOut );
        writeElementScalarFloat(dumpIteration, "Dissipation", dissIn, vtkOut );

    fprintf( vtkOut, "      </CellData>\n");

    /* 
     *
     *  			Write out coordinates
     *
     */
    fprintf( vtkOut, "      <Points>\n");
    fprintf( vtkOut, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    if (fseek( coordIn, dumpIteration * globalNodeNumber * sizeof(float) * 2, SEEK_SET )!=0) {
		fprintf(stderr, "Cannot find read required portion of GeoFlac coordinates output file:  dump iteration=%d, node count=%d\n", dumpIteration, globalNodeNumber );
		exit(1);
    }
	{
		float		coord[globalNodeNumber][2];
		for( dimI=0; dimI<2; dimI++ )
			for( node_gI = 0; node_gI < globalNodeNumber; node_gI++ ) {
				if (fread( &coord[node_gI][dimI], sizeof(float), 1, coordIn )==0)  {
					if (feof(coordIn)) {
						fprintf(stderr, "Error (reached EOF prematurely) while reading GeoFLAC coordinates output file:  dump iteration=%d, node=%d/%d\n", dumpIteration, node_gI, globalNodeNumber );
						exit(1);
					}
					else if(ferror(coordIn)) {
						fprintf(stderr, "Error while reading GeoFLAC coordinates output file:  dump iteration=%d, node=%d/%d\n", dumpIteration, node_gI, globalNodeNumber );
						exit(1);
					}
				}
			}
		for( node_gJ = 0; node_gJ < nodeNumber[1]; node_gJ++ ) 
			for( node_gI = 0; node_gI < nodeNumber[0]; node_gI++ ) {
				unsigned int id=node_gJ + node_gI*nodeNumber[1];
				fprintf( vtkOut, "%g %g 0\n", coord[id][0], coord[id][1] );
				if( node_gJ == 0 )
					fprintf( topoOut, "%g %g\n", coord[id][0], coord[id][1] );
			}
	}
	fprintf( vtkOut, "        </DataArray>\n");
	fprintf( vtkOut, "      </Points>\n");
    fprintf( vtkOut, "    </Piece>\n");
    fprintf( vtkOut, "  </StructuredGrid>\n");
    fprintf( vtkOut, "</VTKFile>\n");
	
    /*
     * Close the output files 
     */
    fclose( vtkOut );
    fclose( topoOut );
}
 
void writeNodeVectorFloat(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile) {
	unsigned int node_gI,node_gJ;
	float		nodeArray[globalNodeNumber][2];
	
	fprintf( vtkFile, "        <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", varName);
	if( fseek( file, dumpIteration * globalNodeNumber * numVelocityVectorComponent * sizeof(float), SEEK_SET )!=0 ) {
		fprintf(stderr, "Cannot find read required portion of GeoFLAC %s output file:  dump iteration=%d, node count=%d\n", varName, dumpIteration, globalNodeNumber );
		exit(1);
	}
	for( node_gI = 0; node_gI < globalNodeNumber; node_gI++ ) {
		if( fread( &nodeArray[node_gI][0], sizeof(float), 1, file )==0 ) {
			if( feof(file) ) {
				fprintf(stderr, "In %s: Error (reached EOF prematurely) while reading %s output file:  dump iteration=%d, node=%d/%d\n",
						__func__, varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			} 
			else if( ferror(file) ) {
				fprintf(stderr, "Error while reading GeoFLAC %s output file:  dump iteration=%d, node=%d/%d\n",
						varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			}
		}
	}
	for( node_gI = 0; node_gI < globalNodeNumber; node_gI++ ) {
		if( fread( &nodeArray[node_gI][1], sizeof(float), 1, file )==0 ) {
			if( feof(file) ) {
				fprintf(stderr, "In %s: Error (reached EOF prematurely) while reading %s output file:  dump iteration=%d, node=%d/%d\n",
						__func__, varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			} 
			else if( ferror(file) ) {
				fprintf(stderr, "Error while reading GeoFLAC %s output file:  dump iteration=%d, node=%d/%d\n",
						varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			}
		}
	}
	for( node_gJ = 0; node_gJ < nodeNumber[1]; node_gJ++ )
		for( node_gI = 0; node_gI < nodeNumber[0]; node_gI++ ) {
			unsigned int id=node_gJ + node_gI*nodeNumber[1];
			fprintf( vtkFile, "%g %g 0\n", nodeArray[id][0], nodeArray[id][1] );
		}
	fprintf( vtkFile, "        </DataArray>\n");
}

void writeNodeScalarDouble(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile) {
	unsigned int node_gI,node_gJ;
	double		nodeScalar[globalNodeNumber];

	fprintf( vtkFile, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", varName);
	if (fseek( file, dumpIteration * globalNodeNumber * sizeof(double), SEEK_SET )!=0) {
		fprintf(stderr, "Cannot find read required portion of GeoFLAC %s output file:  dump iteration=%d, node count=%d\n", varName, dumpIteration, globalNodeNumber );
		exit(1);
	}
	for( node_gI = 0; node_gI < globalNodeNumber; node_gI++ ) {
	    if (fread( &nodeScalar[node_gI], sizeof(double), 1, file )==0)  {
			if (feof(file)) {
				fprintf(stderr, "In %s: Error (reached EOF prematurely) while reading GeoFLAC %s output file:  dump iteration=%d, node=%d/%d\n", 
						__func__, varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			} 
			else if(ferror(file)) {
				fprintf(stderr, "In %s: Error while reading GeoFLAC %s output file:  dump iteration=%d, node=%d/%d\n", 
						__func__, varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			}
	    }
	}
	for( node_gJ = 0; node_gJ < nodeNumber[1]; node_gJ++ ) 
		for( node_gI = 0; node_gI < nodeNumber[0]; node_gI++ ) {
			unsigned int id=node_gJ + node_gI*nodeNumber[1];
			fprintf( vtkFile, "%g\n", nodeScalar[id]);
		}
	fprintf( vtkFile, "        </DataArray>\n");
}

void writeNodeScalarFloat(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile) {
	unsigned int node_gI,node_gJ;
	float		nodeScalar[globalNodeNumber];

	fprintf( vtkFile, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", varName);
	if (fseek( file, dumpIteration * globalNodeNumber * sizeof(float), SEEK_SET )!=0) {
		fprintf(stderr, "Cannot find read required portion of GeoFLAC %s output file:  dump iteration=%d, node count=%d\n", varName, dumpIteration, globalNodeNumber );
		exit(1);
	}
	for( node_gI = 0; node_gI < globalNodeNumber; node_gI++ ) {
	    if (fread( &nodeScalar[node_gI], sizeof(float), 1, file )==0)  {
			if (feof(file)) {
				fprintf(stderr, "In %s: Error (reached EOF prematurely) while reading GeoFLAC %s output file:  dump iteration=%d, node=%d/%d\n", 
						__func__, varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			} 
			else if(ferror(file)) {
				fprintf(stderr, "In %s: Error while reading GeoFLAC %s output file:  dump iteration=%d, node=%d/%d\n", 
						__func__, varName, dumpIteration, node_gI, globalNodeNumber );
				exit(1);
			}
	    }
	}
	for( node_gJ = 0; node_gJ < nodeNumber[1]; node_gJ++ ) 
		for( node_gI = 0; node_gI < nodeNumber[0]; node_gI++ ) {
			unsigned int id=node_gJ + node_gI*nodeNumber[1];
			fprintf( vtkFile, "%g\n", nodeScalar[id]);
		}
	fprintf( vtkFile, "        </DataArray>\n");
}

void writeElementScalarDouble(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile) {
	unsigned int element_gI,element_gJ;
	double		elementScalar[globalElementNumber];

	fprintf( vtkFile, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", varName);
	if (fseek( file, dumpIteration * globalElementNumber * sizeof(double), SEEK_SET )!=0) {
		fprintf(stderr, "Cannot find read required portion of GeoFLAC %s output file:  dump iteration=%d, element count=%d\n", varName, dumpIteration, globalElementNumber );
		exit(1);
	}
	for( element_gI = 0; element_gI < globalElementNumber; element_gI++ ) {
	    if (fread( &elementScalar[element_gI], sizeof(double), 1, file )==0)  {
			if (feof(file)) {
				fprintf(stderr, "In %s: Error (reached EOF prematurely) while reading GeoFLAC %s output file:  dump iteration=%d, element=%d/%d\n",
						__func__, varName, dumpIteration, element_gI, globalElementNumber );
				exit(1);
			} 
			else if(ferror(file)) {
				fprintf(stderr, "In %s: Error while reading GeoFLAC %s output file:  dump iteration=%d, element=%d/%d\n", 
						__func__, varName, dumpIteration, element_gI, globalElementNumber );
				exit(1);
			}
	    }
	}
	for( element_gJ = 0; element_gJ < elementNumber[1]; element_gJ++ ) 
		for( element_gI = 0; element_gI < elementNumber[0]; element_gI++ ) {
			unsigned int id=element_gJ + element_gI*elementNumber[1];
	    fprintf( vtkFile, "%g\n", ((fabs(elementScalar[id])<1.0e-32)?0.0:elementScalar[id]));
	}
	fprintf( vtkFile, "        </DataArray>\n");
 }

void writeElementScalarFloat(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile) {
	unsigned int element_gI,element_gJ;
	float		elementScalar[globalElementNumber];

	fprintf( vtkFile, "        <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", varName);
	if (fseek( file, dumpIteration * globalElementNumber * sizeof(float), SEEK_SET )!=0) {
		fprintf(stderr, "Cannot find read required portion of GeoFLAC %s output file:  dump iteration=%d, element count=%d\n", varName, dumpIteration, globalElementNumber );
		exit(1);
	}
	for( element_gI = 0; element_gI < globalElementNumber; element_gI++ ) {
	    if (fread( &elementScalar[element_gI], sizeof(float), 1, file )==0)  {
			if (feof(file)) {
				fprintf(stderr, "In %s: Error (reached EOF prematurely) while reading GeoFLAC %s output file:  dump iteration=%d, element=%d/%d\n", 
						__func__, varName, dumpIteration, element_gI, globalElementNumber );
				exit(1);
			} 
			else if(ferror(file)) {
				fprintf(stderr, "In %s: Error while reading GeoFLAC %s output file:  dump iteration=%d, element=%d/%d\n", 
						__func__, varName, dumpIteration, element_gI, globalElementNumber );
				exit(1);
			}
	    }
	}
	for( element_gJ = 0; element_gJ < elementNumber[1]; element_gJ++ ) 
		for( element_gI = 0; element_gI < elementNumber[0]; element_gI++ ) {
			unsigned int id=element_gJ + element_gI*elementNumber[1];
	    fprintf( vtkFile, "%g\n", ((fabs(elementScalar[id])<1.0e-32)?0.0:elementScalar[id]));
	}
	fprintf( vtkFile, "        </DataArray>\n");
 }


void writeElementScalarInt(unsigned int dumpIteration, char* varName, FILE* file, FILE* vtkFile) {
	unsigned int element_gI,element_gJ;
	int		elementScalar[globalElementNumber];

	fprintf( vtkFile, "        <DataArray type=\"Int32\" Name=\"%s\" format=\"ascii\">\n", varName);
	if (fseek( file, dumpIteration * globalElementNumber * sizeof(int), SEEK_SET )!=0) {
		fprintf(stderr, "Cannot find read required portion of GeoFLAC %s output file:  dump iteration=%d, element count=%d\n", varName, dumpIteration, globalElementNumber );
		exit(1);
	}
	for( element_gI = 0; element_gI < globalElementNumber; element_gI++ ) {
	    if (fread( &elementScalar[element_gI], sizeof(int), 1, file )==0)  {
			if (feof(file)) {
				fprintf(stderr, "In %s: Error (reached EOF prematurely) while reading GeoFLAC %s output file:  dump iteration=%d, element=%d/%d\n", 
						__func__, varName, dumpIteration, element_gI, globalElementNumber );
				exit(1);
			}
			else if(ferror(file)) {
				fprintf(stderr, "In %s: Error while reading GeoFLAC %s output file:  dump iteration=%d, element=%d/%d\n", 
						__func__, varName, dumpIteration, element_gI, globalElementNumber );
				exit(1);
			}
	    }
	}
	for( element_gJ = 0; element_gJ < elementNumber[1]; element_gJ++ )
		for( element_gI = 0; element_gI < elementNumber[0]; element_gI++ ) {
			unsigned int id=element_gJ + element_gI*elementNumber[1];
	    fprintf( vtkFile, "%d\n", elementScalar[id]);
	}
	fprintf( vtkFile, "        </DataArray>\n");
 }
