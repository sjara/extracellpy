/******************************************************************************
 Read data files from Neuralynx.

 Based on documentation (as of Dec 2010) from:
 http://www.neuralynx.com/static/software/NeuralynxDataFileFormats.pdf

 Santiago Jaramillo - 2011.01.13
******************************************************************************/

/*
Compile OBJ with:  gcc -c -O2 -Wall loadNeuralynx.c
Compile and link with: gcc -O2 -Wall loadNeuralynx.c -o load
*/

#include <stdio.h>
#include <stdlib.h>
#include "loadNeuralynx.h"

#define HEADERSIZE 16384

//#define NCS_INFOSIZE 20        // (64+32+32+32)/8
#define NCS_nSampPerRecord 512
#define NCS_RECORDSIZE 1044      // (64+32+32+32)/8 + (512*16)/8

#define NTT_nParamsPerRecord 8
#define NTT_nSampPerRecord  128  // 32*4
#define NTT_RECORDSIZE 304       // (64+32+32)/8 + (8*32)/8 + (128*16)/8

#define NEV_STRINGSIZE 128
#define NEV_EXTRABITSSIZE 8
#define NEV_RECORDSIZE 184       // (3*16 + 64 + 5*16 + 8*32)/8 + 128

typedef struct {
  TStype qwTimeStamp;
  unsigned int dwChannelNumber;
  unsigned int dwSampleFreq;
  unsigned int dwNumValidSamples;
  SAMPLEtype snSamples[NCS_nSampPerRecord];
} RecordNCS;

typedef struct {
  TStype qwTimeStamp;
  unsigned int dwScNumber;
  unsigned int dwCellNumber;
  unsigned int dnParams[NTT_nParamsPerRecord];
  SAMPLEtype snData[NTT_nSampPerRecord];
} RecordNTT;

typedef struct {
  short nstx;
  short npkt_id;
  short npkt_data_size;
  //short EMPTY;  // Because of struct padding this does not affect the size
  TStype qwTimeStamp;
  short nevent_id;
  unsigned short nttl;
  short ncrc;
  short ndummy1;
  short ndummy2;
  int dnExtra[NEV_EXTRABITSSIZE];
  char EventString[NEV_STRINGSIZE];
} RecordNEV;


DataNCS readNCS(char *fileName) {
  FILE *fid;
  char *header;
  TStype *timestamps;
  SAMPLEtype *samples;
  RecordNCS oneRecord;
  size_t count;
  unsigned long fileSize;
  unsigned long nRecords;
  int indrec, indsamp;
  DataNCS data;

  fid = fopen(fileName,"rb");
  if (!fid) { fprintf(stderr,"ERROR: File not found %s\n",fileName); }

  // -- Find number of records in file --
  fseek(fid,0,SEEK_END);
  fileSize = ftell(fid);
  fseek(fid,0,SEEK_SET);
  nRecords = (fileSize-HEADERSIZE)/NCS_RECORDSIZE;
  //printf("nRecords (inside read_ncs): %lu\n",nRecords);

  // -- Allocate space --
  header = (char *) malloc ( HEADERSIZE );
  timestamps = (TStype*) calloc (nRecords, sizeof(TStype));
  samples = (SAMPLEtype*) calloc (nRecords*NCS_nSampPerRecord, sizeof(SAMPLEtype));

  // -- Read data from file --
  count = fread(header, HEADERSIZE, 1, fid);
  for(indrec=0; indrec<nRecords; indrec++) {
    //count = fread(&oneRecord, sizeof(RecordNCS), 1, fid);
    //timestamps[indrec] = sizeof(RecordNCS);
    count = fread(&oneRecord, NCS_RECORDSIZE, 1, fid);
    timestamps[indrec] = oneRecord.qwTimeStamp;
    for(indsamp=0; indsamp<NCS_nSampPerRecord; indsamp++) {
      samples[indrec*NCS_nSampPerRecord+indsamp] = oneRecord.snSamples[indsamp];
    }
  }
  fclose(fid);

  // -- Update structure to be returned --
  data.header = header;
  data.timestamps = timestamps;
  data.nRecords = nRecords;
  data.samples = samples;
  data.nSamples = nRecords*NCS_nSampPerRecord;

  return data;
}

DataNTT readNTT(char *fileName) {
  FILE *fid;
  char *header;
  TStype *timestamps;
  SAMPLEtype *samples;
  RecordNTT oneRecord;
  size_t count;
  unsigned long fileSize;
  unsigned long nRecords;
  int indrec, indsamp;
  DataNTT data;

  fid = fopen(fileName,"rb");
  if (!fid) { fprintf(stderr,"ERROR: File not found %s\n",fileName); }

  // -- Find number of records in file --
  fseek(fid,0,SEEK_END);
  fileSize = ftell(fid);
  fseek(fid,0,SEEK_SET);
  nRecords = (fileSize-HEADERSIZE)/NTT_RECORDSIZE;
  //printf("nRecords (inside read_ntt): %lu/%lu = %lu\n",(fileSize-HEADERSIZE),NTT_RECORDSIZE,nRecords);

  // -- Allocate space --
  header = (char *) malloc ( HEADERSIZE );
  timestamps = (TStype*) calloc (nRecords, sizeof(TStype));
  samples = (SAMPLEtype*) calloc (nRecords*NTT_nSampPerRecord, sizeof(SAMPLEtype));

  // -- Read data from file --
  count = fread(header, HEADERSIZE, 1, fid);
  for(indrec=0; indrec<nRecords; indrec++) {
    count = fread(&oneRecord, sizeof(RecordNTT), 1, fid);
    timestamps[indrec] = oneRecord.qwTimeStamp;
    for(indsamp=0; indsamp<NTT_nSampPerRecord; indsamp++) {
      samples[indrec*NTT_nSampPerRecord+indsamp] = oneRecord.snData[indsamp];
    }
  }
  fclose(fid);

  // -- Update structure to be returned --
  data.header = header;
  data.timestamps = timestamps;
  data.nEvents = nRecords;
  data.samples = samples;
  data.nSamples = nRecords*NTT_nSampPerRecord;

  return data;
}


DataNEV readNEV(char *fileName) {
  FILE *fid;
  char *header;
  TStype *timestamps;
  short *eventID;
  unsigned short *valueTTL;
  char *eventString;
  RecordNEV oneRecord;
  size_t count;
  unsigned long fileSize;
  unsigned long nRecords;
  int indrec, indchar;
  DataNEV data;

  fid = fopen(fileName,"rb");
  if (!fid) { fprintf(stderr,"ERROR: File not found %s\n",fileName); }

  // -- Find number of records in file --
  fseek(fid,0,SEEK_END);
  fileSize = ftell(fid);
  fseek(fid,0,SEEK_SET);
  nRecords = (fileSize-HEADERSIZE)/NEV_RECORDSIZE;
  //printf("nRecords (inside read_ntt): %lu/%lu = %lu\n",(fileSize-HEADERSIZE),NEV_RECORDSIZE,nRecords);

  // -- Allocate space --
  header = (char *) malloc ( HEADERSIZE );
  timestamps = (TStype*) calloc (nRecords, sizeof(TStype));
  valueTTL = (unsigned short*) calloc (nRecords, sizeof(unsigned short));
  eventID = (short*) calloc (nRecords, sizeof(short));
  eventString = (char*) calloc (nRecords*NEV_STRINGSIZE, sizeof(char));

  /*
    The problem is with structure padding!
    Note that adding and extra EMPTY field does not affect the size.
    Find an alternative way of reading the data.
   */

  printf("========== SizeOf(RecordNEV) : %d ============\n\n",sizeof(RecordNEV));
  printf("========== SizeOf(dnExtra) : %d ============\n\n",sizeof(oneRecord.ncrc));
  // -- Read data from file --
  count = fread(header, HEADERSIZE, 1, fid);
  for(indrec=0; indrec<nRecords; indrec++) {
    count = fread(&oneRecord, sizeof(RecordNEV), 1, fid);
    /*
    count = fread(&oneRecord.nstx, sizeof(oneRecord.nstx), 1, fid);
    count = fread(&oneRecord.nstx, sizeof(oneRecord.nstx), 1, fid);
    count = fread(&oneRecord.npkt_id, sizeof(oneRecord.npkt_id), 1, fid);
    count = fread(&oneRecord.npkt_data_size, sizeof(oneRecord.npkt_data_size), 1, fid);
    count = fread(&oneRecord.qwTimeStamp, sizeof(oneRecord.qwTimeStamp), 1, fid);
    count = fread(&oneRecord.nevent_id, sizeof(oneRecord.nevent_id), 1, fid);
    count = fread(&oneRecord.nttl, sizeof(oneRecord.nttl), 1, fid);
    count = fread(&oneRecord.ncrc, sizeof(oneRecord.ncrc), 1, fid);
    count = fread(&oneRecord.ndummy1, sizeof(oneRecord.ndummy1), 1, fid);
    count = fread(&oneRecord.ndummy2, sizeof(oneRecord.ndummy2), 1, fid);
    count = fread(&oneRecord.dnExtra, sizeof(oneRecord.dnExtra), 1, fid);
    count = fread(&oneRecord.EventString, sizeof(oneRecord.EventString), 1, fid);
    */
    timestamps[indrec] = oneRecord.qwTimeStamp;
    valueTTL[indrec] = oneRecord.nttl;
    eventID[indrec] = oneRecord.nevent_id;
    for(indchar=0; indchar<NEV_STRINGSIZE; indchar++) {
      eventString[indrec*NEV_STRINGSIZE+indchar] = oneRecord.EventString[indchar];
    }
  }
  fclose(fid);

  // -- Update structure to be returned --
  data.header = header;
  data.timestamps = timestamps;
  data.nEvents = nRecords;
  data.valueTTL = valueTTL;
  data.eventID = eventID;
  data.eventString = eventString;
  data.stringSize = NEV_STRINGSIZE;

  return data;
}



int main() {
  char fileName[] = "Events.nev";
  DataNEV data;
  int ind;

  data = readNEV(fileName);
  printf("%s\n",data.header);
  for(ind=0;ind<4;ind++) { printf("%Lu, ",data.timestamps[ind]); }
  printf("... , %Lu\n",data.timestamps[data.nEvents-1]);
  for(ind=0;ind<4;ind++) { printf("%hd, ",data.valueTTL[ind]); }
  printf("... , %hd\n",data.valueTTL[data.nEvents-1]);
  printf("nEvents (main): %lu\n",data.nEvents);

  /*
  char fileName[] = "TT4.ntt";
  DataNTT data;
  int ind;

  data = readNTT(fileName);
  printf("%s\n",data.header);
  for(ind=0;ind<4;ind++) { printf("%Lu, ",data.timestamps[ind]); }
  printf("... , %Lu\n",data.timestamps[data.nEvents-1]);
  for(ind=0;ind<132;ind++) { printf("%hd, ",data.samples[ind]); }
  printf("... , %hd\n",data.samples[data.nSamples-1]);
  printf("nEvents (main): %lu\n",data.nEvents);
  printf("nSamples (main): %lu\n",data.nSamples);
  */

  /*
  char fileName[] = "CSC1.ncs";
  DataNCS data;
  int ind;

  data = readNCS(fileName);
  printf("%s\n",data.header);
  for(ind=0;ind<4;ind++) { printf("%Lu, ",data.timestamps[ind]); }
  printf("... , %Lu\n",data.timestamps[data.nRecords-1]);
  for(ind=0;ind<4;ind++) { printf("%hd, ",data.samples[ind]); }
  printf("... , %hd\n",data.samples[data.nSamples-1]);
  printf("nRecords (main): %lu\n",data.nRecords);
  printf("nSamples (main): %lu\n",data.nSamples);
  */
  return 0;
}

