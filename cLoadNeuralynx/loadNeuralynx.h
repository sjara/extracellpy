/******************************************************************************
 Definitions for reading data files from Neuralynx.

Santiago Jaramillo - 2011-01-13
******************************************************************************/

typedef unsigned long long TStype;
typedef unsigned short SAMPLEtype;

/* -- Continuous data -- */ 
typedef struct {
  char *header;             // 16kB text
  TStype *timestamps;       // For each record, in usec
  unsigned long nRecords;
  SAMPLEtype *samples;      // Continuous waveform
  unsigned long nSamples;
} DataNCS;

DataNCS readNCS(char *fileName);


/* -- Tetrode data -- */
typedef struct {
  char *header;             // 16kB text
  TStype *timestamps;       // For each record/spike, in usec
  unsigned long nEvents;
  SAMPLEtype *samples;      // For each chan, each record/spike
  unsigned long nSamples;
} DataNTT;

DataNTT readNTT(char *fileName);


/* -- Events data -- */
typedef struct {
  char *header;             // 16kB text
  TStype *timestamps;       // For each event, in usec
  unsigned long nEvents;
  unsigned short *valueTTL; // 
  short *eventID;           //
  char *eventString;        // 128 characters for each event
  int stringSize;
} DataNEV;

DataNEV readNEV(char *fileName);
